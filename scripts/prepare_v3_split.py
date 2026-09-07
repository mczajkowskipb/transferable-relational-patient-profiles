#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
import numpy as np

from relational_patient_profiles.artifact import (
    RPPRelation,RPPPrototype,RelationalPatientProfileArtifact,execute_rpp_artifact
)

EXPECTED_RULE_TAG="preliminary-v3-stratified-split-rule-freeze-2026-09-07"
EXPECTED_ARTIFACT_SHA="9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b"

def die(x): raise SystemExit("ERROR: "+x)
def git(repo,*args):
    p=subprocess.run(["git",*args],cwd=repo,capture_output=True,text=True)
    if p.returncode: die(p.stderr)
    return p.stdout.strip()

def load_artifact(path):
    d=json.loads(path.read_text(encoding="utf-8-sig"))
    protos=[]
    for p in d["prototypes"]:
        rels=[]
        for r in p["relations"]:
            rels.append(RPPRelation(
                feature_a=str(r["feature_a"]),feature_b=str(r["feature_b"]),
                direction=1 if r["direction"]==">" else -1,
                weight=float(r["weight"]),
                within_support=r.get("within_support"),contrast=r.get("contrast")
            ))
        protos.append(RPPPrototype(profile_id=p["profile_id"],relations=tuple(rels)))
    art=RelationalPatientProfileArtifact(
        artifact_id=d["artifact_id"],source_dataset=d["source_dataset"],
        feature_namespace=d["feature_namespace"],mapping_hash=d["mapping_hash"],
        preprocessing_hash=d["preprocessing_hash"],config_hash=d["config_hash"],
        software_commit=d["software_commit"],prototypes=tuple(protos),
        min_score=float(d["min_score"]),min_margin=float(d["min_margin"]),
        min_executable_coverage=float(d["min_executable_coverage"]),
        created_prelabel=bool(d["created_prelabel"]),
        transportability_certificate=d.get("transportability_certificate")
    )
    if art.sha256()!=EXPECTED_ARTIFACT_SHA: die("artifact SHA mismatch")
    return art

def rank_key(pid):
    return hashlib.sha256(f"RPPV3_STRAT_SPLIT|{pid}".encode()).hexdigest()

def choose_split(records):
    ranked=sorted(records,key=lambda r:(rank_key(r["participant_id"]),int(r["participant_id"])))
    for pos,r in enumerate(ranked): r["sha_rank"]=pos

    total_i=sum(r["contributes_index"] for r in ranked)
    total_c=sum(r["contributes_comparator"] for r in ranked)
    total_ap=sum(r["assigned_participant"] for r in ranked)
    total_as=sum(r["n_assigned_specimens"] for r in ranked)

    states={(0,0,0,0,0):(0,())}
    for r in ranked:
        new=dict(states)
        pid=r["participant_id"]; rank=r["sha_rank"]
        di=int(r["contributes_index"]); dc=int(r["contributes_comparator"])
        dap=int(r["assigned_participant"]); das=int(r["n_assigned_specimens"])
        for st,(cost,tup) in states.items():
            n,si,sc,sap,sas=st
            if n>=12: continue
            ns=(n+1,si+di,sc+dc,sap+dap,sas+das)
            cand=(cost+rank,tup+(pid,))
            old=new.get(ns)
            if old is None or cand[0]<old[0] or (cand[0]==old[0] and cand[1]<old[1]):
                new[ns]=cand
        states=new

    feasible=[]
    for st,(cost,tup) in states.items():
        n,si,sc,sap,sas=st
        if n!=12: continue
        ri=total_i-si; rc=total_c-sc
        rap=total_ap-sap; ras=total_as-sas
        if si<6 or sc<6: continue
        if ri<8 or rc<8: continue
        if rap<13 or ras<26: continue
        feasible.append((cost,tup,st,(ri,rc,rap,ras)))
    if not feasible: return None,ranked,(total_i,total_c,total_ap,total_as)
    feasible.sort(key=lambda z:(z[0],z[1]))
    return feasible[0],ranked,(total_i,total_c,total_ap,total_as)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--prep-root",required=True); a=ap.parse_args()
    root=Path(a.prep_root).resolve(); repo=root/"repo_patch/_work_repo"

    # Implementation-only repair: after the V3 rule itself has been frozen, later
    # runner commits are allowed as long as the frozen rule commit remains an ancestor.
    rule_commit=git(repo,"rev-list","-n","1",EXPECTED_RULE_TAG)
    head=git(repo,"rev-parse","HEAD")
    anc=subprocess.run(["git","merge-base","--is-ancestor",rule_commit,head],cwd=repo)
    if anc.returncode!=0:
        die(f"frozen V3 rule commit {rule_commit} is not an ancestor of current HEAD {head}")
    if git(repo,"status","--porcelain"): die("repo dirty")

    cache=np.load(root/"SOURCE_AUDIT_V1/GSE19804_COMMON_ENTREZ_V1.npz")
    X=cache["X"].astype(float); gsms=[str(x) for x in cache["gsms"]]; genes=[str(x) for x in cache["gene_ids"]]
    with (root/"_stage_packages_v2/pre_freeze/metadata/GSE19804_participant_map_BLIND.csv").open(encoding="utf-8-sig",newline="") as f:
        m=list(csv.DictReader(f))
    with (repo/"docs/preliminary/SOURCE_PARTICIPANT_SPLIT_V1.csv").open(encoding="utf-8-sig",newline="") as f:
        split=list(csv.DictReader(f))
    art=load_artifact(repo/"docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json")
    hold=[r["participant_id"] for r in split if r["role"]!="ASSIGNMENT_LEARNING"]
    gsm_to_i={g:i for i,g in enumerate(gsms)}
    meta={r["gsm"]:r for r in m}
    rows=[]; row_gsms=[]
    for pid in hold:
        gg=sorted([r["gsm"] for r in m if r["participant_id"]==pid])
        for g in gg: rows.append(gsm_to_i[g]); row_gsms.append(g)
    ex=execute_rpp_artifact(X[np.array(rows,int)],genes,art)
    rec=[]
    for pid in hold:
        inds=[i for i,g in enumerate(row_gsms) if meta[g]["participant_id"]==pid]
        ee=[ex[i] for i in inds]
        ni=sum(e.assignment=="ASSIGNED" and e.best_profile_id=="P0" for e in ee)
        nc=sum(e.assignment=="ASSIGNED" and e.best_profile_id!="P0" for e in ee)
        na=sum(e.assignment=="ASSIGNED" for e in ee)
        rec.append({"participant_id":pid,"contributes_index":int(ni>0),
                    "contributes_comparator":int(nc>0),
                    "assigned_participant":int(na>0),"n_assigned_specimens":int(na)})
    chosen,ranked,tot=choose_split(rec)
    out=root/"SOURCE_FIT_V3_PREP"; out.mkdir(exist_ok=True)
    with (out/"V3_HOLDOUT_ASSIGNMENT_STRATA.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["participant_id","contributes_index","contributes_comparator","assigned_participant","n_assigned_specimens","sha_rank"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(ranked)
    if chosen is None:
        stat={"schema":"V3SplitPreparation/v1","status":"NO_FEASIBLE_ASSIGNMENT_STRATIFIED_SPLIT",
              "total_index_contributors":tot[0],"total_comparator_contributors":tot[1],
              "total_assigned_participants":tot[2],"total_assigned_specimens":tot[3],
              "target_values_opened":False,"target_labels_opened":False}
        (out/"V3_SPLIT_STATUS.json").write_text(json.dumps(stat,indent=2)+"\n")
        print("V3 SPLIT: NO FEASIBLE 12/18 ASSIGNMENT-STRATIFIED SPLIT")
        print(json.dumps(stat,indent=2))
        return
    cost,tup,st,refstats=chosen
    sig=set(tup)
    splitrows=[]
    for r in ranked:
        splitrows.append({"participant_id":r["participant_id"],
                          "role":"SIGNATURE_CONSTRUCTION_V3" if r["participant_id"] in sig else "SOURCE_REFERENCE_V3",
                          "sha_rank":r["sha_rank"],
                          "contributes_index":r["contributes_index"],
                          "contributes_comparator":r["contributes_comparator"],
                          "assigned_participant":r["assigned_participant"],
                          "n_assigned_specimens":r["n_assigned_specimens"]})
    with (out/"SOURCE_PARTICIPANT_SPLIT_V3.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(splitrows[0])); w.writeheader(); w.writerows(splitrows)
    n,si,sc,sap,sas=st; ri,rc,rap,ras=refstats
    stat={"schema":"V3SplitPreparation/v1","status":"SPLIT_READY","signature_n":12,"reference_n":18,
          "signature_index_contributors":si,"signature_comparator_contributors":sc,
          "reference_index_contributors":ri,"reference_comparator_contributors":rc,
          "reference_assigned_participants":rap,"reference_assigned_specimens":ras,
          "objective_rank_sum":cost,"target_values_opened":False,"target_labels_opened":False}
    (out/"V3_SPLIT_STATUS.json").write_text(json.dumps(stat,indent=2)+"\n")
    print("V3 SPLIT READY")
    print(json.dumps(stat,indent=2))
    print("TARGET VALUES OPENED: NO")
    print("TARGET LABELS OPENED: NO")
if __name__=="__main__": main()
