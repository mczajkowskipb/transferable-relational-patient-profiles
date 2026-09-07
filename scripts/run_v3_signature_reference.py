#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,subprocess
from pathlib import Path
import numpy as np
from relational_patient_profiles.artifact import (
    RPPRelation,RPPPrototype,RelationalPatientProfileArtifact,execute_rpp_artifact
)

EXPECTED_SPLIT_TAG="preliminary-v3-source-split-freeze-2026-09-07"
EXPECTED_ARTIFACT_SHA="9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b"

def die(x): raise SystemExit("ERROR: "+x)
def sh(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()
def git(repo,*args):
    p=subprocess.run(["git",*args],cwd=repo,capture_output=True,text=True)
    if p.returncode: die(p.stderr)
    return p.stdout.strip()
def load_v1(repo):
    p=repo/"scripts/run_source_fit_v1.py"
    s=importlib.util.spec_from_file_location("v1src",p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def load_artifact(path):
    d=json.loads(path.read_text(encoding="utf-8-sig")); ps=[]
    for p in d["prototypes"]:
        rs=[RPPRelation(str(r["feature_a"]),str(r["feature_b"]),1 if r["direction"]==">" else -1,
                        float(r["weight"]),r.get("within_support"),r.get("contrast")) for r in p["relations"]]
        ps.append(RPPPrototype(p["profile_id"],tuple(rs)))
    art=RelationalPatientProfileArtifact(
        artifact_id=d["artifact_id"],source_dataset=d["source_dataset"],feature_namespace=d["feature_namespace"],
        mapping_hash=d["mapping_hash"],preprocessing_hash=d["preprocessing_hash"],config_hash=d["config_hash"],
        software_commit=d["software_commit"],prototypes=tuple(ps),min_score=float(d["min_score"]),
        min_margin=float(d["min_margin"]),min_executable_coverage=float(d["min_executable_coverage"]),
        created_prelabel=bool(d["created_prelabel"]),transportability_certificate=d.get("transportability_certificate"))
    if art.sha256()!=EXPECTED_ARTIFACT_SHA: die("artifact SHA mismatch")
    return art

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--prep-root",required=True); a=ap.parse_args()
    root=Path(a.prep_root).resolve(); repo=root/"repo_patch/_work_repo"
    if git(repo,"rev-parse","HEAD")!=git(repo,"rev-list","-n","1",EXPECTED_SPLIT_TAG): die("HEAD is not exact V3 split-freeze commit")
    if git(repo,"status","--porcelain"): die("repo dirty")
    v1=load_v1(repo); cfg=json.loads((repo/"docs/preliminary/V1_CONFIG.json").read_text(encoding="utf-8-sig"))
    art=load_artifact(repo/"docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json")
    splitfile=repo/"docs/preliminary/v3_source_split/SOURCE_PARTICIPANT_SPLIT_V3.csv"
    split=list(csv.DictReader(splitfile.open(encoding="utf-8-sig",newline="")))
    sigp=[r["participant_id"] for r in split if r["role"]=="SIGNATURE_CONSTRUCTION_V3"]
    refp=[r["participant_id"] for r in split if r["role"]=="SOURCE_REFERENCE_V3"]
    if len(sigp)!=12 or len(refp)!=18: die("V3 split sizes wrong")

    cache=np.load(root/"SOURCE_AUDIT_V1/GSE19804_COMMON_ENTREZ_V1.npz")
    X=cache["X"].astype(float); gsms=[str(x) for x in cache["gsms"]]
    genes=np.array([str(x) for x in cache["gene_ids"]]); pools=np.array([str(x) for x in cache["pools"]])
    with (root/"_stage_packages_v2/pre_freeze/metadata/GSE19804_participant_map_BLIND.csv").open(encoding="utf-8-sig",newline="") as f:
        maprows=list(csv.DictReader(f))
    gsm2i={g:i for i,g in enumerate(gsms)}; meta={r["gsm"]:r for r in maprows}
    p2rows={}
    for r in maprows: p2rows.setdefault(r["participant_id"],[]).append(gsm2i[r["gsm"]])
    for p in p2rows: p2rows[p]=sorted(p2rows[p])
    sigrows=np.array([i for p in sigp for i in p2rows[p]],int)
    refrows=np.array([i for p in refp for i in p2rows[p]],int)
    siggs=[gsms[i] for i in sigrows]; refgs=[gsms[i] for i in refrows]
    sigmap=[meta[g] for g in siggs]; refmap=[meta[g] for g in refgs]
    sigex=execute_rpp_artifact(X[sigrows],genes.tolist(),art)
    refex=execute_rpp_artifact(X[refrows],genes.tolist(),art)
    index="P0"

    idxp={m["participant_id"] for e,m in zip(sigex,sigmap) if e.assignment=="ASSIGNED" and e.best_profile_id==index}
    compp={m["participant_id"] for e,m in zip(sigex,sigmap) if e.assignment=="ASSIGNED" and e.best_profile_id!=index}
    out=root/"SOURCE_FIT_V3"; out.mkdir(exist_ok=True)
    if len(idxp)<6 or len(compp)<6:
        stat={"status":"INSUFFICIENT_SOURCE_SIGNATURE_SUPPORT_AFTER_FROZEN_SPLIT",
              "signature_index_participants":len(idxp),"signature_comparator_participants":len(compp),
              "target_values_opened":False,"target_labels_opened":False}
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(stat,indent=2)+"\n")
        print(json.dumps(stat,indent=2)); return

    validx=v1.top_genes_by_mad(X,genes,pools=="VALIDATION",sigrows,cfg["signature"]["validation_feature_budget"])
    valgenes=genes[validx].tolist(); gidx={g:i for i,g in enumerate(genes)}
    pairs=[(a,b) for ii,a in enumerate(valgenes) for b in valgenes[ii+1:]]
    cand=[]; allsig=sorted(sigp,key=int)
    Xsig=X[sigrows]
    for aa,bb in pairs:
        best=None
        for d in (1,-1):
            vals=v1.relation_participant_values(Xsig,gidx,(aa,bb),d,sigex,sigmap,siggs,index)
            if not vals["INDEX"] or not vals["COMPARATOR"]: continue
            sup=float(np.mean(list(vals["INDEX"].values()))); con=v1.contrast_from_group_values(vals)
            if best is None or con>best["contrast"] or (con==best["contrast"] and d==1):
                best={"a":aa,"b":bb,"direction":d,"support":sup,"contrast":con,"vals":vals}
        if best is None: continue
        if best["support"]<cfg["signature"]["min_index_support"] or best["contrast"]<cfg["signature"]["min_observed_contrast"]: continue
        boot=v1.bootstrap_contrast(best["vals"],allsig,cfg["signature"]["bootstrap_replicates"],20260913)
        pos=float(np.mean(boot>0)); med=float(np.median(boot))
        if pos<cfg["signature"]["min_positive_bootstrap_fraction"]: continue
        best["bootstrap_positive_fraction"]=pos; best["bootstrap_median_contrast"]=med; cand.append(best)
    cand.sort(key=lambda z:(-z["bootstrap_median_contrast"],-z["contrast"],z["a"],z["b"],-z["direction"]))
    reuse={}; sel=[]
    for c in cand:
        if reuse.get(c["a"],0)>=cfg["signature"]["max_feature_reuse"] or reuse.get(c["b"],0)>=cfg["signature"]["max_feature_reuse"]: continue
        sel.append(c); reuse[c["a"]]=reuse.get(c["a"],0)+1; reuse[c["b"]]=reuse.get(c["b"],0)+1
        if len(sel)>=cfg["signature"]["max_relations"]: break
    with (out/"source_signature.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["feature_a","feature_b","direction","index_support","observed_contrast","bootstrap_positive_fraction","bootstrap_median_contrast"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for c in sel: w.writerow({"feature_a":c["a"],"feature_b":c["b"],"direction":">" if c["direction"]==1 else "<",
            "index_support":c["support"],"observed_contrast":c["contrast"],
            "bootstrap_positive_fraction":c["bootstrap_positive_fraction"],
            "bootstrap_median_contrast":c["bootstrap_median_contrast"]})
    if len(sel)<cfg["signature"]["min_relations"]:
        stat={"status":"INSUFFICIENT_SOURCE_SIGNATURE_RELATIONS","selected_signature_relations":len(sel),
              "target_values_opened":False,"target_labels_opened":False}
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(stat,indent=2)+"\n")
        print(json.dumps(stat,indent=2)); return

    sigobj={"schema":"ValidationSignatureV3","source_dataset":"GSE19804","index_profile":index,
            "artifact_sha256":art.sha256(),"split_sha256":sh(splitfile),
            "relations":[{"feature_a":c["a"],"feature_b":c["b"],"direction":">" if c["direction"]==1 else "<","weight":1.0} for c in sel]}
    raw=(json.dumps(sigobj,sort_keys=True,separators=(",",":"))+"\n").encode()
    sigobj["sha256"]=hashlib.sha256(raw).hexdigest()
    (out/"VALIDATION_SIGNATURE_V3.json").write_text(json.dumps(sigobj,indent=2)+"\n")

    vs=[]
    for gi in refrows:
        ss=[]
        for c in sel:
            ia=gidx[c["a"]]; ib=gidx[c["b"]]
            ss.append(float(X[gi,ia]>X[gi,ib]) if c["direction"]==1 else float(X[gi,ia]<X[gi,ib]))
        vs.append(float(np.mean(ss)))
    vals=v1.participant_group_values(refex,vs,refmap,refgs,index)
    C=v1.contrast_from_group_values(vals)
    boot=v1.bootstrap_contrast(vals,sorted(refp,key=int),cfg["source_reference"]["bootstrap_replicates"],20260914)
    low=float(np.quantile(boot,cfg["source_reference"]["source_one_sided_alpha"],method="lower"))
    spec=float(np.mean([e.assignment=="ASSIGNED" for e in refex]))
    covp={m["participant_id"] for e,m in zip(refex,refmap) if e.assignment=="ASSIGNED"}
    pcov=len(covp)/len(refp); ni=len(vals["INDEX"]); nc=len(vals["COMPARATOR"])
    adequate=bool(spec>=cfg["source_reference"]["min_specimen_assignment_coverage"] and
                  pcov>=cfg["source_reference"]["min_participant_assignment_coverage"] and
                  ni>=cfg["source_reference"]["min_group_participants"] and
                  nc>=cfg["source_reference"]["min_group_participants"] and
                  C>=cfg["source_reference"]["source_contrast_floor"] and low>0)
    refsum={"C_source":C,"one_sided_95_lower":low,"specimen_assignment_coverage":spec,
            "participant_assignment_coverage":pcov,"index_participants":ni,
            "comparator_participants":nc,"source_reference_adequate":adequate}
    with (out/"source_reference.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(refsum)); w.writeheader(); w.writerow(refsum)
    with (out/"source_reference_assignments.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["gsm","participant_id","assignment","best_profile_id","best_score","margin","executable_coverage","validation_score"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for g,m,e,v in zip(refgs,refmap,refex,vs):
            w.writerow({"gsm":g,"participant_id":m["participant_id"],"assignment":e.assignment,
                        "best_profile_id":e.best_profile_id,"best_score":e.best_score,"margin":e.margin,
                        "executable_coverage":e.executable_coverage,"validation_score":v})
    status={"schema":"V3SourceOutcome/v1","status":"SOURCE_REFERENCE_PASS" if adequate else "INSUFFICIENT_SOURCE_REFERENCE",
            "artifact_sha256":art.sha256(),"signature_sha256":sigobj["sha256"],"split_sha256":sh(splitfile),
            "selected_signature_relations":len(sel),"C_source":C,"source_reference_lower95":low,
            "source_reference_adequate":adequate,"target_values_opened":False,"target_labels_opened":False}
    (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
    with (out/"SOURCE_RESULTS_SHA256.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["filename","bytes","sha256"])
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name!="SOURCE_RESULTS_SHA256.csv": w.writerow([p.name,p.stat().st_size,sh(p)])
    print("V3 SOURCE COMPLETE")
    print(f"Signature relations: {len(sel)}")
    print(f"C_source={C:.6f}; lower95={low:.6f}")
    print(f"Reference support index/comparator={ni}/{nc}")
    print(f"Coverage specimen/participant={spec:.3f}/{pcov:.3f}")
    print(f"STATUS: {status['status']}")
    print("TARGET VALUES OPENED: NO")
    print("TARGET LABELS OPENED: NO")
if __name__=="__main__": main()
