#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, subprocess, sys
from pathlib import Path
import numpy as np

from relational_patient_profiles.artifact import artifact_from_rr_direct, execute_rpp_artifact

EXPECTED_V2_TAG = "preliminary-v2-source-gate-freeze-2026-09-07"

def die(msg): raise SystemExit("ERROR: "+msg)
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
    spec=importlib.util.spec_from_file_location("v1src",p)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def null_complete_discovery(Xdisc, block_rows, rep):
    rng=np.random.default_rng(np.random.SeedSequence([20260919,int(rep)]))
    out=np.empty_like(Xdisc)
    nb=len(block_rows)
    for j in range(Xdisc.shape[1]):
        perm=rng.permutation(nb)
        swaps=rng.random(nb)<0.5
        for d,s in enumerate(perm):
            src=block_rows[s]
            if swaps[d]: src=src[::-1]
            out[block_rows[d],j]=Xdisc[src,j]
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--prep-root",required=True)
    a=ap.parse_args()
    prep=Path(a.prep_root).resolve()
    repo=prep/"repo_patch/_work_repo"
    out=prep/"SOURCE_FIT_V2"; out.mkdir(exist_ok=True)
    audit=prep/"SOURCE_AUDIT_V1"
    stage1=prep/"_stage_packages_v2/pre_freeze"

    if git(repo,"status","--porcelain"): die("repo dirty before V2 run")
    head=git(repo,"rev-parse","HEAD")
    tag=git(repo,"rev-list","-n","1",EXPECTED_V2_TAG)
    if head!=tag: die(f"HEAD {head} != V2 freeze tag commit {tag}")

    v1=load_v1(repo)
    cfg=json.loads((repo/"docs/preliminary/V1_CONFIG.json").read_text(encoding="utf-8-sig"))
    v2cfg=json.loads((repo/"docs/preliminary/V2_SOURCE_GATE_CONFIG.json").read_text(encoding="utf-8-sig"))

    # Verify exact V1 outcome files.
    v1res=prep/"SOURCE_FIT_V1/source_fit_summary.csv"
    v1status=prep/"SOURCE_FIT_V1/SOURCE_RESULT_STATUS.json"
    if not v1res.exists() or not v1status.exists(): die("missing local V1 source result")
    rows=list(csv.DictReader(v1res.open(encoding="utf-8-sig",newline="")))
    r2=next(r for r in rows if int(r["k"])==2)
    r3=next(r for r in rows if int(r["k"])==3)
    if r2["passes"].lower()!="false" or r3["passes"].lower()!="false": die("V1 closure mismatch")
    if not (float(r2["q_observed"])>float(r2["q_null95"])): die("V1 K2 Q gate was not passed")
    if int(r2["min_profile_participant_support_observed"])<8: die("V1 K2 support precondition failed")
    if int(r2["stability_successes"])!=20: die("V1 K2 stability success precondition failed")
    v1stat=json.loads(v1status.read_text())
    if v1stat.get("status")!="NO_STABLE_STRUCTURE": die("V1 status mismatch")

    cache=np.load(audit/"GSE19804_COMMON_ENTREZ_V1.npz")
    X=cache["X"].astype(float)
    gsms=[str(x) for x in cache["gsms"]]
    genes=np.array([str(x) for x in cache["gene_ids"]])
    pools=np.array([str(x) for x in cache["pools"]])

    with (stage1/"metadata/GSE19804_participant_map_BLIND.csv").open(encoding="utf-8-sig",newline="") as f:
        map_rows=list(csv.DictReader(f))
    with (repo/"docs/preliminary/SOURCE_PARTICIPANT_SPLIT_V1.csv").open(encoding="utf-8-sig",newline="") as f:
        split=list(csv.DictReader(f))

    gsm_to_row={g:i for i,g in enumerate(gsms)}
    p_to_rows={}
    for r in map_rows: p_to_rows.setdefault(r["participant_id"],[]).append(gsm_to_row[r["gsm"]])
    for p in p_to_rows: p_to_rows[p]=sorted(p_to_rows[p])

    assign_pids=[r["participant_id"] for r in split if r["role"]=="ASSIGNMENT_LEARNING"]
    sig_pids=[r["participant_id"] for r in split if r["role"]=="SIGNATURE_CONSTRUCTION"]
    ref_pids=[r["participant_id"] for r in split if r["role"]=="SOURCE_REFERENCE"]
    assign_rows=np.array([i for p in assign_pids for i in p_to_rows[p]],int)
    sig_rows=np.array([i for p in sig_pids for i in p_to_rows[p]],int)
    ref_rows=np.array([i for p in ref_pids for i in p_to_rows[p]],int)
    assign_blocks=np.arange(60).reshape(30,2)

    disc_idx=np.flatnonzero(pools=="DISCOVERY")
    val_mask=pools=="VALIDATION"
    Xdisc=X[assign_rows][:,disc_idx]
    disc_genes=genes[disc_idx]

    # Observed full K2 and exact 20 stability values.
    top=v1.top_genes_by_mad(X,genes,pools=="DISCOVERY",assign_rows,60)
    feat=genes[top].tolist()
    full,B=v1.fit_strict(X[assign_rows][:,top],feat,2,
        max_rules=cfg["assignment"]["max_rules_per_profile"],
        min_support=cfg["assignment"]["min_rule_support"],
        min_contrast=cfg["assignment"]["min_rule_contrast"],
        max_iter=cfg["assignment"]["max_iter"])
    full_labels=np.array(full.labels,int)

    obs=[]
    subs=v1.deterministic_subsets(assign_pids)
    for rr,subp in enumerate(subs):
        subrows=np.array([i for p in subp for i in p_to_rows[p]],int)
        sidx=v1.top_genes_by_mad(X,genes,pools=="DISCOVERY",subrows,60)
        sids=genes[sidx].tolist()
        try:
            sf,_=v1.fit_strict(X[subrows][:,sidx],sids,2,
                max_rules=cfg["assignment"]["max_rules_per_profile"],
                min_support=cfg["assignment"]["min_rule_support"],
                min_contrast=cfg["assignment"]["min_rule_contrast"],
                max_iter=cfg["assignment"]["max_iter"])
            pred=v1.forced_labels(X[assign_rows][:,sidx],sids,sf.prototypes)
            aa=v1.ari(full_labels,pred); ok=True
        except Exception:
            aa=0.0; ok=False
        obs.append({"replicate":rr,"ari":float(aa),"fit_ok":ok})
    Sobs=float(np.mean([x["ari"] for x in obs]))
    if abs(Sobs-float(r2["stability_mean_ari"]))>1e-10:
        die(f"Observed stability does not reproduce V1: {Sobs} vs {r2['stability_mean_ari']}")

    # Null stability distribution. Full DISCOVERY pool is permuted.
    nullrows=[]
    for rep in range(v2cfg["null_stability_replicates"]):
        Xn=null_complete_discovery(Xdisc,assign_blocks,rep)
        try:
            full_top_local=v1.top_genes_by_mad(Xn,disc_genes,np.ones(len(disc_genes),dtype=bool),
                                                np.arange(60),60)
            full_ids=disc_genes[full_top_local].tolist()
            nf,_=v1.fit_strict(Xn[:,full_top_local],full_ids,2,
                max_rules=cfg["assignment"]["max_rules_per_profile"],
                min_support=cfg["assignment"]["min_rule_support"],
                min_contrast=cfg["assignment"]["min_rule_contrast"],
                max_iter=cfg["assignment"]["max_iter"])
            nfull=np.array(nf.labels,int)
            aris=[]; succ=0
            for subp in subs:
                # local source-assignment rows for these participant IDs
                pos=[]
                for p in subp:
                    bi=assign_pids.index(p)
                    pos.extend(assign_blocks[bi].tolist())
                pos=np.array(pos,int)
                si=v1.top_genes_by_mad(Xn,disc_genes,np.ones(len(disc_genes),dtype=bool),pos,60)
                ids=disc_genes[si].tolist()
                try:
                    sf,_=v1.fit_strict(Xn[pos][:,si],ids,2,
                        max_rules=cfg["assignment"]["max_rules_per_profile"],
                        min_support=cfg["assignment"]["min_rule_support"],
                        min_contrast=cfg["assignment"]["min_rule_contrast"],
                        max_iter=cfg["assignment"]["max_iter"])
                    pred=v1.forced_labels(Xn[:,si],ids,sf.prototypes)
                    aris.append(v1.ari(nfull,pred)); succ+=1
                except Exception:
                    aris.append(0.0)
            mean=float(np.mean(aris)); fullok=True
        except Exception:
            mean=0.0; succ=0; fullok=False
        nullrows.append({"replicate":rep,"mean_stability_ari":mean,
                         "full_fit_ok":fullok,"subsample_fit_successes":succ})

    nullv=np.array([r["mean_stability_ari"] for r in nullrows],float)
    null95=float(np.quantile(nullv,.95,method="higher"))
    pemp=float((1+np.sum(nullv>=Sobs))/(len(nullv)+1))
    gate=bool(
        full.converged and
        int(r2["min_profile_participant_support_observed"])>=8 and
        float(r2["q_observed"])>float(r2["q_null95"]) and
        all(x["fit_ok"] for x in obs) and
        Sobs>null95 and pemp<=0.05
    )

    with (out/"v2_observed_stability.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["replicate","ari","fit_ok"]); w.writeheader(); w.writerows(obs)
    with (out/"v2_null_stability.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["replicate","mean_stability_ari","full_fit_ok","subsample_fit_successes"])
        w.writeheader(); w.writerows(nullrows)

    gate_summary={
        "schema":"V2SourceGateResult/v1","software_commit":head,
        "S_observed":Sobs,"S_null95":null95,"empirical_p":pemp,
        "v1_q_observed":float(r2["q_observed"]),"v1_q_null95":float(r2["q_null95"]),
        "min_profile_participant_support":int(r2["min_profile_participant_support_observed"]),
        "observed_stability_successes":sum(x["fit_ok"] for x in obs),
        "v2_source_gate_pass":gate,
        "target_values_opened":False,"target_labels_opened":False
    }
    (out/"V2_SOURCE_GATE_RESULT.json").write_text(json.dumps(gate_summary,indent=2)+"\n")

    if not gate:
        status={**gate_summary,"status":"V2_NO_STABLE_STRUCTURE","source_reference_adequate":False}
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print("V2 SOURCE GATE: FAIL")
        print(f"S_obs={Sobs:.6f}  S_null95={null95:.6f}  empirical_p={pemp:.6f}")
        print("STATUS: V2_NO_STABLE_STRUCTURE")
        print("TARGET VALUES OPENED: NO")
        print("TARGET LABELS OPENED: NO")
        return

    # Artifact
    mapping_hash=sh(audit/"COMMON_ENTREZ_UNIVERSE_V1.csv")
    preprocessing_hash=v1.canonical_matrix_hash(X,gsms,genes)
    config_hash=sh(repo/"docs/preliminary/V1_CONFIG.json")
    art=artifact_from_rr_direct(full,artifact_id="GSE19804_RR_DIRECT_V2_K2",
        source_dataset="GSE19804",feature_namespace="Entrez Gene ID",
        mapping_hash=mapping_hash,preprocessing_hash=preprocessing_hash,
        config_hash=config_hash,software_commit=head,
        min_score=cfg["assignment"]["artifact_min_score"],
        min_margin=cfg["assignment"]["artifact_min_margin"],
        min_executable_coverage=cfg["assignment"]["artifact_min_executable_coverage"])
    ad=art.to_dict(); ad["artifact_sha256"]=art.sha256()
    (out/"SOURCE_ASSIGNMENT_ARTIFACT_V2.json").write_text(json.dumps(ad,indent=2)+"\n")

    pscore={}
    for p in art.prototypes:
        vv=[float(r.within_support or 0)*max(float(r.contrast or 0),0) for r in p.relations]
        pscore[p.profile_id]=float(np.mean(vv))
    index_profile=sorted(pscore,key=lambda p:(-pscore[p],p))[0]

    # Signature construction exactly as V1.
    sig_exec=execute_rpp_artifact(X[sig_rows],genes.tolist(),art)
    sig_gsms=[gsms[i] for i in sig_rows]
    meta_by_gsm={r["gsm"]:r for r in map_rows}
    sig_map=[meta_by_gsm[g] for g in sig_gsms]
    idxp=set(); compp=set()
    for ex,m in zip(sig_exec,sig_map):
        if ex.assignment!="ASSIGNED": continue
        (idxp if ex.best_profile_id==index_profile else compp).add(m["participant_id"])
    if len(idxp)<cfg["signature"]["min_group_participants_signature"] or len(compp)<cfg["signature"]["min_group_participants_signature"]:
        status={**gate_summary,"status":"INSUFFICIENT_SOURCE_SIGNATURE",
                "selected_k":2,"index_profile":index_profile,
                "signature_index_participants":len(idxp),
                "signature_comparator_participants":len(compp),
                "source_reference_adequate":False}
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print("V2 SOURCE GATE: PASS")
        print(f"S_obs={Sobs:.6f}  S_null95={null95:.6f}  empirical_p={pemp:.6f}")
        print(f"STATUS: INSUFFICIENT_SOURCE_SIGNATURE groups={len(idxp)}/{len(compp)}")
        print("TARGET VALUES OPENED: NO"); print("TARGET LABELS OPENED: NO"); return

    val_idx=v1.top_genes_by_mad(X,genes,val_mask,sig_rows,cfg["signature"]["validation_feature_budget"])
    val_genes=genes[val_idx].tolist()
    gene_index={g:i for i,g in enumerate(genes)}
    pairs=[(aa,bb) for ii,aa in enumerate(val_genes) for bb in val_genes[ii+1:]]
    candidates=[]; all_sig_pids=sorted(set(sig_pids),key=int)
    Xsig=X[sig_rows]
    for aa,bb in pairs:
        best=None
        for d in (1,-1):
            vals=v1.relation_participant_values(Xsig,gene_index,(aa,bb),d,sig_exec,sig_map,sig_gsms,index_profile)
            if not vals["INDEX"] or not vals["COMPARATOR"]: continue
            support=float(np.mean(list(vals["INDEX"].values())))
            contrast=v1.contrast_from_group_values(vals)
            if best is None or contrast>best["contrast"] or (contrast==best["contrast"] and d==1):
                best={"a":aa,"b":bb,"direction":d,"support":support,"contrast":contrast,"vals":vals}
        if best is None: continue
        if best["support"]<cfg["signature"]["min_index_support"] or best["contrast"]<cfg["signature"]["min_observed_contrast"]: continue
        boot=v1.bootstrap_contrast(best["vals"],all_sig_pids,cfg["signature"]["bootstrap_replicates"],20260913)
        pos=float(np.mean(boot>0)); med=float(np.median(boot))
        if pos<cfg["signature"]["min_positive_bootstrap_fraction"]: continue
        best["bootstrap_positive_fraction"]=pos; best["bootstrap_median_contrast"]=med
        candidates.append(best)
    candidates.sort(key=lambda z:(-z["bootstrap_median_contrast"],-z["contrast"],z["a"],z["b"],-z["direction"]))
    reuse={}; selected=[]
    for c in candidates:
        if reuse.get(c["a"],0)>=cfg["signature"]["max_feature_reuse"] or reuse.get(c["b"],0)>=cfg["signature"]["max_feature_reuse"]: continue
        selected.append(c); reuse[c["a"]]=reuse.get(c["a"],0)+1; reuse[c["b"]]=reuse.get(c["b"],0)+1
        if len(selected)>=cfg["signature"]["max_relations"]: break

    with (out/"source_signature.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["feature_a","feature_b","direction","index_support","observed_contrast","bootstrap_positive_fraction","bootstrap_median_contrast"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for c in selected:
            w.writerow({"feature_a":c["a"],"feature_b":c["b"],"direction":">" if c["direction"]==1 else "<",
                        "index_support":c["support"],"observed_contrast":c["contrast"],
                        "bootstrap_positive_fraction":c["bootstrap_positive_fraction"],
                        "bootstrap_median_contrast":c["bootstrap_median_contrast"]})
    if len(selected)<cfg["signature"]["min_relations"]:
        status={**gate_summary,"status":"INSUFFICIENT_SOURCE_SIGNATURE","selected_k":2,
                "index_profile":index_profile,"selected_signature_relations":len(selected),
                "source_reference_adequate":False}
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print("V2 SOURCE GATE: PASS")
        print(f"S_obs={Sobs:.6f}  S_null95={null95:.6f}  empirical_p={pemp:.6f}")
        print(f"STATUS: INSUFFICIENT_SOURCE_SIGNATURE relations={len(selected)}")
        print("TARGET VALUES OPENED: NO"); print("TARGET LABELS OPENED: NO"); return

    sigobj={"schema":"ValidationSignatureV2","source_dataset":"GSE19804","index_profile":index_profile,
            "software_commit":head,"relations":[{"feature_a":c["a"],"feature_b":c["b"],
            "direction":">" if c["direction"]==1 else "<","weight":1.0} for c in selected]}
    raw=(json.dumps(sigobj,sort_keys=True,separators=(",",":"))+"\n").encode()
    sigobj["sha256"]=hashlib.sha256(raw).hexdigest()
    (out/"VALIDATION_SIGNATURE_V2.json").write_text(json.dumps(sigobj,indent=2)+"\n")

    # Source reference.
    ref_exec=execute_rpp_artifact(X[ref_rows],genes.tolist(),art)
    ref_gsms=[gsms[i] for i in ref_rows]; ref_map=[meta_by_gsm[g] for g in ref_gsms]
    vs=[]
    for gi in ref_rows:
        sats=[]
        for c in selected:
            ia=gene_index[c["a"]]; ib=gene_index[c["b"]]
            sats.append(float(X[gi,ia]>X[gi,ib]) if c["direction"]==1 else float(X[gi,ia]<X[gi,ib]))
        vs.append(float(np.mean(sats)))
    vals=v1.participant_group_values(ref_exec,vs,ref_map,ref_gsms,index_profile)
    C=v1.contrast_from_group_values(vals)
    boot=v1.bootstrap_contrast(vals,sorted(set(ref_pids),key=int),cfg["source_reference"]["bootstrap_replicates"],20260914)
    lower=float(np.quantile(boot,cfg["source_reference"]["source_one_sided_alpha"],method="lower"))
    specimen_cov=float(np.mean([e.assignment=="ASSIGNED" for e in ref_exec]))
    covered=set(m["participant_id"] for e,m in zip(ref_exec,ref_map) if e.assignment=="ASSIGNED")
    participant_cov=len(covered)/len(ref_pids)
    ni=len(vals["INDEX"]); nc=len(vals["COMPARATOR"])
    adequate=bool(specimen_cov>=cfg["source_reference"]["min_specimen_assignment_coverage"] and
                  participant_cov>=cfg["source_reference"]["min_participant_assignment_coverage"] and
                  ni>=cfg["source_reference"]["min_group_participants"] and
                  nc>=cfg["source_reference"]["min_group_participants"] and
                  C>=cfg["source_reference"]["source_contrast_floor"] and lower>0)
    refsum={"selected_k":2,"index_profile":index_profile,"C_source":C,
            "one_sided_95_lower":lower,"specimen_assignment_coverage":specimen_cov,
            "participant_assignment_coverage":participant_cov,"index_participants":ni,
            "comparator_participants":nc,"source_reference_adequate":adequate}
    with (out/"source_reference.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(refsum)); w.writeheader(); w.writerow(refsum)

    status={**gate_summary,
            "status":"SOURCE_REFERENCE_PASS" if adequate else "INSUFFICIENT_SOURCE_REFERENCE",
            "selected_k":2,"index_profile":index_profile,"artifact_sha256":art.sha256(),
            "signature_sha256":sigobj["sha256"],"selected_signature_relations":len(selected),
            "C_source":C,"source_reference_lower95":lower,"source_reference_adequate":adequate}
    (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")

    # Result manifest.
    with (out/"SOURCE_RESULTS_SHA256.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["filename","bytes","sha256"])
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name!="SOURCE_RESULTS_SHA256.csv":
                w.writerow([p.name,p.stat().st_size,sh(p)])

    print("V2 SOURCE GATE: PASS")
    print(f"S_obs={Sobs:.6f}  S_null95={null95:.6f}  empirical_p={pemp:.6f}")
    print(f"Index profile: {index_profile}")
    print(f"Signature relations: {len(selected)}")
    print(f"C_source={C:.6f}  one-sided lower95={lower:.6f}")
    print(f"Reference support index/comparator={ni}/{nc}")
    print(f"Assignment coverage specimen/participant={specimen_cov:.3f}/{participant_cov:.3f}")
    print(f"STATUS: {status['status']}")
    print("TARGET VALUES OPENED: NO")
    print("TARGET LABELS OPENED: NO")

if __name__=="__main__": main()
