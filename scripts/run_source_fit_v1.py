#!/usr/bin/env python3
from __future__ import annotations

import argparse, csv, hashlib, json, math, subprocess
from dataclasses import replace
from pathlib import Path
from typing import Sequence
import numpy as np

from relational_patient_profiles.rr_direct import RRDirectResult, SparseRelationalPrototype
from relational_patient_profiles.artifact import (
    artifact_from_rr_direct, execute_rpp_artifact, RelationalPatientProfileArtifact
)

EXPECTED_PROTOCOL_FREEZE = "8eb98b22555816dc2771582797e73f068cff062d"
EXPECTED_IMPL_TAG = "preliminary-v1-source-implementation-2026-09-07"

def die(msg): raise SystemExit("ERROR: " + msg)

def sha256_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024), b""): h.update(b)
    return h.hexdigest()

def canonical_matrix_hash(X, gsms, genes):
    h=hashlib.sha256()
    h.update(np.asarray(X,dtype="<f4",order="C").tobytes(order="C"))
    h.update(b"\0GSMS\0"); h.update("\n".join(map(str,gsms)).encode())
    h.update(b"\0GENES\0"); h.update("\n".join(map(str,genes)).encode())
    return h.hexdigest()

def git(repo: Path, *args):
    p=subprocess.run(["git",*args], cwd=repo, capture_output=True, text=True)
    if p.returncode: die(f"git {' '.join(args)} failed:\n{p.stderr}")
    return p.stdout.strip()

def mad(v, axis=0):
    v=np.asarray(v,float)
    med=np.nanmedian(v,axis=axis,keepdims=True)
    return np.nanmedian(np.abs(v-med),axis=axis)

def ari(labels_a, labels_b):
    a=np.asarray(labels_a); b=np.asarray(labels_b)
    if a.size!=b.size: raise ValueError("ARI length mismatch")
    n=a.size
    if n<2: return 1.0
    ua, ia=np.unique(a,return_inverse=True); ub, ib=np.unique(b,return_inverse=True)
    cont=np.zeros((len(ua),len(ub)),dtype=np.int64)
    np.add.at(cont,(ia,ib),1)
    def c2(x): return x*(x-1)//2
    sum_comb=int(np.sum(c2(cont)))
    row=cont.sum(1); col=cont.sum(0)
    sum_row=int(np.sum(c2(row))); sum_col=int(np.sum(c2(col)))
    total=c2(n)
    expected=(sum_row*sum_col/total) if total else 0.0
    max_index=0.5*(sum_row+sum_col)
    den=max_index-expected
    if den==0: return 1.0
    return float((sum_comb-expected)/den)

def top_genes_by_mad(X, gene_ids, eligible_mask, rows, n=60):
    rows=np.asarray(rows,int)
    idx=np.flatnonzero(eligible_mask)
    vals=mad(X[rows][:,idx],axis=0)
    order=sorted(range(len(idx)), key=lambda q:(-float(vals[q]), int(gene_ids[idx[q]])))
    keep=idx[np.asarray(order[:n],int)]
    return keep

def pair_list(p):
    return [(i,j) for i in range(p) for j in range(i+1,p)]

def binary_forward_reverse(Xs, pairs):
    B=np.column_stack([(Xs[:,i] > Xs[:,j]).astype(np.uint8) for i,j in pairs])
    R=np.column_stack([(Xs[:,i] < Xs[:,j]).astype(np.uint8) for i,j in pairs])
    return B,R

def initial_labels(B,k):
    n=B.shape[0]
    prevalence=B.mean(0)
    central=np.mean(np.abs(B-prevalence),axis=1)
    seeds=[int(np.argmax(central))]
    while len(seeds)<k:
        dmin=np.min([np.mean(B != B[s],axis=1) for s in seeds],axis=0)
        dmin[seeds]=-1
        seeds.append(int(np.argmax(dmin)))
    D=np.column_stack([np.mean(B != B[s],axis=1) for s in seeds])
    return np.argmin(D,axis=1).astype(int)

def build_prototypes_strict(B,R,labels,pairs,feature_ids,k,max_rules,min_support,min_contrast):
    gp=B.mean(0); gr=R.mean(0); out=[]
    for c in range(k):
        mask=labels==c
        if mask.sum()==0: raise RuntimeError("empty cluster")
        p=B[mask].mean(0); rp=R[mask].mean(0)
        direction=(p>=rp).astype(np.uint8)
        support=np.where(direction==1,p,rp)
        other=np.where(direction==1,gp,gr)
        contrast=support-other
        cand=[j for j in range(B.shape[1]) if support[j]>=min_support and contrast[j]>=min_contrast]
        if not cand:
            raise RuntimeError(f"profile {c}: no rule meets frozen support/contrast thresholds")
        cand.sort(key=lambda j:(-float(contrast[j]),-float(support[j]),j))
        chosen=cand[:min(max_rules,len(cand))]
        rules=[]
        for j in chosen:
            aa,bb=pairs[j]
            rules.append((str(feature_ids[aa]),str(feature_ids[bb]),int(direction[j]),
                          float(support[j]),float(contrast[j])))
        out.append(SparseRelationalPrototype(cluster=c,rules=tuple(rules)))
    return tuple(out)

def score_prototypes(X,feature_ids,prototypes):
    X=np.asarray(X,float); index={str(f):j for j,f in enumerate(feature_ids)}
    scores=np.full((len(X),len(prototypes)),-np.inf,float)
    for ci,p in enumerate(prototypes):
        num=np.zeros(len(X)); den=np.zeros(len(X))
        for a,b,d,s,c in p.rules:
            xa=X[:,index[a]]; xb=X[:,index[b]]
            ok=np.isfinite(xa)&np.isfinite(xb)
            obs=(xa[ok]>xb[ok]) if d==1 else (xa[ok]<xb[ok])
            w=max(1e-12,s*max(c,1e-6))
            num[ok]+=w*obs; den[ok]+=w
        good=den>0; scores[good,ci]=num[good]/den[good]
    return scores

def repair_empty(labels,scores,k):
    labels=labels.copy()
    for c in range(k):
        if np.any(labels==c): continue
        counts=np.bincount(labels,minlength=k); donor=int(np.argmax(counts))
        di=np.flatnonzero(labels==donor)
        m=scores[di,donor]-np.max(np.delete(scores[di],donor,axis=1),axis=1)
        labels[di[int(np.argmin(m))]]=c
    return labels

def fit_strict(X60, feature_ids, k, max_rules=10, min_support=.8, min_contrast=.1, max_iter=50):
    pairs=pair_list(X60.shape[1])
    B,R=binary_forward_reverse(X60,pairs)
    labels=initial_labels(B,k)
    converged=False
    for it in range(1,max_iter+1):
        protos=build_prototypes_strict(B,R,labels,pairs,feature_ids,k,max_rules,min_support,min_contrast)
        scores=score_prototypes(X60,feature_ids,protos)
        nl=repair_empty(np.argmax(scores,axis=1).astype(int),scores,k)
        if np.array_equal(nl,labels):
            labels=nl; converged=True; break
        labels=nl
    protos=build_prototypes_strict(B,R,labels,pairs,feature_ids,k,max_rules,min_support,min_contrast)
    scores=score_prototypes(X60,feature_ids,protos)
    order=np.sort(scores,axis=1)
    margin=order[:,-1]-order[:,-2]
    return RRDirectResult(tuple(map(int,labels)),tuple(protos),tuple(map(str,feature_ids)),
                          len(pairs),it,converged,tuple(map(float,margin))),B

def q_stat(B,labels):
    labels=np.asarray(labels)
    same=[]; diff=[]
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            sim=1.0-float(np.mean(B[i]!=B[j]))
            (same if labels[i]==labels[j] else diff).append(sim)
    if not same or not diff: return float("nan")
    return float(np.mean(same)-np.mean(diff))

def forced_labels(X,feature_ids,protos):
    return np.argmax(score_prototypes(X,feature_ids,protos),axis=1).astype(int)

def rows_by_participant(map_rows, gsm_order, participants):
    gsm_to_row={g:i for i,g in enumerate(gsm_order)}
    p2=[]
    for pid in participants:
        gs=sorted([r["gsm"] for r in map_rows if r["participant_id"]==pid])
        idx=[gsm_to_row[g] for g in gs]
        if len(idx)!=2: die(f"participant {pid} does not have 2 source specimens")
        p2.append(idx)
    return np.array(p2,dtype=int)

def deterministic_subsets(pids):
    out=[]
    for r in range(20):
        s=sorted(pids,key=lambda p:hashlib.sha256(f"RPPV1_STABILITY|{r}|{p}".encode()).hexdigest())
        out.append(s[:24])
    return out

def null_permute_selected(Xsel, block_rows, k, rep, gene_ids):
    out=np.empty_like(Xsel)
    nblocks=len(block_rows)
    for gcol,gid in enumerate(gene_ids):
        rng=np.random.default_rng(np.random.SeedSequence([20260911,int(k),int(rep),int(gid)]))
        perm=rng.permutation(nblocks)
        swaps=rng.random(nblocks)<0.5
        for dest_b,src_b in enumerate(perm):
            src=block_rows[src_b].copy()
            if swaps[dest_b]: src=src[::-1]
            out[block_rows[dest_b],gcol]=Xsel[src,gcol]
    return out

def artifact_dict_with_sha(art):
    d=art.to_dict(); d["artifact_sha256"]=art.sha256(); return d

def participant_group_values(exec_rows, validation_scores, map_rows, gsm_order, index_profile):
    gsm_to_meta={r["gsm"]:r for r in map_rows}
    by={}
    for i,gsm in enumerate(gsm_order):
        ex=exec_rows[i]
        if ex.assignment!="ASSIGNED": continue
        grp="INDEX" if ex.best_profile_id==index_profile else "COMPARATOR"
        pid=gsm_to_meta[gsm]["participant_id"]
        by.setdefault((pid,grp),[]).append(float(validation_scores[i]))
    vals={"INDEX":{}, "COMPARATOR":{}}
    for (pid,g),x in by.items():
        vals[g][pid]=float(np.mean(x))
    return vals

def contrast_from_group_values(vals):
    if not vals["INDEX"] or not vals["COMPARATOR"]: return 0.0
    return float(np.mean(list(vals["INDEX"].values()))-np.mean(list(vals["COMPARATOR"].values())))

def bootstrap_contrast(vals, all_pids, n, seed):
    rng=np.random.default_rng(seed); pids=list(all_pids); arr=[]
    for _ in range(n):
        draw=rng.choice(pids,size=len(pids),replace=True)
        iv=[]; cv=[]
        for pid in draw:
            if pid in vals["INDEX"]: iv.append(vals["INDEX"][pid])
            if pid in vals["COMPARATOR"]: cv.append(vals["COMPARATOR"][pid])
        arr.append(float(np.mean(iv)-np.mean(cv)) if iv and cv else 0.0)
    return np.array(arr,float)

def relation_participant_values(X, gene_index, pair, direction, exec_rows, map_rows, gsm_order, index_profile):
    a,b=pair
    va=X[:,gene_index[a]]; vb=X[:,gene_index[b]]
    sat=(va>vb) if direction==1 else (va<vb)
    gsm_to_meta={r["gsm"]:r for r in map_rows}
    vals={"INDEX":{}, "COMPARATOR":{}}
    tmp={}
    for i,gsm in enumerate(gsm_order):
        ex=exec_rows[i]
        if ex.assignment!="ASSIGNED": continue
        grp="INDEX" if ex.best_profile_id==index_profile else "COMPARATOR"
        pid=gsm_to_meta[gsm]["participant_id"]
        tmp.setdefault((pid,grp),[]).append(float(sat[i]))
    for (pid,g),x in tmp.items(): vals[g][pid]=float(np.mean(x))
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--prep-root",required=True)
    args=ap.parse_args()
    prep=Path(args.prep_root).resolve()
    repo=prep/"repo_patch/_work_repo"
    audit=prep/"SOURCE_AUDIT_V1"
    stage1=prep/"_stage_packages_v2/pre_freeze"
    out=prep/"SOURCE_FIT_V1"; out.mkdir(exist_ok=True)

    if git(repo,"status","--porcelain"): die("repo dirty before source fit")
    head=git(repo,"rev-parse","HEAD")
    tag=git(repo,"rev-list","-n","1",EXPECTED_IMPL_TAG)
    if head!=tag: die(f"HEAD {head} is not source implementation tag commit {tag}")

    cfg=json.loads((repo/"docs/preliminary/V1_CONFIG.json").read_text(encoding="utf-8-sig"))
    impl=json.loads((repo/"docs/preliminary/SOURCE_EXECUTION_IMPLEMENTATION_V1.json").read_text(encoding="utf-8-sig"))
    cache=np.load(audit/"GSE19804_COMMON_ENTREZ_V1.npz")
    X=cache["X"].astype(float); gsms=[str(x) for x in cache["gsms"]]; genes=np.array([str(x) for x in cache["gene_ids"]])
    pools=np.array([str(x) for x in cache["pools"]])

    with (stage1/"metadata/GSE19804_participant_map_BLIND.csv").open(encoding="utf-8-sig",newline="") as f:
        map_rows=list(csv.DictReader(f))
    with (repo/"docs/preliminary/SOURCE_PARTICIPANT_SPLIT_V1.csv").open(encoding="utf-8-sig",newline="") as f:
        split=list(csv.DictReader(f))
    role={r["participant_id"]:r["role"] for r in split}
    gsm_to_row={g:i for i,g in enumerate(gsms)}
    p_to_rows={}
    for r in map_rows:
        p_to_rows.setdefault(r["participant_id"],[]).append(gsm_to_row[r["gsm"]])
    for p in p_to_rows: p_to_rows[p]=sorted(p_to_rows[p])

    assign_pids=[r["participant_id"] for r in split if r["role"]=="ASSIGNMENT_LEARNING"]
    sig_pids=[r["participant_id"] for r in split if r["role"]=="SIGNATURE_CONSTRUCTION"]
    ref_pids=[r["participant_id"] for r in split if r["role"]=="SOURCE_REFERENCE"]
    assign_rows=np.array([i for p in assign_pids for i in p_to_rows[p]],int)
    sig_rows=np.array([i for p in sig_pids for i in p_to_rows[p]],int)
    ref_rows=np.array([i for p in ref_pids for i in p_to_rows[p]],int)

    disc_mask=pools=="DISCOVERY"; val_mask=pools=="VALIDATION"
    full_top_idx=top_genes_by_mad(X,genes,disc_mask,assign_rows,60)
    full_feature_ids=genes[full_top_idx].tolist()
    X_assign60=X[assign_rows][:,full_top_idx]
    assign_block_rows=np.arange(len(assign_rows)).reshape(len(assign_pids),2)

    krows=[]; full_results={}
    full_labels={}
    for k in cfg["assignment"]["candidate_k"]:
        rec={"k":k}
        try:
            fit,B=fit_strict(X_assign60,full_feature_ids,k,
                             max_rules=cfg["assignment"]["max_rules_per_profile"],
                             min_support=cfg["assignment"]["min_rule_support"],
                             min_contrast=cfg["assignment"]["min_rule_contrast"],
                             max_iter=cfg["assignment"]["max_iter"])
            q=q_stat(B,fit.labels)
            rec.update(full_fit_ok=True,converged=bool(fit.converged),q_observed=q)
            full_results[k]=fit; full_labels[k]=np.array(fit.labels,int)
            # participant support
            support={}
            for c in range(k):
                ps=set()
                for pp,inds in zip(assign_pids,assign_block_rows):
                    if np.any(full_labels[k][inds]==c): ps.add(pp)
                support[c]=len(ps)
            rec["min_profile_participant_support_observed"]=min(support.values())
            rec["profile_participant_support"]=json.dumps(support,sort_keys=True)

            # stability
            aris=[]; ok=0
            for subp in deterministic_subsets(assign_pids):
                sub_global_rows=np.array([i for p in subp for i in p_to_rows[p]],int)
                sub_idx=top_genes_by_mad(X,genes,disc_mask,sub_global_rows,60)
                sub_ids=genes[sub_idx].tolist()
                try:
                    sf,_=fit_strict(X[sub_global_rows][:,sub_idx],sub_ids,k,
                                    max_rules=cfg["assignment"]["max_rules_per_profile"],
                                    min_support=cfg["assignment"]["min_rule_support"],
                                    min_contrast=cfg["assignment"]["min_rule_contrast"],
                                    max_iter=cfg["assignment"]["max_iter"])
                    pred=forced_labels(X[assign_rows][:,sub_idx],sub_ids,sf.prototypes)
                    aris.append(ari(full_labels[k],pred)); ok+=1
                except Exception:
                    aris.append(0.0)
            rec["stability_successes"]=ok
            rec["stability_mean_ari"]=float(np.mean(aris))

            # 199 nulls
            qnull=[]; nfail=0
            for rep in range(cfg["assignment"]["source_structure_null_replicates"]):
                Xn=null_permute_selected(X_assign60,assign_block_rows,k,rep,full_feature_ids)
                try:
                    nf,nB=fit_strict(Xn,full_feature_ids,k,
                                     max_rules=cfg["assignment"]["max_rules_per_profile"],
                                     min_support=cfg["assignment"]["min_rule_support"],
                                     min_contrast=cfg["assignment"]["min_rule_contrast"],
                                     max_iter=cfg["assignment"]["max_iter"])
                    qnull.append(q_stat(nB,nf.labels))
                except Exception:
                    qnull.append(0.0); nfail+=1
            null95=float(np.quantile(np.array(qnull),.95,method="higher"))
            rec["null_fit_failures"]=nfail
            rec["q_null95"]=null95
            rec["q_excess"]=float(q-null95)
            rec["passes"]=bool(
                fit.converged and
                ok>=18 and
                rec["stability_mean_ari"]>=cfg["assignment"]["k_min_mean_stability_ari"] and
                min(support.values())>=cfg["assignment"]["min_profile_participant_support"] and
                q>null95
            )
        except Exception as e:
            rec.update(full_fit_ok=False,passes=False,error=str(e))
        krows.append(rec)

    with (out/"source_fit_summary.csv").open("w",newline="",encoding="utf-8") as f:
        fields=sorted({k for r in krows for k in r})
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(krows)

    eligible=[r for r in krows if r.get("passes")]
    status={"schema":"SourceFitStatusV1","source":"GSE19804","software_commit":head,
            "target_values_opened":False,"target_labels_opened":False}
    if not eligible:
        status.update(status="NO_STABLE_STRUCTURE",source_reference_adequate=False)
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print("SOURCE FIT COMPLETE: NO_STABLE_STRUCTURE")
        print("Do NOT open Stage 2.")
        return

    eligible.sort(key=lambda r:(-r["q_excess"],-r["stability_mean_ari"],r["k"]))
    chosen=eligible[0]; k=int(chosen["k"]); fit=full_results[k]

    mapping_hash=sha256_file(audit/"COMMON_ENTREZ_UNIVERSE_V1.csv")
    preprocessing_hash=canonical_matrix_hash(X,gsms,genes)
    config_hash=sha256_file(repo/"docs/preliminary/V1_CONFIG.json")
    art=artifact_from_rr_direct(
        fit, artifact_id=f"GSE19804_RR_DIRECT_V1_K{k}", source_dataset="GSE19804",
        feature_namespace="Entrez Gene ID", mapping_hash=mapping_hash,
        preprocessing_hash=preprocessing_hash, config_hash=config_hash,
        software_commit=head,
        min_score=cfg["assignment"]["artifact_min_score"],
        min_margin=cfg["assignment"]["artifact_min_margin"],
        min_executable_coverage=cfg["assignment"]["artifact_min_executable_coverage"]
    )
    (out/"SOURCE_ASSIGNMENT_ARTIFACT_V1.json").write_text(
        json.dumps(artifact_dict_with_sha(art),indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

    # index profile
    pscore={}
    for p in art.prototypes:
        vals=[float(r.within_support or 0)*max(float(r.contrast or 0),0) for r in p.relations]
        pscore[p.profile_id]=float(np.mean(vals))
    index_profile=sorted(pscore,key=lambda p:(-pscore[p],p))[0]

    # execute frozen artifact on signature subset
    sig_exec=execute_rpp_artifact(X[sig_rows],genes.tolist(),art)
    sig_gsms=[gsms[i] for i in sig_rows]
    sig_map=[next(r for r in map_rows if r["gsm"]==g) for g in sig_gsms]

    # group support
    idxp=set(); compp=set()
    for ex,m in zip(sig_exec,sig_map):
        if ex.assignment!="ASSIGNED": continue
        (idxp if ex.best_profile_id==index_profile else compp).add(m["participant_id"])
    if len(idxp)<cfg["signature"]["min_group_participants_signature"] or len(compp)<cfg["signature"]["min_group_participants_signature"]:
        status.update(status="INSUFFICIENT_SOURCE_SIGNATURE",selected_k=k,index_profile=index_profile,
                      signature_index_participants=len(idxp),signature_comparator_participants=len(compp),
                      source_reference_adequate=False)
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print(f"SOURCE FIT COMPLETE: INSUFFICIENT_SOURCE_SIGNATURE (groups {len(idxp)}/{len(compp)})")
        print("Do NOT open Stage 2.")
        return

    val_idx=top_genes_by_mad(X,genes,val_mask,sig_rows,cfg["signature"]["validation_feature_budget"])
    val_genes=genes[val_idx].tolist()
    Xsig=X[sig_rows]
    gene_index={g:i for i,g in enumerate(genes)}
    pairs=[(a,b) for ii,a in enumerate(val_genes) for b in val_genes[ii+1:]]
    candidates=[]
    all_sig_pids=sorted(set(sig_pids),key=int)
    for a,b in pairs:
        best=None
        for d in (1,-1):
            vals=relation_participant_values(Xsig,gene_index,(a,b),d,sig_exec,sig_map,sig_gsms,index_profile)
            if not vals["INDEX"] or not vals["COMPARATOR"]: continue
            support=float(np.mean(list(vals["INDEX"].values())))
            contrast=contrast_from_group_values(vals)
            if best is None or contrast>best["contrast"] or (contrast==best["contrast"] and d==1):
                best={"a":a,"b":b,"direction":d,"support":support,"contrast":contrast,"vals":vals}
        if best is None: continue
        if best["support"]<cfg["signature"]["min_index_support"] or best["contrast"]<cfg["signature"]["min_observed_contrast"]:
            continue
        boot=bootstrap_contrast(best["vals"],all_sig_pids,cfg["signature"]["bootstrap_replicates"],20260913)
        pos=float(np.mean(boot>0)); med=float(np.median(boot))
        if pos<cfg["signature"]["min_positive_bootstrap_fraction"]: continue
        best["bootstrap_positive_fraction"]=pos; best["bootstrap_median_contrast"]=med
        candidates.append(best)
    candidates.sort(key=lambda z:(-z["bootstrap_median_contrast"],-z["contrast"],z["a"],z["b"],-z["direction"]))
    reuse={}; selected=[]
    for c in candidates:
        if reuse.get(c["a"],0)>=cfg["signature"]["max_feature_reuse"] or reuse.get(c["b"],0)>=cfg["signature"]["max_feature_reuse"]:
            continue
        selected.append(c); reuse[c["a"]]=reuse.get(c["a"],0)+1; reuse[c["b"]]=reuse.get(c["b"],0)+1
        if len(selected)>=cfg["signature"]["max_relations"]: break

    sig_csv=out/"source_signature.csv"
    with sig_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["feature_a","feature_b","direction","index_support","observed_contrast",
                                      "bootstrap_positive_fraction","bootstrap_median_contrast"])
        w.writeheader()
        for c in selected:
            w.writerow({"feature_a":c["a"],"feature_b":c["b"],"direction":">" if c["direction"]==1 else "<",
                        "index_support":c["support"],"observed_contrast":c["contrast"],
                        "bootstrap_positive_fraction":c["bootstrap_positive_fraction"],
                        "bootstrap_median_contrast":c["bootstrap_median_contrast"]})

    if len(selected)<cfg["signature"]["min_relations"]:
        status.update(status="INSUFFICIENT_SOURCE_SIGNATURE",selected_k=k,index_profile=index_profile,
                      selected_signature_relations=len(selected),source_reference_adequate=False)
        (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        print(f"SOURCE FIT COMPLETE: INSUFFICIENT_SOURCE_SIGNATURE ({len(selected)} relations)")
        print("Do NOT open Stage 2.")
        return

    sig_obj={"schema":"ValidationSignatureV1","source_dataset":"GSE19804","index_profile":index_profile,
             "software_commit":head,"relations":[
                 {"feature_a":c["a"],"feature_b":c["b"],"direction":">" if c["direction"]==1 else "<","weight":1.0}
                 for c in selected]}
    raw=(json.dumps(sig_obj,sort_keys=True,separators=(",",":"),ensure_ascii=False)+"\n").encode()
    sig_obj["sha256"]=hashlib.sha256(raw).hexdigest()
    (out/"VALIDATION_SIGNATURE_V1.json").write_text(json.dumps(sig_obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

    # reference execution
    ref_exec=execute_rpp_artifact(X[ref_rows],genes.tolist(),art)
    ref_gsms=[gsms[i] for i in ref_rows]
    ref_map=[next(r for r in map_rows if r["gsm"]==g) for g in ref_gsms]
    vs=[]
    for gi in ref_rows:
        sats=[]
        for c in selected:
            aa=gene_index[c["a"]]; bb=gene_index[c["b"]]
            xa=X[gi,aa]; xb=X[gi,bb]
            sats.append(float(xa>xb) if c["direction"]==1 else float(xa<xb))
        vs.append(float(np.mean(sats)))
    vals=participant_group_values(ref_exec,vs,ref_map,ref_gsms,index_profile)
    C=contrast_from_group_values(vals)
    boot=bootstrap_contrast(vals,sorted(set(ref_pids),key=int),cfg["source_reference"]["bootstrap_replicates"],20260914)
    lower=float(np.quantile(boot,cfg["source_reference"]["source_one_sided_alpha"],method="lower"))
    specimen_cov=float(np.mean([e.assignment=="ASSIGNED" for e in ref_exec]))
    covered_p=set(m["participant_id"] for e,m in zip(ref_exec,ref_map) if e.assignment=="ASSIGNED")
    part_cov=len(covered_p)/len(ref_pids)
    ni=len(vals["INDEX"]); nc=len(vals["COMPARATOR"])
    adequate=bool(specimen_cov>=cfg["source_reference"]["min_specimen_assignment_coverage"] and
                  part_cov>=cfg["source_reference"]["min_participant_assignment_coverage"] and
                  ni>=cfg["source_reference"]["min_group_participants"] and
                  nc>=cfg["source_reference"]["min_group_participants"] and
                  C>=cfg["source_reference"]["source_contrast_floor"] and lower>0)

    ref_summary={"selected_k":k,"index_profile":index_profile,"C_source":C,"one_sided_95_lower":lower,
                 "specimen_assignment_coverage":specimen_cov,"participant_assignment_coverage":part_cov,
                 "index_participants":ni,"comparator_participants":nc,"source_reference_adequate":adequate}
    with (out/"source_reference.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(ref_summary)); w.writeheader(); w.writerow(ref_summary)

    # source assignments, no phenotype labels
    with (out/"source_reference_assignments.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["gsm","participant_id","assignment","best_profile_id","best_score","margin","executable_coverage","validation_score"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for g,m,e,v in zip(ref_gsms,ref_map,ref_exec,vs):
            w.writerow({"gsm":g,"participant_id":m["participant_id"],"assignment":e.assignment,
                        "best_profile_id":e.best_profile_id,"best_score":e.best_score,"margin":e.margin,
                        "executable_coverage":e.executable_coverage,"validation_score":v})

    status.update(status="SOURCE_REFERENCE_PASS" if adequate else "INSUFFICIENT_SOURCE_REFERENCE",
                  selected_k=k,index_profile=index_profile,artifact_sha256=art.sha256(),
                  signature_sha256=sig_obj["sha256"],selected_signature_relations=len(selected),
                  C_source=C,source_reference_lower95=lower,source_reference_adequate=adequate)
    (out/"SOURCE_RESULT_STATUS.json").write_text(json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

    # manifest
    with (out/"SOURCE_RESULTS_SHA256.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["filename","bytes","sha256"])
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name!="SOURCE_RESULTS_SHA256.csv":
                w.writerow([p.name,p.stat().st_size,sha256_file(p)])

    print("SOURCE FIT COMPLETE")
    print(f"Selected K: {k}")
    print(f"Index profile: {index_profile}")
    print(f"Signature relations: {len(selected)}")
    print(f"C_source: {C:.6f}")
    print(f"One-sided 95% lower bound: {lower:.6f}")
    print(f"Reference support: index={ni}, comparator={nc}")
    print(f"Assignment coverage: specimen={specimen_cov:.3f}, participant={part_cov:.3f}")
    print(f"STATUS: {status['status']}")
    print(f"Output: {out}")
    print("TARGET VALUES OPENED: NO")
    print("TARGET LABELS OPENED: NO")
    if not adequate:
        print("Do NOT open Stage 2.")
    else:
        print("Do NOT open Stage 2 yet; freeze the source artifact/signature first.")

if __name__=="__main__":
    main()
