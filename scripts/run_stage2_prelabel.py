#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,gzip,hashlib,importlib.util,json,math,os,re,subprocess,sys,zipfile
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

from relational_patient_profiles.artifact import (
    RPPRelation,RPPPrototype,RelationalPatientProfileArtifact,execute_rpp_artifact
)

EXPECTED_IMPL_TAG="preliminary-stage2-partial-recovery-2026-09-07"
EXPECTED_STAGE2_SHA="94e648dfa9ec9b9cd16a6eaea067f1d823e7d74eb128cc27cfc6de473081e2df"
EXPECTED_ARTIFACT_SHA="9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b"
EXPECTED_SIGNATURE_SHA="64eb3631d6a25f0fdf0deede0e0c2da17c63fb25d82ee4612aaa8bd12470128d"
TARGETS={"GSE27262":"GPL570","GSE32863":"GPL6884"}

def die(x): raise SystemExit("ERROR: "+str(x))
def sha(p:Path):
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
    s=importlib.util.spec_from_file_location("v1src",p)
    m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

def load_artifact(path):
    d=json.loads(path.read_text(encoding="utf-8-sig"))
    ps=[]
    for p in d["prototypes"]:
        rs=[]
        for r in p["relations"]:
            rs.append(RPPRelation(str(r["feature_a"]),str(r["feature_b"]),
                1 if r["direction"]==">" else -1,float(r["weight"]),
                r.get("within_support"),r.get("contrast")))
        ps.append(RPPPrototype(p["profile_id"],tuple(rs)))
    art=RelationalPatientProfileArtifact(
        artifact_id=d["artifact_id"],source_dataset=d["source_dataset"],
        feature_namespace=d["feature_namespace"],mapping_hash=d["mapping_hash"],
        preprocessing_hash=d["preprocessing_hash"],config_hash=d["config_hash"],
        software_commit=d["software_commit"],prototypes=tuple(ps),
        min_score=float(d["min_score"]),min_margin=float(d["min_margin"]),
        min_executable_coverage=float(d["min_executable_coverage"]),
        created_prelabel=bool(d["created_prelabel"]),
        transportability_certificate=d.get("transportability_certificate"))
    if art.sha256()!=EXPECTED_ARTIFACT_SHA: die("artifact SHA mismatch")
    return art

def parse_annotation(path):
    p2g={}
    with gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="") as f:
        reader=None
        for line in f:
            if line.startswith("ID\t"):
                header=next(csv.reader([line],delimiter="\t"))
                reader=csv.DictReader(f,fieldnames=header,delimiter="\t"); break
        if reader is None or "Gene ID" not in header: die(f"bad annotation {path}")
        for r in reader:
            probe=(r.get("ID") or "").strip(); gid=(r.get("Gene ID") or "").strip()
            if probe and re.fullmatch(r"\d+",gid): p2g[probe]=gid
    return p2g

def safe_extract(zp:Path,dest:Path):
    dest.mkdir(exist_ok=True)
    with zipfile.ZipFile(zp) as z:
        for info in z.infolist():
            target=(dest/info.filename).resolve()
            if not str(target).startswith(str(dest.resolve())): die("unsafe zip member")
        z.extractall(dest)

def verify_extracted_matches_zip(zp:Path,dest:Path):
    """Verify an already-opened Stage-2 extraction against the exact pre-hashed ZIP."""
    with zipfile.ZipFile(zp) as z:
        expected=set()
        for info in z.infolist():
            if info.is_dir():
                continue
            expected.add(info.filename.replace("\\","/"))
            p=(dest/info.filename)
            if not p.is_file():
                die(f"missing previously extracted Stage-2 member: {info.filename}")
            h1=hashlib.sha256()
            with z.open(info,"r") as f:
                for b in iter(lambda:f.read(1024*1024),b""):
                    h1.update(b)
            if h1.hexdigest()!=sha(p):
                die(f"previously extracted Stage-2 member differs from frozen ZIP: {info.filename}")
        actual=set()
        for p in dest.rglob("*"):
            if p.is_file():
                actual.add(p.relative_to(dest).as_posix())
        extra=actual-expected
        if extra:
            die("unexpected extra files in Stage-2 extraction: "+str(sorted(extra)[:5]))

def locate_target_matrix(root:Path,dataset:str):
    cand=[]
    for p in root.rglob("*"):
        if not p.is_file(): continue
        n=p.name.lower()
        if dataset.lower() in n and (n.endswith(".gz") or n.endswith(".tsv") or n.endswith(".txt")):
            if "expression" in n or "matrix" in n: cand.append(p)
    if not cand: die(f"target matrix not found for {dataset}")
    cand.sort(key=lambda p:(0 if "expression_only" in p.name.lower() else 1,len(str(p))))
    return cand[0]

def map_expression(path,p2g,common_genes):
    opener=gzip.open if path.suffix.lower()==".gz" else open
    by={}
    with opener(path,"rt",encoding="utf-8",errors="replace",newline="") as f:
        rd=csv.reader(f,delimiter="\t"); header=next(rd)
        if header[0].strip().strip('"')!="ID_REF": die(f"bad expression header {path}")
        gsms=[x.strip().strip('"') for x in header[1:]]
        for row in rd:
            if not row: continue
            probe=row[0].strip().strip('"'); gid=p2g.get(probe)
            if gid is None or gid not in common_genes: continue
            vals=np.array([float(x) if x not in ("","NA","NaN","nan") else np.nan for x in row[1:]],float)
            by.setdefault(gid,[]).append(vals)
    gi={g:i for i,g in enumerate(common_genes)}
    X=np.full((len(gsms),len(common_genes)),np.nan,float)
    for g,arrs in by.items():
        X[:,gi[g]]=np.nanmedian(np.vstack(arrs),axis=0)
    order=np.argsort(np.array(gsms,dtype=str))
    return X[order], [gsms[i] for i in order], len(by)

def strict_signature_scores(X,genes,sig):
    idx={str(g):i for i,g in enumerate(genes)}
    out=np.full(len(X),np.nan,float)
    for i,row in enumerate(X):
        vals=[]; ok=True
        for r in sig["relations"]:
            a,b=str(r["feature_a"]),str(r["feature_b"])
            if a not in idx or b not in idx: ok=False; break
            xa,xb=row[idx[a]],row[idx[b]]
            if not(np.isfinite(xa) and np.isfinite(xb)): ok=False; break
            vals.append(float(xa>xb) if r["direction"]==">" else float(xa<xb))
        if ok: out[i]=float(np.mean(vals))
    return out

def exec_fingerprint(execs):
    # Stable comparison for C-invariant assignments.
    return [(e.assignment,e.reason,e.best_profile_id,e.best_score,e.margin,
             e.executable_coverage,tuple(e.profile_scores),tuple(e.profile_coverages)) for e in execs]

def participant_sort_key(x):
    """Canonical deterministic order for opaque participant identifiers.

    Purely numeric source IDs retain the historical numeric ordering used by V3.
    Non-numeric target IDs are ordered lexically and are never coerced to int.
    """
    s=str(x)
    return (0,int(s),s) if s.isdigit() else (1,s.casefold(),s)

def group_arrays(execs,vs,maprows,gsms,index_profile):
    meta={r["gsm"]:r for r in maprows}
    pids=sorted({r["participant_id"] for r in maprows},key=participant_sort_key)
    pi={p:i for i,p in enumerate(pids)}
    acc_i={p:[] for p in pids}; acc_c={p:[] for p in pids}
    assigned_spec=0; eval_assigned=0; covered=set()
    for e,v,g in zip(execs,vs,gsms):
        p=meta[g]["participant_id"]
        if e.assignment!="ASSIGNED": continue
        assigned_spec+=1; covered.add(p)
        if not np.isfinite(v): continue
        eval_assigned+=1
        if e.best_profile_id==index_profile: acc_i[p].append(float(v))
        else: acc_c[p].append(float(v))
    ai=np.full(len(pids),np.nan); ac=np.full(len(pids),np.nan)
    for p in pids:
        if acc_i[p]: ai[pi[p]]=float(np.mean(acc_i[p]))
        if acc_c[p]: ac[pi[p]]=float(np.mean(acc_c[p]))
    return pids,ai,ac,{
        "specimen_assignment_coverage":assigned_spec/len(execs) if execs else 0.0,
        "participant_assignment_coverage":len(covered)/len(pids) if pids else 0.0,
        "validation_evaluability":eval_assigned/assigned_spec if assigned_spec else 0.0,
        "index_participants":int(np.isfinite(ai).sum()),
        "comparator_participants":int(np.isfinite(ac).sum())
    }

def contrast(ai,ac):
    return float(np.nanmean(ai)-np.nanmean(ac)) if np.isfinite(ai).any() and np.isfinite(ac).any() else float("nan")

def bootstrap_indices(n,dataset_num,B=2000):
    rng=np.random.default_rng(np.random.SeedSequence([20260915,int(dataset_num)]))
    return rng.choice(np.arange(n),size=(B,n),replace=True)

def draw_contrast(ai,ac,inds):
    A=ai[inds]; C=ac[inds]
    ca=np.isfinite(A).sum(1); cc=np.isfinite(C).sum(1)
    sa=np.nansum(A,1); sc=np.nansum(C,1)
    out=np.zeros(len(inds),float)
    ok=(ca>0)&(cc>0)
    out[ok]=sa[ok]/ca[ok]-sc[ok]/cc[ok]
    return out

def structural_decision(ai,ac,support,source_boot,target_inds,cfg,Csource):
    Ct=contrast(ai,ac)
    adequate=(
        support["specimen_assignment_coverage"]>=cfg["retention"]["min_specimen_assignment_coverage"] and
        support["participant_assignment_coverage"]>=cfg["retention"]["min_participant_assignment_coverage"] and
        support["index_participants"]>=cfg["retention"]["min_group_participants"] and
        support["comparator_participants"]>=cfg["retention"]["min_group_participants"] and
        support["validation_evaluability"]>=cfg["retention"]["min_validation_evaluability"]
    )
    td=draw_contrast(ai,ac,target_inds)
    Dd=td-cfg["retention"]["rho"]*source_boot
    loD=float(np.quantile(Dd,.05,method="lower"))
    hiD=float(np.quantile(Dd,.95,method="higher"))
    loC=float(np.quantile(td,.05,method="lower"))
    hiC=float(np.quantile(td,.95,method="higher"))
    D=Ct-cfg["retention"]["rho"]*Csource if np.isfinite(Ct) else float("nan")
    if not adequate:
        state="INSUFFICIENT_SUPPORT"
    elif Ct>=cfg["retention"]["target_absolute_contrast_floor"] and loD>0:
        state="PASS"
    elif hiC<=0 or hiD<0:
        state="FAIL"
    else:
        state="INSUFFICIENT_SUPPORT"
    return {"C_target":Ct,"D":D,"R":Ct/Csource if np.isfinite(Ct) and Csource>0 else float("nan"),
            "C_lower95_one_sided":loC,"C_upper95_one_sided":hiC,
            "D_lower95_one_sided":loD,"D_upper95_one_sided":hiD,
            "decision":state,**support}

def participant_blocks(maprows,gsms):
    meta={r["gsm"]:r for r in maprows}; pos={g:i for i,g in enumerate(gsms)}
    by={}
    for r in maprows: by.setdefault(r["participant_id"],[]).append(r["gsm"])
    classes={}
    for pid,gg in by.items():
        gg=sorted(gg)
        vals=sorted({meta[g]["pair_status"] for g in gg})
        cls="|".join(vals)+f"|n{len(gg)}"
        classes.setdefault(cls,[]).append(np.array([pos[g] for g in gg],int))
    return classes

def randomise_genes(X,genes,maprows,gsms,selected,seed,dataset_num,rep):
    Y=X.copy(); gi={g:i for i,g in enumerate(genes)}; classes=participant_blocks(maprows,gsms)
    for gid in selected:
        j=gi[gid]
        rng=np.random.default_rng(np.random.SeedSequence([int(seed),int(dataset_num),int(rep),int(gid)]))
        for cls,blocks in sorted(classes.items()):
            nb=len(blocks)
            perm=rng.permutation(nb); swaps=rng.random(nb)<0.5
            for d,s in enumerate(perm):
                src=blocks[s].copy()
                if len(src)==2 and swaps[d]: src=src[::-1]
                Y[blocks[d],j]=X[src,j]
    return Y

def relation_states(X,genes,relations):
    gi={g:i for i,g in enumerate(genes)}; states=[]
    for r in relations:
        a,b=str(r["feature_a"]),str(r["feature_b"])
        xa,xb=X[:,gi[a]],X[:,gi[b]]
        if r["direction"]==">": z=xa>xb
        else: z=xa<xb
        states.append(z)
    return np.column_stack(states)

def kmedoids(Z,gsms,K=2,max_iter=100):
    Z=np.asarray(Z,float)
    # finite-only diagnostic input is required for the frozen 60 common genes.
    if not np.isfinite(Z).all(): die("nonfinite value in k-medoids diagnostic")
    G=np.sum(Z*Z,axis=1)
    D=np.maximum(G[:,None]+G[None,:]-2*Z@Z.T,0.0)
    med=[0]
    while len(med)<K:
        d=np.min(D[:,med],axis=1); d[med]=-1
        med.append(int(np.argmax(d)))
    med=np.array(med,int)
    lab=None
    for _ in range(max_iter):
        dist=D[:,med]; nl=np.argmin(dist,axis=1)
        nmed=med.copy()
        for c in range(K):
            mem=np.flatnonzero(nl==c)
            if len(mem)==0:
                cand=[i for i in range(len(Z)) if i not in nmed]
                nmed[c]=max(cand,key=lambda i:(np.min(D[i,nmed]),-i))
            else:
                costs=D[np.ix_(mem,mem)].sum(1)
                nmed[c]=int(mem[np.argmin(costs)])
        if np.array_equal(nmed,med):
            lab=nl; break
        med=nmed
    if lab is None: lab=np.argmin(D[:,med],axis=1)
    return lab.astype(int),med.tolist()

def ari(a,b):
    a=np.asarray(a); b=np.asarray(b); n=len(a)
    ua,ia=np.unique(a,return_inverse=True); ub,ib=np.unique(b,return_inverse=True)
    C=np.zeros((len(ua),len(ub)),dtype=np.int64); np.add.at(C,(ia,ib),1)
    c2=lambda x:x*(x-1)//2
    sr=int(np.sum(c2(C.sum(1)))); sc=int(np.sum(c2(C.sum(0)))); s=int(np.sum(c2(C))); t=c2(n)
    exp=sr*sc/t if t else 0.; mx=.5*(sr+sc); den=mx-exp
    return 1.0 if den==0 else float((s-exp)/den)

def flatten_exec(dataset,execs,vs,maprows,gsms,outpath):
    meta={r["gsm"]:r for r in maprows}
    with outpath.open("w",newline="",encoding="utf-8") as f:
        fields=["dataset","gsm","participant_id","pair_status","assignment","reason","best_profile_id",
                "best_score","margin","executable_coverage","profile_scores","profile_coverages","validation_score"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for g,e,v in zip(gsms,execs,vs):
            w.writerow({"dataset":dataset,"gsm":g,"participant_id":meta[g]["participant_id"],
                "pair_status":meta[g]["pair_status"],"assignment":e.assignment,"reason":e.reason,
                "best_profile_id":e.best_profile_id,"best_score":e.best_score,"margin":e.margin,
                "executable_coverage":e.executable_coverage,
                "profile_scores":json.dumps(e.profile_scores),"profile_coverages":json.dumps(e.profile_coverages),
                "validation_score":v})

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--prep-root",required=True);a=ap.parse_args()
    root=Path(a.prep_root).resolve(); repo=root/"repo_patch/_work_repo"
    head=git(repo,"rev-parse","HEAD"); tag=git(repo,"rev-list","-n","1",EXPECTED_IMPL_TAG)
    if head!=tag: die("HEAD is not exact Stage-2 implementation-freeze commit")
    if git(repo,"status","--porcelain"): die("repo dirty at Stage-2 execution start")

    # Locate + verify exact sealed Stage-2 ZIP before opening.
    name="POST_PROTOCOL_TARGET_VALUES_2026-09-07.zip"
    candidates=[root/name,root.parent/name]
    zp=next((p for p in candidates if p.exists()),None)
    if zp is None: die(f"{name} not found in prep root or parent")
    if sha(zp)!=EXPECTED_STAGE2_SHA: die("Stage-2 ZIP SHA-256 mismatch")

    opened=root/"STAGE2_OPENED_TARGET_VALUES"
    results=root/"STAGE2_PRELABEL_RESULTS"; results.mkdir(exist_ok=True)
    access_path=results/"STAGE2_ACCESS_RECORD.json"

    if opened.exists():
        # Stage 2 was already opened by the first, failed prelabel run. Preserve that
        # access event; verify bytes rather than deleting/re-extracting the payload.
        if not access_path.is_file():
            die("Stage-2 extraction exists but the original access record is missing")
        prior=json.loads(access_path.read_text(encoding="utf-8-sig"))
        if prior.get("stage2_sha256")!=EXPECTED_STAGE2_SHA:
            die("original Stage-2 access record SHA mismatch")
        if prior.get("stage3_labels_opened") is not False:
            die("original access record does not keep Stage-3 labels sealed")
        verify_extracted_matches_zip(zp,opened)

        # The original crash occurred before authentic_transfer / assignment output.
        # Refuse to overwrite an existing endpoint result silently.
        endpoint_files=[
            "authentic_transfer.csv","pilot_summary.csv","conditional_null.csv",
            "false_reassurance_summary.csv","GSE27262_prelabel_assignments.csv",
            "GSE32863_prelabel_assignments.csv","STAGE2_PRELABEL_STATUS.json"
        ]
        present=[x for x in endpoint_files if (results/x).exists()]
        if present:
            die("pre-existing Stage-2 endpoint outputs found; refusing silent rerun: "+str(present))

        resume={"schema":"Stage2ResumeRecord/v1",
                "resumed_at_utc":datetime.now(timezone.utc).isoformat(),
                "original_access_record":"STAGE2_ACCESS_RECORD.json",
                "stage2_sha256":EXPECTED_STAGE2_SHA,
                "repair_commit":head,
                "reason":"resume after archived partial Stage-2 outputs and opaque-ID runner repair; scientific artifact/signature/thresholds unchanged",
                "stage3_labels_opened":False,"colorectal_opened":False}
        (results/"STAGE2_RESUME_RECORD.json").write_text(json.dumps(resume,indent=2)+"\n")
    else:
        safe_extract(zp,opened)
        access={"schema":"Stage2AccessRecord/v1","opened_at_utc":datetime.now(timezone.utc).isoformat(),
                "stage2_zip":str(zp),"stage2_sha256":sha(zp),"implementation_commit":head,
                "stage3_labels_opened":False,"colorectal_opened":False}
        access_path.write_text(json.dumps(access,indent=2)+"\n")

    # Fail if payload unexpectedly contains label-like files.
    bad=[p for p in opened.rglob("*") if p.is_file() and any(x in p.name.lower() for x in ("label","evaluation_label","phenotype"))]
    if bad: die("label-like file found inside Stage-2 payload: "+str(bad[:3]))

    v1=load_v1(repo)
    cfg=json.loads((repo/"docs/preliminary/V1_CONFIG.json").read_text(encoding="utf-8-sig"))
    art=load_artifact(repo/"docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json")
    sig=json.loads((repo/"docs/preliminary/v3_source_results/VALIDATION_SIGNATURE_V3.json").read_text())
    if sig["sha256"]!=EXPECTED_SIGNATURE_SHA: die("signature SHA mismatch")
    source_status=json.loads((repo/"docs/preliminary/v3_source_results/SOURCE_RESULT_STATUS.json").read_text())
    Csource=float(source_status["C_source"])

    cache=np.load(root/"SOURCE_AUDIT_V1/GSE19804_COMMON_ENTREZ_V1.npz")
    Xsrc=cache["X"].astype(float); src_gsms=[str(x) for x in cache["gsms"]]
    genes=[str(x) for x in cache["gene_ids"]]; pools=np.array([str(x) for x in cache["pools"]])
    gene_set=set(genes)

    # Source metadata, frozen V3 reference, exact source bootstrap reconstruction.
    with (root/"_stage_packages_v2/pre_freeze/metadata/GSE19804_participant_map_BLIND.csv").open(encoding="utf-8-sig",newline="") as f:
        smap=list(csv.DictReader(f))
    with (repo/"docs/preliminary/v3_source_split/SOURCE_PARTICIPANT_SPLIT_V3.csv").open(encoding="utf-8-sig",newline="") as f:
        v3split=list(csv.DictReader(f))
    refp=[r["participant_id"] for r in v3split if r["role"]=="SOURCE_REFERENCE_V3"]
    sg2i={g:i for i,g in enumerate(src_gsms)}; p2={}
    for r in smap:p2.setdefault(r["participant_id"],[]).append(sg2i[r["gsm"]])
    for p in p2:p2[p]=sorted(p2[p])
    refrows=np.array([i for p in refp for i in p2[p]],int)
    refgs=[src_gsms[i] for i in refrows]; smeta={r["gsm"]:r for r in smap}; refmap=[smeta[g] for g in refgs]
    refexec=execute_rpp_artifact(Xsrc[refrows],genes,art)
    refvs=strict_signature_scores(Xsrc[refrows],genes,sig)
    _,sai,sac,ssup=group_arrays(refexec,refvs,refmap,refgs,sig["index_profile"])
    # Rebuild exact V3 group dict then exact V1 bootstrap implementation.
    vals={"INDEX":{},"COMPARATOR":{}}
    spids=sorted(refp,key=int)
    for p,x,y in zip(spids,sai,sac):
        if np.isfinite(x):vals["INDEX"][p]=float(x)
        if np.isfinite(y):vals["COMPARATOR"][p]=float(y)
    source_boot=v1.bootstrap_contrast(vals,spids,cfg["source_reference"]["bootstrap_replicates"],20260914)
    source_lower=float(np.quantile(source_boot,cfg["source_reference"]["source_one_sided_alpha"],method="lower"))
    if abs(v1.contrast_from_group_values(vals)-Csource)>1e-12: die("C_source reconstruction mismatch")
    if abs(source_lower-float(source_status["source_reference_lower95"]))>1e-12: die("source lower95 reconstruction mismatch")

    # Source-frozen top60 and scaling for B diagnostic.
    with (repo/"docs/preliminary/SOURCE_PARTICIPANT_SPLIT_V1.csv").open(encoding="utf-8-sig",newline="") as f:
        osplit=list(csv.DictReader(f))
    assignp=[r["participant_id"] for r in osplit if r["role"]=="ASSIGNMENT_LEARNING"]
    assignrows=np.array([i for p in assignp for i in p2[p]],int)
    top60=v1.top_genes_by_mad(Xsrc,np.array(genes),pools=="DISCOVERY",assignrows,60)
    top60genes=[genes[i] for i in top60]
    src_mu=np.mean(Xsrc[assignrows][:,top60],axis=0); src_sd=np.std(Xsrc[assignrows][:,top60],axis=0,ddof=0)
    src_sd=np.where(src_sd>0,src_sd,1.0)
    s_source=float(np.median(np.std(Xsrc,axis=1,ddof=0)))

    # frozen relation sets
    artifact_rel=[]
    for p in art.prototypes:
        for r in p.relations:
            artifact_rel.append({"feature_a":r.feature_a,"feature_b":r.feature_b,"direction":">" if r.direction==1 else "<"})
    all_rel=artifact_rel+sig["relations"]
    coregenes=sorted({x for r in artifact_rel for x in (str(r["feature_a"]),str(r["feature_b"]))},
                     key=lambda g:hashlib.sha256(f"D_DROPOUT|{g}".encode()).hexdigest())
    siggenes=sorted({x for r in sig["relations"] for x in (str(r["feature_a"]),str(r["feature_b"]))},
                    key=lambda g:hashlib.sha256(f"C_DAMAGE|{g}".encode()).hexdigest())
    if set(coregenes)&set(siggenes): die("assignment core and validation signature genes overlap")

    # annotations
    ann570=parse_annotation(root/"_stage_packages_v2/pre_freeze/annotations/GPL570.annot.gz")
    ann6884=parse_annotation(root/"_stage_packages_v2/pre_freeze/annotations/GPL6884.annot.gz")
    anns={"GPL570":ann570,"GPL6884":ann6884}

    authentic=[]; Brows=[]; Crows=[]; Drows=[]; nullrows=[]; false_rows=[]
    target_cache={}
    for dataset,platform in TARGETS.items():
        dnum=int(dataset.replace("GSE",""))
        matrix=locate_target_matrix(opened,dataset)
        X,gsms,nmapped=map_expression(matrix,anns[platform],genes)
        mpath=root/f"_stage_packages_v2/pre_freeze/metadata/{dataset}_participant_map_BLIND.csv"
        with mpath.open(encoding="utf-8-sig",newline="") as f: maprows=list(csv.DictReader(f))
        if set(gsms)!={r["gsm"] for r in maprows}: die(f"{dataset} GSM mismatch target matrix vs blind map")
        # Reorder map is not needed; group funcs use GSM key.
        ex=execute_rpp_artifact(X,genes,art); vs=strict_signature_scores(X,genes,sig)
        pids,ai,ac,sup=group_arrays(ex,vs,maprows,gsms,sig["index_profile"])
        tinds=bootstrap_indices(len(pids),dnum,2000)
        dec=structural_decision(ai,ac,sup,source_boot,tinds,cfg,Csource)
        dec.update({"dataset":dataset,"platform":platform,"mapped_common_genes":nmapped,
                    "n_specimens":len(gsms),"n_participants":len(pids)})
        authentic.append(dec)
        flatten_exec(dataset,ex,vs,maprows,gsms,results/f"{dataset}_prelabel_assignments.csv")
        target_cache[dataset]=(X,gsms,maprows,ex,vs,pids,tinds)

        # B
        # Frozen B specification: initialise PCG64(20260916) independently for
        # each target in lexical GSM order. No dataset-specific B substream.
        z_rng=np.random.default_rng(20260916)
        z=z_rng.standard_normal(len(gsms)); u=z_rng.standard_normal(len(gsms))
        original_states=relation_states(X,genes,all_rel)
        gi={g:i for i,g in enumerate(genes)}
        t60=[gi[g] for g in top60genes]
        Z0=(X[:,t60]-src_mu)/src_sd
        lab0,med0=kmedoids(Z0,gsms,2)
        accepted=np.array([e.assignment=="ASSIGNED" for e in ex])
        rpplab=np.array([e.best_profile_id if e.assignment=="ASSIGNED" else "UNASSIGNED" for e in ex],dtype=object)
        for t in (0,0.25,0.5,1,2):
            Y=np.exp(t*z)[:,None]*X + (t*s_source*u)[:,None]
            st=relation_states(Y,genes,all_rel)
            flips=int(np.sum(st!=original_states))
            if flips!=0: die(f"{dataset} B relation flips at t={t}: {flips}")
            exb=execute_rpp_artifact(Y,genes,art)
            if exec_fingerprint(exb)!=exec_fingerprint(ex): die(f"{dataset} B assignment changed at t={t}")
            Z=(Y[:,t60]-src_mu)/src_sd
            lab,med=kmedoids(Z,gsms,2)
            Brows.append({"dataset":dataset,"t":t,"strict_relation_flips":flips,
                "strict_relation_flip_rate":0.0,
                "accepted_rpp_vs_value_cluster_ari":ari(rpplab[accepted],lab[accepted]) if accepted.sum()>1 else float("nan"),
                "value_cluster_intact_vs_perturbed_ari":ari(lab0,lab),
                "accepted_fraction":float(accepted.mean()),"value_medoids":json.dumps(med)})

        # C nested damage
        for dose in (0,0.1,0.25,0.5,0.75,1.0):
            reps=[0] if dose==0 else range(20)
            nsel=int(math.ceil(dose*len(siggenes)))
            selected=siggenes[:nsel]
            for rep in reps:
                Y=X if dose==0 else randomise_genes(X,genes,maprows,gsms,selected,20260917,dnum,rep)
                exc=execute_rpp_artifact(Y,genes,art)
                if exec_fingerprint(exc)!=exec_fingerprint(ex): die(f"{dataset} C assignment invariant failed dose={dose} rep={rep}")
                vsc=strict_signature_scores(Y,genes,sig)
                _,ci,cc,csup=group_arrays(exc,vsc,maprows,gsms,sig["index_profile"])
                cd=structural_decision(ci,cc,csup,source_boot,tinds,cfg,Csource)
                Crows.append({"dataset":dataset,"dose":dose,"replicate":rep,"n_disrupted_signature_genes":nsel,**cd})

        # D dropout
        prev_profile_cov=None
        for dose in (0,0.1,0.25,0.5,0.75,1.0):
            nrm=int(math.ceil(dose*len(coregenes))); removed=set(coregenes[:nrm])
            keep=[i for i,g in enumerate(genes) if g not in removed]
            kd=[genes[i] for i in keep]
            exd=execute_rpp_artifact(X[:,keep],kd,art)
            if dose==1.0 and any(e.assignment=="ASSIGNED" for e in exd): die(f"{dataset} D full dropout did not unassign all")
            profmean=np.mean(np.array([e.profile_coverages for e in exd],float),axis=0)
            if prev_profile_cov is not None and np.any(profmean>prev_profile_cov+1e-12): die(f"{dataset} D profile coverage increased")
            prev_profile_cov=profmean
            _,di,dc,dsup=group_arrays(exd,vs,maprows,gsms,sig["index_profile"])
            dd=structural_decision(di,dc,dsup,source_boot,tinds,cfg,Csource)
            reasons={}
            for e in exd: reasons[e.reason]=reasons.get(e.reason,0)+1
            Drows.append({"dataset":dataset,"dose":dose,"n_removed_core_genes":nrm,
                "mean_profile_coverages":json.dumps([float(x) for x in profmean]),
                "assigned_fraction":float(np.mean([e.assignment=="ASSIGNED" for e in exd])),
                "unassigned_reasons":json.dumps(reasons,sort_keys=True),**dd})

        # Conditional null - randomize all unique signature genes.
        passn=0
        for rep in range(999):
            Y=randomise_genes(X,genes,maprows,gsms,siggenes,20260918,dnum,rep)
            vn=strict_signature_scores(Y,genes,sig)
            _,ni,nc,nsup=group_arrays(ex,vn,maprows,gsms,sig["index_profile"])
            nd=structural_decision(ni,nc,nsup,source_boot,tinds,cfg,Csource)
            ispass=nd["decision"]=="PASS"; passn+=int(ispass)
            nullrows.append({"dataset":dataset,"replicate":rep,**nd})
        bt=binomtest(passn,999)
        ci=bt.proportion_ci(confidence_level=.95,method="exact")
        false_rows.append({"dataset":dataset,"passes":passn,"replicates":999,
            "false_reassurance":passn/999,"exact95_lower":float(ci.low),"exact95_upper":float(ci.high),
            "calibration_acceptable":bool(ci.high<=0.10)})

    # Write outputs.
    def write_rows(name,rows):
        p=results/name
        keys=[]
        for r in rows:
            for k in r:
                if k not in keys: keys.append(k)
        with p.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

    write_rows("authentic_transfer.csv",authentic)
    write_rows("perturbation_B.csv",Brows)
    write_rows("perturbation_C.csv",Crows)
    write_rows("perturbation_D.csv",Drows)
    write_rows("conditional_null.csv",nullrows)
    write_rows("false_reassurance_summary.csv",false_rows)
    write_rows("baseline_summary.csv",[{"status":"DEFERRED_NOT_RUN","reason":"non-gating explanatory control; primary endpoint unchanged"}])
    write_rows("source_full_pipeline_null_summary.csv",[{"status":"DEFERRED_NOT_RUN","reason":"non-gating source diagnostic; cannot rescue target result"}])

    pilot=[]
    for r in authentic:
        fr=next(x for x in false_rows if x["dataset"]==r["dataset"])
        pilot.append({"dataset":r["dataset"],"authentic_decision":r["decision"],"C_target":r["C_target"],
                      "D":r["D"],"R":r["R"],"false_reassurance":fr["false_reassurance"],
                      "false_reassurance_upper95":fr["exact95_upper"],
                      "calibration_acceptable":fr["calibration_acceptable"]})
    write_rows("pilot_summary.csv",pilot)

    manifest=results/"STAGE2_PRELABEL_SHA256.csv"
    with manifest.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["filename","bytes","sha256"])
        for p in sorted(results.iterdir()):
            if p.is_file() and p.name!=manifest.name:
                w.writerow([p.name,p.stat().st_size,sha(p)])

    summary={"schema":"Stage2PrelabelStatus/v1","stage2_opened":True,"stage3_labels_opened":False,
             "colorectal_opened":False,"implementation_commit":head,
             "artifact_sha256":art.sha256(),"signature_sha256":sig["sha256"],
             "C_source":Csource,"targets":pilot,
             "prelabel_manifest_scope":"all result files except STAGE2_PRELABEL_SHA256.csv itself",
             "status":"STAGE2_PRELABEL_COMPLETE"}
    (results/"STAGE2_PRELABEL_STATUS.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    # rebuild final manifest including status
    with manifest.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["filename","bytes","sha256"])
        for p in sorted(results.iterdir()):
            if p.is_file() and p.name!=manifest.name:
                w.writerow([p.name,p.stat().st_size,sha(p)])

    print("STAGE2 PRELABEL COMPLETE")
    for p in pilot:
        print(f"{p['dataset']}: decision={p['authentic_decision']} C_target={p['C_target']:.6f} "
              f"D={p['D']:.6f} R={p['R']:.3f} false_reassurance={p['false_reassurance']:.4f} "
              f"upper95={p['false_reassurance_upper95']:.4f} calibrated={p['calibration_acceptable']}")
    print(f"Output: {results}")
    print("TARGET VALUES OPENED: YES (Stage 2, pre-hashed payload)")
    print("TARGET LABELS OPENED: NO")
    print("COLorectal OPENED: NO")

if __name__=="__main__": main()
