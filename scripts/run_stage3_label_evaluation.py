#!/usr/bin/env python3
from __future__ import annotations

import argparse,csv,hashlib,itertools,json,math,subprocess,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

EXPECTED_IMPL_TAG="preliminary-stage3-label-evaluation-implementation-freeze-2026-09-07"
EXPECTED_STAGE3_SHA="72fff8d5fe2c1785018915bf1c567d58b20184d71d5122e2fdebb001d344d30a"
TARGETS=("GSE27262","GSE32863")
SAMPLE_ALIASES=("gsm","sample_id","sample","geo_accession")
LABEL_ALIASES=("evaluation_label","label","class","phenotype","tissue","status")

def die(x): raise SystemExit("ERROR: "+str(x))

def sha(p:Path):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def git(repo,*args):
    p=subprocess.run(["git",*args],cwd=repo,capture_output=True,text=True)
    if p.returncode: die(p.stderr)
    return p.stdout.strip()

def safe_extract(zp:Path,dest:Path):
    dest.mkdir(exist_ok=False)
    with zipfile.ZipFile(zp) as z:
        for info in z.infolist():
            target=(dest/info.filename).resolve()
            if not str(target).startswith(str(dest.resolve())):
                die("unsafe Stage-3 ZIP member")
        z.extractall(dest)

def ari(a,b):
    a=np.asarray(list(a),dtype=object); b=np.asarray(list(b),dtype=object)
    if len(a)!=len(b): die("ARI length mismatch")
    n=len(a)
    if n<2: return 1.0
    ua,ia=np.unique(a,return_inverse=True)
    ub,ib=np.unique(b,return_inverse=True)
    C=np.zeros((len(ua),len(ub)),dtype=np.int64)
    np.add.at(C,(ia,ib),1)
    c2=lambda x:x*(x-1)//2
    s=int(np.sum(c2(C)))
    sr=int(np.sum(c2(C.sum(1))))
    sc=int(np.sum(c2(C.sum(0))))
    total=c2(n)
    exp=(sr*sc/total) if total else 0.0
    mx=0.5*(sr+sc)
    den=mx-exp
    return 1.0 if den==0 else float((s-exp)/den)

def entropy(labels):
    c=Counter(labels); n=sum(c.values())
    if n==0: return float("nan")
    h=0.0
    for v in c.values():
        p=v/n
        if p>0: h-=p*math.log(p)
    return h

def nmi(a,b):
    a=list(a); b=list(b)
    if len(a)!=len(b): die("NMI length mismatch")
    n=len(a)
    if n==0: return float("nan")
    ca=Counter(a); cb=Counter(b); joint=Counter(zip(a,b))
    mi=0.0
    for (x,y),v in joint.items():
        pxy=v/n; px=ca[x]/n; py=cb[y]/n
        mi+=pxy*math.log(pxy/(px*py))
    ha=entropy(a); hb=entropy(b)
    den=0.5*(ha+hb)
    return 1.0 if den==0 else float(mi/den)

def detect_column(fieldnames,aliases,kind):
    lower={c.strip().lower():c for c in fieldnames}
    hits=[lower[a] for a in aliases if a in lower]
    if len(hits)!=1:
        die(f"expected exactly one {kind} column from {aliases}; got {hits} in {fieldnames}")
    return hits[0]

def locate_label_file(root,dataset):
    exact=f"{dataset}_evaluation_labels.csv"
    hits=[p for p in root.rglob(exact) if p.is_file()]
    if len(hits)!=1: die(f"expected exactly one {exact}; found {len(hits)}")
    return hits[0]

def read_labels(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        if not rd.fieldnames: die(f"empty label header: {path}")
        gsmcol=detect_column(rd.fieldnames,SAMPLE_ALIASES,"sample identifier")
        labcol=detect_column(rd.fieldnames,LABEL_ALIASES,"evaluation label")
        out={}
        for r in rd:
            gsm=(r.get(gsmcol) or "").strip()
            lab=(r.get(labcol) or "").strip()
            if not gsm or not lab: die(f"missing GSM/label in {path}")
            if gsm in out: die(f"duplicate GSM {gsm} in {path}")
            out[gsm]=lab
    if not out: die(f"no labels in {path}")
    return out,gsmcol,labcol

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def forced_profile(row,profile_ids):
    vals=json.loads(row["profile_scores"])
    if len(vals)!=len(profile_ids): die("profile_scores length mismatch")
    cand=[]
    for pid,v in zip(profile_ids,vals):
        if v is not None and math.isfinite(float(v)):
            cand.append((float(v),pid))
    if not cand: return None
    # score descending, profile lexical ascending on tie
    cand.sort(key=lambda z:(-z[0],z[1]))
    return cand[0][1]

def optimal_binary_mapping(pred,true):
    profiles=sorted(set(pred)); labels=sorted(set(true))
    if len(profiles)!=2 or len(labels)!=2:
        return None,float("nan"),float("nan")
    options=[]
    for perm in itertools.permutations(labels):
        mp=dict(zip(profiles,perm))
        mapped=[mp[x] for x in pred]
        correct=sum(x==y for x,y in zip(mapped,true))
        key=json.dumps(mp,sort_keys=True,separators=(",",":"))
        options.append((-correct,key,mp,mapped))
    options.sort(key=lambda x:(x[0],x[1]))
    _,_,mp,mapped=options[0]
    acc=float(np.mean([x==y for x,y in zip(mapped,true)]))
    recalls=[]
    for lab in labels:
        inds=[i for i,y in enumerate(true) if y==lab]
        recalls.append(float(np.mean([mapped[i]==lab for i in inds])) if inds else float("nan"))
    bal=float(np.mean(recalls))
    return mp,acc,bal

def write_csv(path,rows,fieldnames=None):
    if fieldnames is None:
        fieldnames=[]
        for r in rows:
            for k in r:
                if k not in fieldnames: fieldnames.append(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

def stage2_summary(repo):
    p=repo/"docs/preliminary/stage2_prelabel_results"
    authentic=read_csv(p/"authentic_transfer.csv")
    false=read_csv(p/"false_reassurance_summary.csv")
    B=read_csv(p/"perturbation_B.csv")
    C=read_csv(p/"perturbation_C.csv")
    D=read_csv(p/"perturbation_D.csv")
    return authentic,false,B,C,D

def build_final_markdown(repo,label_rows):
    source=json.loads((repo/"docs/preliminary/v3_source_results/SOURCE_RESULT_STATUS.json").read_text(encoding="utf-8-sig"))
    authentic,false,B,C,D=stage2_summary(repo)
    lm={r["dataset"]:r for r in label_rows}
    lines=[]
    lines.append("# Final development-exposed lung preliminary evidence")
    lines.append("")
    lines.append("## Audit chronology")
    lines.append("")
    lines.append("- V1 remained closed as `NO_STABLE_STRUCTURE`.")
    lines.append("- V2 replaced the absolute stability cliff with a prospectively frozen pre-target null-calibrated stability gate; it passed that gate but stopped for insufficient signature-subset index support.")
    lines.append("- V3 froze an assignment-stratified 12/18 holdout allocation before validation-gene analysis and then obtained `SOURCE_REFERENCE_PASS`.")
    lines.append("- Stage 2 executed both authentic targets and froze all structural results before Stage-3 labels were opened.")
    lines.append("- The evidence is therefore development-exposed preliminary evidence, not untouched confirmation.")
    lines.append("")
    lines.append("## Source reference")
    lines.append("")
    lines.append(f"- Frozen artifact SHA-256: `{source['artifact_sha256']}`")
    lines.append(f"- Frozen validation signature: {source['selected_signature_relations']} relations; SHA-256 `{source['signature_sha256']}`")
    lines.append(f"- `C_source = {float(source['C_source']):.3f}`; one-sided 95% lower bound `{float(source['source_reference_lower95']):.3f}`.")
    lines.append("")
    lines.append("## Authentic cross-cohort structural retention")
    lines.append("")
    for r in authentic:
        ds=r["dataset"]
        fr=next(x for x in false if x["dataset"]==ds)
        lines.append(
            f"- **{ds}: {r['decision']}** — C_target={float(r['C_target']):.3f}, "
            f"D={float(r['D']):.3f}, R={float(r['R']):.3f}; "
            f"assignment coverage={float(r['specimen_assignment_coverage']):.3f} specimen / "
            f"{float(r['participant_assignment_coverage']):.3f} participant; "
            f"conditional false reassurance={int(fr['passes'])}/{int(fr['replicates'])}, "
            f"exact 95% upper bound={float(fr['exact95_upper']):.4f}."
        )
    lines.append("")
    lines.append("Both authentic targets therefore met the frozen structural PASS rule, and neither target produced a PASS in 999 conditional validation-block randomisations.")
    lines.append("")
    lines.append("## Mechanistic controls")
    lines.append("")
    for ds in TARGETS:
        br=[r for r in B if r["dataset"]==ds]
        all_zero=all(int(r["strict_relation_flips"])==0 for r in br)
        min_geom=min(float(r["value_cluster_intact_vs_perturbed_ari"]) for r in br if float(r["t"])>0)
        lines.append(f"- **B / {ds}:** frozen assignment/signature relation flips were {'zero at every dose' if all_zero else 'non-zero'}; target value-space clustering agreement with the intact geometry fell as low as ARI={min_geom:.3f}.")
        cr=[r for r in C if r["dataset"]==ds]
        dose_vals=sorted({float(r["dose"]) for r in cr})
        parts=[]
        for d in dose_vals:
            rr=[x for x in cr if float(x["dose"])==d]
            parts.append(f"{d:g}: {sum(x['decision']=='PASS' for x in rr)}/{len(rr)} PASS")
        lines.append(f"- **C / {ds}:** " + "; ".join(parts) + ".")
        dr=sorted([r for r in D if r["dataset"]==ds],key=lambda x:float(x["dose"]))
        nonpass=next((r for r in dr if r["decision"]!="PASS"),None)
        if nonpass:
            lines.append(f"- **D / {ds}:** first non-PASS state occurred at core-gene deletion dose {float(nonpass['dose']):g}, with assigned fraction {float(nonpass['assigned_fraction']):.3f}; complete core loss yielded no assigned samples.")
    lines.append("")
    lines.append("## Post-freeze external-label agreement")
    lines.append("")
    for ds in TARGETS:
        r=lm[ds]
        lines.append(
            f"- **{ds}:** accepted-only ARI={float(r['accepted_ari']):.3f}, "
            f"NMI={float(r['accepted_nmi']):.3f} on {int(r['accepted_n'])}/{int(r['total_n'])} specimens; "
            f"forced ARI={float(r['forced_ari']):.3f}, NMI={float(r['forced_nmi']):.3f} "
            f"on {int(r['forced_n'])}/{int(r['total_n'])} specimens."
        )
        if r.get("accepted_mapped_accuracy","") not in ("","nan"):
            lines.append(
                f"  Post-hoc descriptive binary mapping: accepted accuracy={float(r['accepted_mapped_accuracy']):.3f}, "
                f"balanced accuracy={float(r['accepted_mapped_balanced_accuracy']):.3f}; "
                f"forced accuracy={float(r['forced_mapped_accuracy']):.3f}, "
                f"balanced accuracy={float(r['forced_mapped_balanced_accuracy']):.3f}."
            )
    lines.append("")
    lines.append("## Interpretation ceiling")
    lines.append("")
    lines.append("These results support feasibility of frozen source-defined relational execution, independent source-side structural contrast, and retention of that contrast in two development-exposed lung targets under explicit abstention and conditional false-reassurance controls. They do **not** establish a universal biological subtype, clinical utility, H2 family-level predictive validity, or independent colorectal confirmation.")
    lines.append("")
    return "\n".join(lines)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--prep-root",required=True);a=ap.parse_args()
    root=Path(a.prep_root).resolve(); repo=root/"repo_patch/_work_repo"

    head=git(repo,"rev-parse","HEAD")
    tag=git(repo,"rev-list","-n","1",EXPECTED_IMPL_TAG)
    if head!=tag: die("HEAD is not exact frozen Stage-3 implementation commit")
    if git(repo,"status","--porcelain"): die("repository dirty before label access")

    # Verify frozen Stage-2 status before opening labels.
    s2=json.loads((repo/"docs/preliminary/stage2_prelabel_results/STAGE2_PRELABEL_STATUS.json").read_text(encoding="utf-8-sig"))
    if s2.get("status")!="STAGE2_PRELABEL_COMPLETE": die("Stage-2 prelabel package incomplete")
    if s2.get("stage3_labels_opened") is not False: die("Stage-2 status does not say labels sealed")

    name="POST_ASSIGNMENT_LABELS_2026-09-07.zip"
    zp=next((p for p in (root/name,root.parent/name) if p.exists()),None)
    if zp is None: die(f"{name} not found in prep root or parent")
    if sha(zp)!=EXPECTED_STAGE3_SHA: die("Stage-3 ZIP SHA-256 mismatch")

    opened=root/"STAGE3_OPENED_LABELS"
    if opened.exists(): die("Stage-3 extraction directory already exists; refusing duplicate label access")
    safe_extract(zp,opened)

    out=root/"STAGE3_LABEL_RESULTS"; out.mkdir(exist_ok=False)
    access={"schema":"Stage3LabelAccessRecord/v1",
            "opened_at_utc":datetime.now(timezone.utc).isoformat(),
            "stage3_zip":str(zp),"stage3_sha256":sha(zp),
            "implementation_commit":head,
            "prelabel_commit":"f633d7cc4ae3c26e794798f3e4ec3194ad7e8aa0",
            "labels_may_change_prelabel_decisions":False}
    (out/"STAGE3_LABEL_ACCESS_RECORD.json").write_text(json.dumps(access,indent=2)+"\n",encoding="utf-8")

    artifact=json.loads((repo/"docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json").read_text(encoding="utf-8-sig"))
    profile_ids=[p["profile_id"] for p in artifact["prototypes"]]

    agreement=[]; contingency=[]; joined=[]
    for ds in TARGETS:
        labels,gsmcol,labcol=read_labels(locate_label_file(opened,ds))
        ass=read_csv(repo/f"docs/preliminary/stage2_prelabel_results/{ds}_prelabel_assignments.csv")
        if len({r["gsm"] for r in ass})!=len(ass): die(f"duplicate assignment GSM in {ds}")
        agsms={r["gsm"] for r in ass}
        if agsms!=set(labels): 
            missing=sorted(agsms-set(labels)); extra=sorted(set(labels)-agsms)
            die(f"{ds} label GSM set mismatch; missing={missing[:5]} extra={extra[:5]}")
        if len(set(labels.values()))<2: die(f"{ds} has fewer than two evaluation labels")

        accepted_pred=[]; accepted_true=[]
        forced_pred=[]; forced_true=[]
        for r in ass:
            lab=labels[r["gsm"]]
            fp=forced_profile(r,profile_ids)
            row=dict(r); row["evaluation_label"]=lab; row["forced_profile_id"]=fp
            joined.append(row)
            if r["assignment"]=="ASSIGNED":
                accepted_pred.append(r["best_profile_id"]); accepted_true.append(lab)
            if fp is not None:
                forced_pred.append(fp); forced_true.append(lab)

        if not accepted_pred: die(f"{ds}: no accepted assignments")
        if not forced_pred: die(f"{ds}: no forced assignments")

        amap,aacc,abal=optimal_binary_mapping(accepted_pred,accepted_true)
        fmap,facc,fbal=optimal_binary_mapping(forced_pred,forced_true)
        row={
            "dataset":ds,"total_n":len(ass),
            "label_classes":json.dumps(sorted(set(labels.values()))),
            "accepted_n":len(accepted_pred),"accepted_fraction":len(accepted_pred)/len(ass),
            "accepted_ari":ari(accepted_pred,accepted_true),
            "accepted_nmi":nmi(accepted_pred,accepted_true),
            "accepted_mapping":json.dumps(amap,sort_keys=True) if amap else "",
            "accepted_mapped_accuracy":aacc,
            "accepted_mapped_balanced_accuracy":abal,
            "forced_n":len(forced_pred),"forced_fraction":len(forced_pred)/len(ass),
            "forced_ari":ari(forced_pred,forced_true),
            "forced_nmi":nmi(forced_pred,forced_true),
            "forced_mapping":json.dumps(fmap,sort_keys=True) if fmap else "",
            "forced_mapped_accuracy":facc,
            "forced_mapped_balanced_accuracy":fbal,
            "label_file_sample_column":gsmcol,"label_file_label_column":labcol,
        }
        agreement.append(row)

        for pop,pred,true in (("accepted",accepted_pred,accepted_true),("forced",forced_pred,forced_true)):
            cnt=Counter(zip(pred,true))
            for prof in sorted(set(pred)):
                for lab in sorted(set(true)):
                    contingency.append({"dataset":ds,"population":pop,
                        "profile_id":prof,"evaluation_label":lab,"count":cnt.get((prof,lab),0)})

    write_csv(out/"label_agreement.csv",agreement)
    write_csv(out/"label_contingency.csv",contingency)
    write_csv(out/"target_assignments_with_labels.csv",joined)

    # Combined final table.
    authentic=read_csv(repo/"docs/preliminary/stage2_prelabel_results/authentic_transfer.csv")
    false=read_csv(repo/"docs/preliminary/stage2_prelabel_results/false_reassurance_summary.csv")
    lm={r["dataset"]:r for r in agreement}
    final=[]
    for r in authentic:
        ds=r["dataset"]; fr=next(x for x in false if x["dataset"]==ds); la=lm[ds]
        final.append({
            "dataset":ds,"structural_decision":r["decision"],
            "C_target":r["C_target"],"D":r["D"],"R":r["R"],
            "specimen_assignment_coverage":r["specimen_assignment_coverage"],
            "participant_assignment_coverage":r["participant_assignment_coverage"],
            "false_reassurance":fr["false_reassurance"],
            "false_reassurance_exact95_upper":fr["exact95_upper"],
            "accepted_label_ari":la["accepted_ari"],"accepted_label_nmi":la["accepted_nmi"],
            "accepted_label_n":la["accepted_n"],
            "forced_label_ari":la["forced_ari"],"forced_label_nmi":la["forced_nmi"],
            "forced_label_n":la["forced_n"],
            "accepted_posthoc_mapped_accuracy":la["accepted_mapped_accuracy"],
            "forced_posthoc_mapped_accuracy":la["forced_mapped_accuracy"]
        })
    write_csv(out/"FINAL_PRELIMINARY_EVIDENCE_TABLE.csv",final)
    (out/"FINAL_PRELIMINARY_EVIDENCE_SUMMARY.md").write_text(build_final_markdown(repo,agreement),encoding="utf-8")

    status={"schema":"Stage3LabelStatus/v1","status":"STAGE3_LABEL_EVALUATION_COMPLETE",
            "labels_opened":True,"implementation_commit":head,
            "prelabel_commit":"f633d7cc4ae3c26e794798f3e4ec3194ad7e8aa0",
            "prelabel_decisions_modified":False,
            "targets":agreement}
    (out/"STAGE3_LABEL_STATUS.json").write_text(json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

    manifest=out/"STAGE3_RESULTS_SHA256.csv"
    with manifest.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["filename","bytes","sha256"])
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name!=manifest.name:
                w.writerow([p.name,p.stat().st_size,sha(p)])

    print("STAGE3 LABEL EVALUATION COMPLETE")
    for r in agreement:
        print(
            f"{r['dataset']}: accepted ARI={r['accepted_ari']:.6f} NMI={r['accepted_nmi']:.6f} "
            f"n={r['accepted_n']}/{r['total_n']} | forced ARI={r['forced_ari']:.6f} "
            f"NMI={r['forced_nmi']:.6f} n={r['forced_n']}/{r['total_n']}"
        )
        if math.isfinite(float(r["accepted_mapped_accuracy"])):
            print(
                f"  post-hoc mapped accuracy accepted={r['accepted_mapped_accuracy']:.6f} "
                f"forced={r['forced_mapped_accuracy']:.6f}"
            )
    print(f"Output: {out}")
    print("PRELABEL STRUCTURAL DECISIONS MODIFIED: NO")
    print("TARGET LABELS OPENED: YES")
    print("COLORECTAL OPENED: NO")

if __name__=="__main__": main()
