#!/usr/bin/env python3
"""Public replay of the final frozen lung RPP artifact/signature from NCBI GEO files.

This script is intentionally separate from the historical stage-aware runners.
It does not recreate the original blinding chronology. Instead, it downloads/uses
the public GEO Series Matrix representation and verifies that the final frozen
scientific objects reproduce the committed source/target quantities.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import re
from pathlib import Path

import numpy as np

REQUIRED = {
    "GSE19804": "GSE19804_series_matrix.txt.gz",
    "GSE27262": "GSE27262_series_matrix.txt.gz",
    "GSE32863": "GSE32863_series_matrix.txt.gz",
    "GPL570": "GPL570.annot.gz",
    "GPL6884": "GPL6884.annot.gz",
}
TARGET_PLATFORM = {"GSE27262": "GPL570", "GSE32863": "GPL6884"}
EXPECTED_COMMON_GENES = 17062
TOL = 1e-10

def die(msg: str) -> None:
    raise SystemExit("ERROR: " + msg)

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        die(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def parse_annotation(path: Path) -> dict[str, str]:
    p2g: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        header = None
        for line in f:
            if line.startswith("ID\t"):
                header = next(csv.reader([line], delimiter="\t"))
                break
        if header is None or "Gene ID" not in header:
            die(f"cannot find ID/Gene ID header in {path}")
        rd = csv.DictReader(f, fieldnames=header, delimiter="\t")
        for row in rd:
            probe = (row.get("ID") or "").strip()
            gid = (row.get("Gene ID") or "").strip()
            if probe and re.fullmatch(r"\d+", gid):
                p2g[probe] = gid
    return p2g

def _meta_values(line: str) -> list[str]:
    row = next(csv.reader([line], delimiter="\t"))
    return [x.strip() for x in row[1:]]

def read_series_matrix(
    path: Path,
    probe_to_gene: dict[str, str],
    common_genes: list[str],
) -> tuple[np.ndarray, list[str], dict[str, str]]:
    """Map one GEO Series Matrix directly to the frozen Entrez namespace."""
    common = set(common_genes)
    by_gene: dict[str, list[np.ndarray]] = {}
    geo_accessions: list[str] | None = None
    titles: list[str] | None = None
    table_gsms: list[str] | None = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        in_table = False
        for line in f:
            if not in_table:
                if line.startswith("!Sample_geo_accession"):
                    geo_accessions = _meta_values(line)
                elif line.startswith("!Sample_title"):
                    titles = _meta_values(line)
                elif line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    header = next(csv.reader([next(f)], delimiter="\t"))
                    table_gsms = [x.strip().strip('"') for x in header[1:]]
                continue

            if line.startswith("!series_matrix_table_end"):
                break
            row = next(csv.reader([line], delimiter="\t"))
            if not row:
                continue
            probe = row[0].strip().strip('"')
            gid = probe_to_gene.get(probe)
            if gid is None or gid not in common:
                continue
            vals = np.array(
                [
                    float(x) if x not in ("", "NA", "NaN", "nan", "null") else np.nan
                    for x in row[1:]
                ],
                dtype=float,
            )
            by_gene.setdefault(gid, []).append(vals)

    if table_gsms is None:
        die(f"Series Matrix table not found in {path}")
    if geo_accessions is None or titles is None:
        die(f"Sample_geo_accession / Sample_title metadata not found in {path}")
    geo_accessions = [x.strip().strip('"') for x in geo_accessions]
    titles = [x.strip().strip('"') for x in titles]
    if len(geo_accessions) != len(titles):
        die(f"metadata length mismatch in {path}")
    title_by_gsm = dict(zip(geo_accessions, titles))
    if set(table_gsms) != set(title_by_gsm):
        die(f"Series Matrix table GSMs do not match metadata GSMs in {path}")

    gene_index = {g: i for i, g in enumerate(common_genes)}
    X = np.full((len(table_gsms), len(common_genes)), np.nan, dtype=float)
    for gid, arrays in by_gene.items():
        X[:, gene_index[gid]] = np.nanmedian(np.vstack(arrays), axis=0)

    # Historical target mapping sorted GSM lexically.
    order = np.argsort(np.asarray(table_gsms, dtype=str))
    gsms = [table_gsms[i] for i in order]
    X = X[order]
    titles_sorted = {g: title_by_gsm[g] for g in gsms}
    return X, gsms, titles_sorted

def participant_and_label(dataset: str, title: str) -> tuple[str, str]:
    if dataset == "GSE32863":
        m = re.search(r"([A-Za-z0-9]+)_([NT])(?:\s|\(|$)", title)
    else:
        m = re.search(r"([0-9]+)([NT])(?:\s|\(|$)", title)
    if not m:
        die(f"cannot derive participant/tissue code from {dataset} title: {title!r}")
    pid, code = m.group(1), m.group(2)
    return pid, "tumor" if code == "T" else "adjacent_normal"

def build_metadata(dataset: str, gsms: list[str], titles: dict[str, str]) -> tuple[list[dict], dict[str, str]]:
    tmp = []
    counts: dict[str, int] = {}
    labels: dict[str, str] = {}
    for gsm in gsms:
        pid, lab = participant_and_label(dataset, titles[gsm])
        tmp.append((gsm, pid))
        labels[gsm] = lab
        counts[pid] = counts.get(pid, 0) + 1
    rows = [
        {
            "dataset": dataset,
            "gsm": gsm,
            "participant_id": pid,
            "pair_status": "paired" if counts[pid] == 2 else "singleton",
        }
        for gsm, pid in tmp
    ]
    return rows, labels

def csv_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def almost(a: float, b: float, tol: float = TOL) -> bool:
    return abs(float(a) - float(b)) <= tol

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/public_geo")
    ap.add_argument("--output-dir", default="public_replay_outputs")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    data = Path(args.data_dir).resolve()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    for key, filename in REQUIRED.items():
        if not (data / filename).is_file():
            die(f"missing {filename}; run scripts/download_public_geo_inputs.py first")

    p570 = parse_annotation(data / REQUIRED["GPL570"])
    p6884 = parse_annotation(data / REQUIRED["GPL6884"])
    g570 = set(p570.values())
    g6884 = set(p6884.values())
    common_genes = sorted(g570 & g6884, key=int)
    if len(common_genes) != EXPECTED_COMMON_GENES:
        die(f"common gene universe is {len(common_genes)}, expected {EXPECTED_COMMON_GENES}")

    s2 = load_module(repo / "scripts" / "run_stage2_prelabel.py", "rpp_stage2_public_replay")
    s3 = load_module(repo / "scripts" / "run_stage3_label_evaluation.py", "rpp_stage3_public_replay")
    v1 = load_module(repo / "scripts" / "run_source_fit_v1.py", "rpp_source_v1_public_replay")

    cfg = json.loads((repo / "docs/preliminary/V1_CONFIG.json").read_text(encoding="utf-8-sig"))
    art = s2.load_artifact(repo / "docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json")
    sig = json.loads((repo / "docs/preliminary/v3_source_results/VALIDATION_SIGNATURE_V3.json").read_text(encoding="utf-8-sig"))
    source_status = json.loads((repo / "docs/preliminary/v3_source_results/SOURCE_RESULT_STATUS.json").read_text(encoding="utf-8-sig"))
    C_source_expected = float(source_status["C_source"])

    # Source reconstruction.
    Xs, source_gsms, source_titles = read_series_matrix(
        data / REQUIRED["GSE19804"], p570, common_genes
    )
    source_meta, _ = build_metadata("GSE19804", source_gsms, source_titles)
    source_gsm_to_row = {g: i for i, g in enumerate(source_gsms)}

    split = csv_rows(repo / "docs/preliminary/v3_source_split/SOURCE_PARTICIPANT_SPLIT_V3.csv")
    ref_pids = [r["participant_id"] for r in split if r["role"] == "SOURCE_REFERENCE_V3"]
    by_pid: dict[str, list[int]] = {}
    for r in source_meta:
        by_pid.setdefault(r["participant_id"], []).append(source_gsm_to_row[r["gsm"]])
    missing_ref = [p for p in ref_pids if p not in by_pid]
    if missing_ref:
        die(f"source reference participants missing from public metadata: {missing_ref}")

    ref_rows = np.array([i for p in ref_pids for i in sorted(by_pid[p])], dtype=int)
    ref_gsms = [source_gsms[i] for i in ref_rows]
    meta_by_gsm = {r["gsm"]: r for r in source_meta}
    ref_meta = [meta_by_gsm[g] for g in ref_gsms]
    ref_exec = s2.execute_rpp_artifact(Xs[ref_rows], common_genes, art)
    ref_vs = s2.strict_signature_scores(Xs[ref_rows], common_genes, sig)
    source_pids, sai, sac, _ = s2.group_arrays(
        ref_exec, ref_vs, ref_meta, ref_gsms, sig["index_profile"]
    )

    vals = {"INDEX": {}, "COMPARATOR": {}}
    for p, a, c in zip(source_pids, sai, sac):
        if np.isfinite(a):
            vals["INDEX"][p] = float(a)
        if np.isfinite(c):
            vals["COMPARATOR"][p] = float(c)
    C_source = v1.contrast_from_group_values(vals)
    source_boot = v1.bootstrap_contrast(
        vals,
        source_pids,
        cfg["source_reference"]["bootstrap_replicates"],
        20260914,
    )
    source_lower = float(
        np.quantile(
            source_boot,
            cfg["source_reference"]["source_one_sided_alpha"],
            method="lower",
        )
    )

    if not almost(C_source, C_source_expected):
        die(f"C_source replay mismatch: {C_source} vs {C_source_expected}")
    if not almost(source_lower, float(source_status["source_reference_lower95"])):
        die(f"source lower95 replay mismatch: {source_lower}")

    expected_auth = {
        r["dataset"]: r
        for r in csv_rows(repo / "docs/preliminary/stage2_prelabel_results/authentic_transfer.csv")
    }
    expected_lab = {
        r["dataset"]: r
        for r in csv_rows(repo / "docs/preliminary/stage3_label_results/label_agreement.csv")
    }

    output_rows = []
    checks = {
        "common_genes": len(common_genes),
        "C_source": C_source,
        "source_lower95": source_lower,
        "targets": {},
        "all_checks_passed": True,
    }

    for dataset, platform in TARGET_PLATFORM.items():
        p2g = p570 if platform == "GPL570" else p6884
        X, gsms, titles = read_series_matrix(data / REQUIRED[dataset], p2g, common_genes)
        meta, labels = build_metadata(dataset, gsms, titles)

        ex = s2.execute_rpp_artifact(X, common_genes, art)
        vs = s2.strict_signature_scores(X, common_genes, sig)
        pids, ai, ac, support = s2.group_arrays(ex, vs, meta, gsms, sig["index_profile"])
        tinds = s2.bootstrap_indices(len(pids), int(dataset.replace("GSE", "")), 2000)
        dec = s2.structural_decision(ai, ac, support, source_boot, tinds, cfg, C_source)

        accepted_pred = []
        accepted_true = []
        forced_pred = []
        forced_true = []
        profile_ids = [p.profile_id for p in art.prototypes]

        for gsm, e in zip(gsms, ex):
            lab = labels[gsm]
            if e.assignment == "ASSIGNED":
                accepted_pred.append(e.best_profile_id)
                accepted_true.append(lab)
            row = {"profile_scores": json.dumps(e.profile_scores)}
            forced = s3.forced_profile(row, profile_ids)
            if forced is not None:
                forced_pred.append(forced)
                forced_true.append(lab)

        accepted_ari = s3.ari(accepted_pred, accepted_true)
        accepted_nmi = s3.nmi(accepted_pred, accepted_true)
        forced_ari = s3.ari(forced_pred, forced_true)
        forced_nmi = s3.nmi(forced_pred, forced_true)

        ea = expected_auth[dataset]
        el = expected_lab[dataset]
        target_checks = {
            "decision": dec["decision"] == ea["decision"],
            "C_target": almost(dec["C_target"], float(ea["C_target"])),
            "R": almost(dec["R"], float(ea["R"])),
            "accepted_ari": almost(accepted_ari, float(el["accepted_ari"])),
            "accepted_nmi": almost(accepted_nmi, float(el["accepted_nmi"])),
            "forced_ari": almost(forced_ari, float(el["forced_ari"])),
            "forced_nmi": almost(forced_nmi, float(el["forced_nmi"])),
        }
        checks["targets"][dataset] = target_checks
        checks["all_checks_passed"] = checks["all_checks_passed"] and all(target_checks.values())

        output_rows.append({
            "dataset": dataset,
            "structural_decision": dec["decision"],
            "C_target": dec["C_target"],
            "D": dec["D"],
            "R": dec["R"],
            "specimen_assignment_coverage": support["specimen_assignment_coverage"],
            "participant_assignment_coverage": support["participant_assignment_coverage"],
            "accepted_n": len(accepted_pred),
            "accepted_ari": accepted_ari,
            "accepted_nmi": accepted_nmi,
            "forced_n": len(forced_pred),
            "forced_ari": forced_ari,
            "forced_nmi": forced_nmi,
        })

    out_csv = out / "public_replay_summary.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(output_rows[0]))
        w.writeheader()
        w.writerows(output_rows)

    (out / "public_replay_checks.json").write_text(
        json.dumps(checks, indent=2) + "\n", encoding="utf-8"
    )

    print(f"Common Entrez universe: {len(common_genes)}")
    print(f"C_source={C_source:.10f}; lower95={source_lower:.10f}")
    for r in output_rows:
        print(
            f"{r['dataset']}: {r['structural_decision']} "
            f"C_target={r['C_target']:.10f} R={r['R']:.10f} "
            f"accepted ARI={r['accepted_ari']:.10f} NMI={r['accepted_nmi']:.10f}"
        )
    print("ALL CHECKS PASSED:", checks["all_checks_passed"])
    print("Output:", out)

    if not checks["all_checks_passed"]:
        raise SystemExit(2)

if __name__ == "__main__":
    main()
