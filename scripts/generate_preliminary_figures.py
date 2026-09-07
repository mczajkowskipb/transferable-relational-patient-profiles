#!/usr/bin/env python3
"""Regenerate compact SVG figures from committed preliminary-result CSV files."""

from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
PRE = REPO / "docs" / "preliminary"
OUT = PRE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def save():
    plt.tight_layout()
    plt.savefig(CURRENT, format="svg", bbox_inches="tight")
    plt.close()

# 1. Structural retention.
auth = rows(PRE / "stage2_prelabel_results" / "authentic_transfer.csv")
names = [r["dataset"] for r in auth]
vals = [float(r["C_target"]) for r in auth]
CURRENT = OUT / "structural_retention_overview.svg"
plt.figure(figsize=(6.2, 4.2))
plt.bar(names, vals)
plt.axhline(0.5 * 0.6053475935828878, linestyle="--", label="0.5 × C_source")
plt.axhline(0.10, linestyle=":", label="absolute C_target floor")
plt.ylabel("Independent validation contrast")
plt.title("Frozen structural retention in external lung cohorts")
plt.ylim(0, 1.0)
plt.legend()
for i, r in enumerate(auth):
    plt.text(i, float(r["C_target"]) + 0.025, f"PASS\nR={float(r['R']):.2f}", ha="center")
save()

# 2. Validation damage.
c = rows(PRE / "stage2_prelabel_results" / "perturbation_C.csv")
CURRENT = OUT / "signature_damage_pass_fraction.svg"
plt.figure(figsize=(6.2, 4.2))
for ds in sorted(set(r["dataset"] for r in c)):
    rr = [r for r in c if r["dataset"] == ds]
    doses = sorted(set(float(r["dose"]) for r in rr))
    frac = []
    for d in doses:
        z = [r for r in rr if float(r["dose"]) == d]
        frac.append(sum(r["decision"] == "PASS" for r in z) / len(z))
    plt.plot(doses, frac, marker="o", label=ds)
plt.xlabel("Fraction of frozen validation-signature genes disrupted")
plt.ylabel("PASS fraction")
plt.ylim(-0.03, 1.03)
plt.title("Structural PASS collapses under validation-specific damage")
plt.legend()
save()

# 3. Core dropout.
d = rows(PRE / "stage2_prelabel_results" / "perturbation_D.csv")
CURRENT = OUT / "core_dropout_assigned_fraction.svg"
plt.figure(figsize=(6.2, 4.2))
for ds in sorted(set(r["dataset"] for r in d)):
    rr = sorted([r for r in d if r["dataset"] == ds], key=lambda r: float(r["dose"]))
    plt.plot(
        [float(r["dose"]) for r in rr],
        [float(r["assigned_fraction"]) for r in rr],
        marker="o",
        label=ds,
    )
plt.xlabel("Fraction of frozen assignment-core genes removed")
plt.ylabel("Assigned specimen fraction")
plt.ylim(-0.03, 1.03)
plt.title("Core feature loss produces abstention")
plt.legend()
save()

# 4. Geometry perturbation.
b = rows(PRE / "stage2_prelabel_results" / "perturbation_B.csv")
CURRENT = OUT / "geometry_perturbation_ari.svg"
plt.figure(figsize=(6.2, 4.2))
for ds in sorted(set(r["dataset"] for r in b)):
    rr = sorted([r for r in b if r["dataset"] == ds], key=lambda r: float(r["t"]))
    plt.plot(
        [float(r["t"]) for r in rr],
        [float(r["value_cluster_intact_vs_perturbed_ari"]) for r in rr],
        marker="o",
        label=ds,
    )
plt.xlabel("Order-preserving geometry perturbation t")
plt.ylabel("ARI: intact vs perturbed value-space clustering")
plt.title("Absolute-value geometry changes while frozen relations remain unchanged")
plt.legend()
save()

print("Wrote figures to", OUT)
