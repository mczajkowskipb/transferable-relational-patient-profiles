> **Legacy document notice — 2026-09-07.** Retained for provenance. This file predates the current structural-retention grant design and must not be treated as its current protocol. See [pre-pilot audit](../preliminary/AUDIT_2026_09_07.md). The independent-signature pilot is not implemented or frozen yet.

# 21. NOVELTY MATRIX - REVIEWER DEFENCE

The table is intentionally phrased as a **combination-of-properties gap**, not a claim that every individual component is unprecedented.

| Approach | Unsupervised group discovery | Group definition itself is executable within-sample rules | Frozen single-sample execution without target reclustering | Explicit abstention / rejection | Source-derived transportability certificate + shift boundary |
|---|---:|---:|---:|---:|---:|
| TSP / supervised Relative Expression Analysis | No | Predictor rules, not discovered groups | Yes | Sometimes | No |
| REO / rank transform -> conventional clustering | Yes | No; clustering occurs after transformation | Usually no direct frozen group definition | Rarely | No |
| Order-preserving biclustering | Yes | Order constraints describe biclusters | Not generally a frozen patient-profile assignment object | No | No |
| Interpretable / rule-based clustering | Yes | Sometimes | Not generally the central source-target object | Rarely | No |
| Cluster -> supervised classifier | Discovery yes; transfer step supervised | Classifier is post-hoc, not the discovered group definition | Yes | Possible | Usually no |
| Certified robustness / shift certificates for predictive models | Usually supervised | No: certificate concerns a predictor, not an unsupervised group definition | Yes for predictor | Depends on method | Yes, but for a different predictive object/shift model |
| **Proposed RPP framework** | **Yes** | **Yes: sparse within-sample partial-order profile learned jointly with membership** | **Yes** | **Yes: NO_STABLE_STRUCTURE / UNASSIGNED** | **Yes, conditional on declared shift classes and prospectively tested** |

## One-sentence novelty statement

**The novelty is the joint formulation of an unsupervised patient group as a sparse executable within-sample relational profile that is itself the frozen single-sample transfer object, augmented with explicit abstention and an independently source-calibrated transportability radius whose validity is tested under controlled and real cohort/platform shift.**

## Reviewer boundary

Do not claim:
- first rule-based clustering;
- first use of relative expression orderings;
- first interpretable clustering;
- first single-sample predictor;
- universal domain invariance;
- first certified-robustness method under distribution shift.

Claim the exact combination and test it against matched post-hoc and transform-then-cluster baselines.
