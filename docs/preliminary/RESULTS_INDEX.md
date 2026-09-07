# Preliminary results index

This page is the shortest route through the complete development-exposed lung pilot.

## Final evidence

- [`stage3_label_results/FINAL_PRELIMINARY_EVIDENCE_SUMMARY.md`](stage3_label_results/FINAL_PRELIMINARY_EVIDENCE_SUMMARY.md) — human-readable final summary.
- [`stage3_label_results/FINAL_PRELIMINARY_EVIDENCE_TABLE.csv`](stage3_label_results/FINAL_PRELIMINARY_EVIDENCE_TABLE.csv) — compact final machine-readable table.
- [`stage3_label_results/label_agreement.csv`](stage3_label_results/label_agreement.csv) — accepted/forced ARI, NMI and post-hoc mapped accuracy.
- [`stage2_prelabel_results/authentic_transfer.csv`](stage2_prelabel_results/authentic_transfer.csv) — structural endpoint frozen before labels.
- [`stage2_prelabel_results/false_reassurance_summary.csv`](stage2_prelabel_results/false_reassurance_summary.csv) — 999 conditional nulls per target.
- [`v3_source_results/SOURCE_RESULT_STATUS.json`](v3_source_results/SOURCE_RESULT_STATUS.json) — final independent source-reference status.

## Frozen scientific objects

- [`v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json`](v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json) — frozen executable assignment artifact.
- [`v3_source_results/VALIDATION_SIGNATURE_V3.json`](v3_source_results/VALIDATION_SIGNATURE_V3.json) — frozen feature-disjoint validation signature.
- [`v3_source_split/SOURCE_PARTICIPANT_SPLIT_V3.csv`](v3_source_split/SOURCE_PARTICIPANT_SPLIT_V3.csv) — exact assignment-stratified V3 holdout allocation.
- [`SOURCE_PARTICIPANT_SPLIT_V1.csv`](SOURCE_PARTICIPANT_SPLIT_V1.csv) — original source split.

## Source-development chronology

### V1

- [`v1_source_results/V1_SOURCE_OUTCOME.md`](v1_source_results/V1_SOURCE_OUTCOME.md)
- [`v1_source_results/source_fit_summary.csv`](v1_source_results/source_fit_summary.csv)

Formal result: `NO_STABLE_STRUCTURE`.

### V2

- [`V2_SOURCE_GATE_PROTOCOL.md`](V2_SOURCE_GATE_PROTOCOL.md)
- [`v2_source_results/V2_SOURCE_GATE_RESULT.json`](v2_source_results/V2_SOURCE_GATE_RESULT.json)
- [`v2_source_results/SOURCE_RESULT_STATUS.json`](v2_source_results/SOURCE_RESULT_STATUS.json)

The null-calibrated source stability gate passed, but the original 12-person signature subset had insufficient index-profile support.

Formal result: `INSUFFICIENT_SOURCE_SIGNATURE`.

### V3

- [`V3_SOURCE_PROTOCOL.md`](V3_SOURCE_PROTOCOL.md)
- [`v3_source_split/V3_SPLIT_STATUS.json`](v3_source_split/V3_SPLIT_STATUS.json)
- [`v3_source_results/source_signature.csv`](v3_source_results/source_signature.csv)
- [`v3_source_results/source_reference.csv`](v3_source_results/source_reference.csv)
- [`v3_source_results/SOURCE_RESULT_STATUS.json`](v3_source_results/SOURCE_RESULT_STATUS.json)

Formal result: `SOURCE_REFERENCE_PASS`.

## Target execution frozen before labels

- [`stage2_prelabel_results/GSE27262_prelabel_assignments.csv`](stage2_prelabel_results/GSE27262_prelabel_assignments.csv)
- [`stage2_prelabel_results/GSE32863_prelabel_assignments.csv`](stage2_prelabel_results/GSE32863_prelabel_assignments.csv)
- [`stage2_prelabel_results/authentic_transfer.csv`](stage2_prelabel_results/authentic_transfer.csv)
- [`stage2_prelabel_results/perturbation_B.csv`](stage2_prelabel_results/perturbation_B.csv)
- [`stage2_prelabel_results/perturbation_C.csv`](stage2_prelabel_results/perturbation_C.csv)
- [`stage2_prelabel_results/perturbation_D.csv`](stage2_prelabel_results/perturbation_D.csv)
- [`stage2_prelabel_results/conditional_null.csv`](stage2_prelabel_results/conditional_null.csv)
- [`stage2_prelabel_results/false_reassurance_summary.csv`](stage2_prelabel_results/false_reassurance_summary.csv)
- [`stage2_prelabel_results/STAGE2_PRELABEL_SHA256.csv`](stage2_prelabel_results/STAGE2_PRELABEL_SHA256.csv)

Both authentic targets were `PASS` before label access.

## Preserved implementation failures

The repository intentionally retains these records:

- [`STAGE2_RUNNER_REPAIR_2026-09-07.md`](STAGE2_RUNNER_REPAIR_2026-09-07.md)
- [`STAGE2_PARTIAL_OUTPUT_CORRECTION_2026-09-07.md`](STAGE2_PARTIAL_OUTPUT_CORRECTION_2026-09-07.md)
- [`stage2_partial_archive_after_failed_runs/`](stage2_partial_archive_after_failed_runs/)

The partial-output archive is part of the audit trail and must not be replaced by the final rerun.

## Post-freeze labels

- [`stage3_label_results/label_contingency.csv`](stage3_label_results/label_contingency.csv)
- [`stage3_label_results/target_assignments_with_labels.csv`](stage3_label_results/target_assignments_with_labels.csv)
- [`stage3_label_results/STAGE3_RESULTS_SHA256.csv`](stage3_label_results/STAGE3_RESULTS_SHA256.csv)

Labels are descriptive/secondary and did not alter the structural endpoint.

## Figures

- [`figures/structural_retention_overview.svg`](figures/structural_retention_overview.svg)
- [`figures/signature_damage_pass_fraction.svg`](figures/signature_damage_pass_fraction.svg)
- [`figures/core_dropout_assigned_fraction.svg`](figures/core_dropout_assigned_fraction.svg)
- [`figures/geometry_perturbation_ari.svg`](figures/geometry_perturbation_ari.svg)

Regenerate with:

```bash
python scripts/generate_preliminary_figures.py
```
