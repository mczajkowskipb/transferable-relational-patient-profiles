# Stage-3 label-evaluation implementation freeze

**Date:** 2026-09-07  
**Required prelabel result commit:** `f633d7cc4ae3c26e794798f3e4ec3194ad7e8aa0`  
**Required prelabel tag:** `preliminary-stage2-prelabel-freeze-2026-09-07`

This implementation is frozen and tagged before the Stage-3 evaluation-label
ZIP is opened.

Stage 3 is strictly descriptive. No label may alter:

- the frozen assignment artifact;
- K or profile definitions;
- the V3 validation signature;
- C_source;
- source/target thresholds;
- PASS/FAIL/INSUFFICIENT_SUPPORT decisions;
- B/C/D perturbation results;
- conditional false-reassurance results;
- any prelabel target assignment.

## 1. Stage-3 payload

Allowed payload:

`POST_ASSIGNMENT_LABELS_2026-09-07.zip`

Expected SHA-256:

`72fff8d5fe2c1785018915bf1c567d58b20184d71d5122e2fdebb001d344d30a`

Only the two target evaluation-label files are used:

- `GSE27262_evaluation_labels.csv`
- `GSE32863_evaluation_labels.csv`

The GSE19804 label file present in the same package is not used for the primary
Stage-3 target evaluation.

## 2. Label import

For each target CSV:

- sample identifier column is detected case-insensitively from:
  `gsm`, `sample_id`, `sample`, `geo_accession`;
- evaluation label column is detected case-insensitively from:
  `evaluation_label`, `label`, `class`, `phenotype`, `tissue`, `status`;
- if detection is ambiguous or absent, stop rather than guess;
- labels are stripped of surrounding whitespace but otherwise preserved exactly;
- no biological recoding or relabelling is permitted;
- duplicate GSM rows or missing labels are errors;
- the label GSM set must exactly equal the frozen prelabel assignment GSM set.

## 3. Primary descriptive agreement

Report at specimen level.

### Accepted-only
Population: rows with frozen `assignment == ASSIGNED`.

Metrics:
- Adjusted Rand Index (ARI);
- Normalized Mutual Information (NMI), arithmetic normalization:
  `MI / ((H(profile)+H(label))/2)`;
- exact contingency counts profile × evaluation label;
- accepted fraction and n.

### Forced
Population: all target specimens for which at least one frozen profile score is
finite.

The forced profile is the profile with highest already-frozen score; ties are
resolved by lexical profile ID. No score is recomputed after label access.

Report:
- ARI;
- NMI;
- exact contingency counts;
- n.

These metrics are secondary/descriptive and cannot rescue or overturn the
structural endpoint.

## 4. Post-hoc binary mapped accuracy

If and only if the target has exactly two evaluation labels and the artifact has
exactly two profiles, report descriptive optimal one-to-one profile↔label
mapping for the accepted and forced populations separately.

The mapping maximizes the number of correctly mapped specimens. Ties are
resolved lexically by the serialized mapping.

Report:
- mapped accuracy;
- mapped balanced accuracy;
- selected mapping.

This is explicitly post-hoc descriptive reporting, not a frozen classifier
accuracy endpoint.

## 5. Final evidence summary

After label evaluation, create:

- `label_agreement.csv`
- `label_contingency.csv`
- `target_assignments_with_labels.csv`
- `FINAL_PRELIMINARY_EVIDENCE_TABLE.csv`
- `FINAL_PRELIMINARY_EVIDENCE_SUMMARY.md`
- `STAGE3_LABEL_STATUS.json`
- `STAGE3_RESULTS_SHA256.csv`

The Markdown summary must retain the development chronology:

- V1 closed as `NO_STABLE_STRUCTURE`;
- V2 passed null-calibrated source stability but stopped for signature support;
- V3 used a transparently frozen assignment-stratified holdout allocation and
  obtained `SOURCE_REFERENCE_PASS`;
- both authentic lung targets were then executed and frozen before labels;
- label agreement is reported only after the Stage-2 prelabel freeze.

The lung evidence may be described as development-exposed preliminary evidence,
not untouched confirmation and not proof of a universal biological subtype.
