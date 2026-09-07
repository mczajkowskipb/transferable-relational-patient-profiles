# V3 source supplement — assignment-stratified independent validation split

**Status:** prospectively defined after V2 source closure and before any Stage-2 target-value access.  
**V1 protocol freeze:** `8eb98b22555816dc2771582797e73f068cff062d`  
**V2 source-gate freeze:** `f3e3f62a2eaf91a6f284806e0b0cb048ffb6b74e`  
**V2 source outcome:** `1c70e6ee390c214405e050e7c202b5a8a9a7c472`

## 1. Why V3 exists

V1 remains closed as `NO_STABLE_STRUCTURE`.

V2 replaced the arbitrary absolute stability cliff by a source-null-calibrated gate and passed strongly:

- observed mean stability ARI = 0.5895076868;
- null 95th percentile = 0.3201601421;
- empirical one-sided p = 0.005;
- K=2 Q = 0.1018639845 versus Q-null95 = 0.0223006541;
- minimum K=2 profile participant support = 20.

V2 then stopped **before validation-signature selection** because the original unstratified 12-person SIGNATURE_CONSTRUCTION subset contained only 5 participants contributing an accepted P0/index specimen, versus the frozen minimum of 6. Comparator support was 12.

V2 remains `INSUFFICIENT_SOURCE_SIGNATURE`; it is not relabelled as success.

V3 addresses only the sampling-design problem exposed by V2: an unstratified random holdout can make an otherwise supported frozen profile unavailable for independent signature construction. It does **not** lower the signature threshold, change the artifact, change K, alter assignment thresholds, use validation-gene values to choose the split, or use target data.

## 2. Frozen assignment artifact

V3 reuses without refitting the exact V2 K=2 artifact:

`docs/preliminary/v2_source_results/SOURCE_ASSIGNMENT_ARTIFACT_V2.json`

Expected artifact SHA-256:

`9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b`

The frozen index profile is `P0`.

The original 30 ASSIGNMENT_LEARNING participants remain unchanged and never enter signature/reference allocation.

## 3. Assignment-only stratification of the remaining 30 participants

The original 12 SIGNATURE_CONSTRUCTION and 18 SOURCE_REFERENCE participants are pooled into one 30-participant holdout pool **only for reallocation of roles**.

Before accessing any VALIDATION-pool gene value for V3 signature construction:

1. execute the frozen V2 artifact on the holdout participants;
2. use only frozen assignment outputs derived from artifact relations;
3. for each participant record:
   - contributes_index: at least one accepted P0 specimen;
   - contributes_comparator: at least one accepted non-P0 specimen;
   - assigned_participant: at least one accepted specimen;
   - n_assigned_specimens: 0, 1 or 2.

No validation relation, validation-gene MAD, validation contrast, phenotype label, target value or colorectal datum may enter this allocation.

## 4. Exact allocation rule

Retain the original role sizes:

- 12 participants for `SIGNATURE_CONSTRUCTION_V3`;
- 18 participants for `SOURCE_REFERENCE_V3`.

For participant `p`, define deterministic rank:

`SHA256("RPPV3_STRAT_SPLIT|<participant_id>")`.

Among all 12-person signature subsets satisfying the constraints below, choose the subset with:

1. minimum sum of zero-based SHA ranks;
2. tie-break: lexicographically smallest tuple of selected participant IDs.

Required constraints:

### Signature subset
- at least 6 participants contribute an accepted index/P0 specimen;
- at least 6 participants contribute an accepted comparator/P1 specimen.

### Reference subset
- at least 8 participants contribute an accepted index/P0 specimen;
- at least 8 participants contribute an accepted comparator/P1 specimen;
- at least 13/18 participants have at least one accepted specimen (>=0.70 participant assignment coverage);
- at least 26/36 specimens are accepted (>=0.70 specimen assignment coverage).

Mixed participants may count toward both profile-support requirements because C_source is defined by participant-by-group contributions.

If no 12/18 allocation satisfies all constraints, V3 terminates as `NO_FEASIBLE_ASSIGNMENT_STRATIFIED_SPLIT`. Thresholds are not relaxed.

The exact selected split and its SHA-256 are committed and tagged **before** V3 validation-signature construction.

## 5. Independent validation signature

After the exact V3 split is committed:

- use only `SIGNATURE_CONSTRUCTION_V3` participants and VALIDATION genes;
- use the unchanged V1 rules:
  - top 60 VALIDATION genes by MAD in the signature subset;
  - all 1,770 strict pairs;
  - orient toward positive participant-balanced index-minus-comparator contrast;
  - index support >=0.70;
  - observed contrast >=0.20;
  - positive contrast in >=80% of 500 participant-block bootstrap draws;
  - rank by bootstrap-median contrast, then observed contrast, then lexical pair/orientation;
  - at most 10 relations;
  - each gene reused at most twice;
  - at least 5 relations.

No `SOURCE_REFERENCE_V3` participant may influence signature construction.

## 6. V3 source-reference contrast

Use the frozen artifact and the frozen V3 validation signature on `SOURCE_REFERENCE_V3`.

`C_source` remains:

mean participant-level validation score in index profile minus mean participant-level validation score in comparator.

Unchanged adequacy rules:

- specimen assignment coverage >=0.70;
- participant assignment coverage >=0.70;
- at least 8 index contributors;
- at least 8 comparator contributors;
- C_source >=0.15;
- one-sided 95% participant-bootstrap lower bound >0.

Use 2,000 participant-block bootstrap draws.

## 7. Interpretation

The V3 source-reference estimate is conditional on an assignment-stratified holdout allocation. Because allocation uses assignment outputs, not validation values, it preserves independence of the validation feature block from role assignment conditional on the frozen artifact, but it is not an unconditional random-sample estimate of source prevalence.

This is acceptable for a development-exposed technical-anchor pilot. It must not be described as untouched confirmation.

## 8. Firewall

Until V3 source outcome is committed and tagged:

- Stage-2 target expression values remain unopened;
- Stage-3 evaluation labels remain unopened;
- colorectal data remain unopened.

If V3 reaches `SOURCE_REFERENCE_PASS`, Stage 2 becomes eligible for opening only after the final artifact/signature/C_source package is committed and tagged.
