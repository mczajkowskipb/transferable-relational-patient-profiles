# Stage-2 prelabel execution implementation freeze

**Prepared:** 2026-09-07  
**Required source-package tag:** `preliminary-v3-source-package-freeze-2026-09-07`  
**Required source-package commit:** `6820763d00d879a24acfda0141adb48384680d6d`

This document freezes operational details that were not fully specified in the
earlier scientific protocol. It is committed and tagged **before the Stage-2
target-value ZIP is opened**.

The scientific thresholds and mechanisms remain those of the frozen V1
protocol plus the transparent V2/V3 source-development supplements. No target
label is available at this stage.

## 1. Stage-2 payload

The only Stage-2 payload that may be opened is:

`POST_PROTOCOL_TARGET_VALUES_2026-09-07.zip`

Expected SHA-256:

`94e648dfa9ec9b9cd16a6eaea067f1d823e7d74eb128cc27cfc6de473081e2df`

The Stage-3 evaluation-label payload remains unopened.

## 2. Target mapping

For each target:

- GSE27262 uses GPL570;
- GSE32863 uses GPL6884;
- retain only probes with exactly one numeric Entrez Gene ID;
- aggregate multiple probes per Entrez ID by within-specimen median;
- represent the exact 17,062-gene common universe frozen at Stage 1;
- genes absent from a Series Matrix are represented as missing, never replaced.

No target statistic changes the feature namespace, artifact, signature or
thresholds.

## 3. Frozen objects

Use unchanged:

- assignment artifact SHA-256
  `9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b`;
- validation signature SHA-256
  `64eb3631d6a25f0fdf0deede0e0c2da17c63fb25d82ee4612aaa8bd12470128d`;
- index profile `P0`;
- V3 C_source `0.6053475935828878`.

The source-reference bootstrap vector is regenerated from the frozen V3 split,
artifact and signature using the already frozen seed 20260914 and must
reproduce the frozen V3 source lower bound before target analysis proceeds.

## 4. Target bootstrap substreams

For each target create one fixed 2,000 x N participant bootstrap index matrix
and reuse it for authentic, C, D and conditional-null comparisons.

RNG: NumPy PCG64 via:

- GSE27262: `SeedSequence([20260915, 27262])`;
- GSE32863: `SeedSequence([20260915, 32863])`.

For every row, sample N participant indices with replacement.

This makes target resampling independent between targets and common across
within-target method/perturbation comparisons. The already frozen source
bootstrap vector is paired row-wise with every target bootstrap vector when
forming `D = C_target - 0.5*C_source`.

One-sided bounds use NumPy quantiles:
- lower: q=0.05, `method="lower"`;
- upper: q=0.95, `method="higher"`.

## 5. Validation evaluability

A validation score is evaluable only when all frozen signature relations have
both features finite. No relation is silently dropped. The score is the
equal-weight mean strict satisfaction across all ten relations.

## 6. Participant-balanced target contrast

For each participant and each frozen assignment group separately, average
validation score across accepted, evaluable specimens assigned to that group.
Then average participant-level group values with equal participant weight.

Support/decision thresholds remain exactly frozen:
- specimen assignment coverage >=0.70;
- participant assignment coverage >=0.70;
- >=8 index contributors;
- >=8 comparator contributors;
- validation evaluability >=0.90 among accepted specimens;
- C_target absolute PASS floor 0.10;
- rho=0.50;
- PASS/FAIL/INSUFFICIENT_SUPPORT rules unchanged.

## 7. Mechanistic B implementation

The B transform remains:

`x'_ij = exp(t*z_i)*x_ij + t*s_source*u_i`

with t in {0,.25,.5,1,2}; z/u are generated once per target in lexical GSM
order from PCG64(20260916), and `s_source` is the median within-specimen SD
across all mapped common source genes.

The required frozen assignment/signature strict-relation flip count is zero.

For the target value-space diagnostic only:

1. reconstruct the source-frozen top-60 DISCOVERY genes from the original 30
   ASSIGNMENT_LEARNING participants;
2. freeze source mean/SD on those 60 genes;
3. source-standardize target values with those parameters;
4. cluster with deterministic Euclidean k-medoids, K=2;
5. initialize first medoid as lexical GSM #1 and second as the farthest sample;
6. update medoids by minimum within-cluster total distance, lexical index on ties.

No target hyperparameter is selected.

Report:
- ARI between accepted frozen RPP assignments and value-space clustering;
- ARI between intact and perturbed value-space clusterings;
- strict relation flip rate.
Label ARI is deferred until Stage 3.

## 8. Mechanistic C implementation

Unique signature genes are ordered by the frozen SHA rule
`SHA256("C_DAMAGE|<gene_id>")`.

For each disrupted gene and replicate, participant blocks are permuted within
their blind `pair_status` class. Blocks of different specimen count are never
mixed. Two-specimen source blocks are independently reversed with probability
0.5.

Deterministic substream:

`SeedSequence([20260917, <dataset_numeric>, <replicate>, <gene_id>])`.

The same gene/replicate randomisation is reused across nested doses.

Assignment-core and validation-signature genes are asserted disjoint.
Frozen assignment outputs must be identical at all C doses.

## 9. Mechanistic D implementation

Unique assignment-core genes are ordered by
`SHA256("D_DROPOUT|<gene_id>")`.

At dose f remove the first ceil(f*m) genes from the executable feature table.
No replacement or imputation is permitted. Validation values remain intact.

At f=1 all samples must be UNASSIGNED; otherwise execution stops as an
implementation failure.

## 10. Conditional false reassurance

For each target use 999 full-signature randomisations. Every unique frozen
signature gene is randomised with the same blind participant-block mechanism.

Substream:

`SeedSequence([20260918, <dataset_numeric>, <replicate>, <gene_id>])`.

Assignments remain fixed. Each replicate is evaluated with the full frozen
structural decision rule and the fixed target/source bootstrap draws.

The exact two-sided 95% Clopper-Pearson interval is reported for PASS/999.
Calibration is acceptable only when its upper bound <=0.10.

## 11. Deferred non-gating diagnostics

The 99-replicate full-source pipeline null and explanatory frozen centroid
baselines do not gate Stage-2 access or the primary structural endpoint.
If not completed in this workflow they are explicitly reported as deferred /
missing, as permitted by the frozen protocol; they cannot rescue the primary
result.

## 12. Stage-3 firewall

Before Stage 3:

- authentic target assignments, scores, margins, coverages and validation
  scores are written;
- authentic structural outputs and B/C/D/conditional-null outputs are written;
- a SHA-256 manifest is written;
- all prelabel Stage-2 results are committed and tagged.

Only after that tag exists may the Stage-3 label ZIP be opened.
