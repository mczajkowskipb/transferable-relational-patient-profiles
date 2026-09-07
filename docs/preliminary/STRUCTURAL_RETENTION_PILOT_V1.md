# Structural-retention pilot V1 — FROZEN PROTOCOL CANDIDATE

**Protocol version:** V1  
**Prepared:** 2026-09-07  
**Preflight parent commit:** `9b090b62d0d7024c86344e7d031b6fd5cc840e98`  
**Target expression values must remain sealed until this protocol, `V1_CONFIG.json`, the source split, Stage-1 input manifest and sealed-payload hashes are committed and tagged. Evaluation labels remain sealed until frozen target assignments are written and hashed.**

## 1. Scientific purpose

This preliminary pilot tests whether frozen assignment and an independently constructed molecular validation contrast measure a different inferential object from partition agreement and whether the proposed endpoint behaves correctly under controlled mechanisms of transport, structural failure and non-executability.

The pilot is not designed to establish a universal biological subtype, clinical utility or general transportability. Positive structural retention supports persistence of one independently characterised source molecular property after unchanged source-defined assignment. Biological interpretation remains secondary.

The four prespecified real-data/mechanistic states are:

- **A — authentic frozen execution:** unchanged source artifact on both lung targets;
- **B — order-preserving geometry change:** value-space geometry may change while strict within-sample order relations remain unchanged;
- **C — validation-specific damage:** assignment information is held fixed while reserved validation information is progressively disrupted;
- **D — assignment-core feature loss:** executability is progressively removed and the system must abstain rather than substitute target-specific features.

All V1 outcomes are retained. No target, dose, seed, K, profile, threshold or signature may be replaced because its result is unfavourable.

## 2. Historical transition and dataset roles

`docs/prospective/EXTERNAL_VALIDATION_FREEZE_v2.md` dated 2026-08-25 remains an immutable historical record.

At the commit/tag that freezes this V1 protocol, **GSE27262 and GSE32863 are prospectively re-designated as development-exposed preliminary/technical-anchor targets**. They may no longer be described as untouched confirmatory targets once Stage 2 is opened.

Dataset roles:

- source: **GSE19804**, GPL570;
- preliminary target 1: **GSE27262**, GPL570;
- preliminary target 2: **GSE32863**, GPL6884;
- excluded from V1 development and outcome access: **GSE39582, GSE14333, GSE33113**.

The colorectal family remains reserved for the principal future heterogeneous-disease confirmation.

## 3. Three-stage firewall

### Stage 1 — allowed before protocol freeze

Permitted:

- GSE19804 expression-only source matrix;
- phenotype-blind participant maps for all three lung datasets;
- GPL570/GPL6884 annotations;
- provenance, checksums and metadata needed for mapping;
- no target expression values;
- no evaluation labels.

### Stage 2 — allowed only after V1 protocol commit/tag

Open the exact pre-hashed payload containing expression-only GSE27262 and GSE32863 target matrices. Evaluation labels remain sealed.

The complete source artifact and validation signature must already be frozen before authentic target execution. Target values may only be used for execution, structural-retention evaluation and the prespecified B/C/D/null analyses. They cannot modify source fitting, mapping, K, profile selection, validation-signature selection, assignment thresholds or artifact choice.

Before Stage 3, write and SHA-256 hash target assignments, profile scores, margins, executable coverage and structural outputs.

### Stage 3 — allowed only after target assignment hashes exist

Open the separately pre-hashed evaluation-label payload. Labels may then be used only for secondary ARI/NMI or descriptive phenotype comparisons. They cannot cause any V1 rerun or model change.

## 4. Participant unit and source partition

Participant, not specimen, is the independent resampling/splitting unit. Paired specimens from one participant never cross source partitions.

The 60 GSE19804 participants are ordered by the hexadecimal SHA-256 digest of:

`20260907|<participant_id>`

The first 30 are **ASSIGNMENT_LEARNING**, the next 12 **SIGNATURE_CONSTRUCTION**, and the remaining 18 **SOURCE_REFERENCE**. Exact IDs are frozen in `SOURCE_PARTICIPANT_SPLIT_V1.csv`.

No phenotype stratification is used.

This split solves the artifact-level estimand ambiguity:

1. learn and freeze the assignment object on ASSIGNMENT_LEARNING;
2. construct the independent validation signature on SIGNATURE_CONSTRUCTION using the already frozen assignment object;
3. estimate `C_source` on SOURCE_REFERENCE using the same frozen assignment object and the same frozen signature;
4. do **not** refit after source-reference evaluation;
5. execute that identical object on every target.

Cross-fitted or resampling quantities describe the learning procedure; they are not substituted for the artifact-level `C_source` denominator.

## 5. Platform mapping and feature namespace

Primary namespace: **Entrez Gene ID**.

Using GPL570/GPL6884 annotation only:

1. retain a platform probe only when `Gene ID` is exactly one numeric Entrez identifier;
2. discard missing and multi-ID annotations;
3. define the eligible common universe as the intersection of such Gene IDs on GPL570 and GPL6884;
4. when multiple probes map to one eligible Gene ID, aggregate them within each specimen by the median;
5. no target expression value or evaluation label may influence mapping.

Feature-pool separation is metadata-deterministic. For each eligible Gene ID compute:

`SHA256("RPPV1_POOL|<gene_id>")`.

If the first digest byte is even, assign the gene to **DISCOVERY**; if odd, assign it to **VALIDATION**.

Validation genes therefore cannot participate in discovery screening, K selection, RR_DIRECT induction, index-profile choice or assignment. Feature disjointness prevents direct information reuse; it is not interpreted as statistical independence from shared batch, cell composition or continuous biology.

## 6. Source assignment learning

Only ASSIGNMENT_LEARNING participants and DISCOVERY genes may influence assignment learning.

For each candidate `K ∈ {2,3}`:

- rank DISCOVERY genes by MAD on ASSIGNMENT_LEARNING specimens only;
- retain the top 60; ties are broken by numeric Entrez ID;
- use all 1,770 unordered pairs among those genes;
- strict relation semantics: `A>B` and `A<B`; ties satisfy neither;
- RR_DIRECT maximum 10 relations/profile;
- minimum rule support 0.80;
- minimum rule contrast 0.10;
- maximum 50 iterations;
- the legacy “no eligible rule → fall back to arbitrary candidates” path is **not allowed for the V1 primary fit**. If a profile lacks an eligible rule set, that K is ineligible.

Frozen artifact execution thresholds are fixed before target access:

- minimum score = 0.60;
- minimum winning margin = 0.05;
- minimum executable relation coverage = 0.80.

These thresholds are not target-calibrated.

### 6.1 Source-only K gate

For each K compute two prespecified source-only diagnostics.

**Resampling stability.** Run 20 deterministic 80% participant subsamples (24/30 participants). Fit the complete K-specific assignment procedure on each subsample and force-execute it on all ASSIGNMENT_LEARNING specimens. Mean ARI against the full-assignment fit is the stability statistic. A failed resample contributes zero stability. At least 18/20 resamples must fit successfully and mean stability must be ≥0.60.

**Null-calibrated relational separation.** For the full fit, form the binary candidate-relation matrix and compute:

`Q_K = mean(within-cluster relation similarity) − mean(between-cluster relation similarity)`

where relation similarity is one minus normalized Hamming distance. Generate 199 source-null matrices by independently permuting each discovery gene across participant blocks, with an independent random within-pair swap for each gene/block mapping. This preserves each gene's marginal values and participant-block size while disrupting multigene assignment structure without labels. Refit the K-specific pipeline for every null. The observed `Q_K` must exceed the 95th percentile of its null distribution.

Each final profile must additionally occur in at least 8 distinct assignment participants.

If neither K passes, V1 source status is `NO_STABLE_STRUCTURE` and authentic target retention cannot be declared PASS.

If both pass, select K by:

1. largest `Q_K − null95_K`;
2. then larger mean resampling stability;
3. then smaller K.

No evaluation label is used.

### 6.2 Index profile

For each final profile compute the mean of `support × positive contrast` across its frozen rules. The profile with the largest value is the **index profile**; ties are broken by lexical profile ID. The comparator is all other accepted frozen profiles combined with equal participant weighting at the contrast stage. The index profile cannot be replaced after target access.

## 7. Independent validation signature

Only SIGNATURE_CONSTRUCTION participants and VALIDATION genes may select the signature.

1. Execute the already frozen assignment artifact on signature participants.
2. Require at least 6 distinct participants contributing an accepted index-profile specimen and at least 6 contributing an accepted comparator specimen.
3. Rank VALIDATION genes by MAD within the signature subset and retain 60; ties use numeric Entrez ID.
4. Evaluate all 1,770 strict pairs.
5. Orient each relation so its participant-balanced index-minus-comparator contrast is positive.
6. A candidate must have:
   - index-group satisfaction ≥0.70;
   - observed contrast ≥0.20;
   - positive contrast in ≥80% of 500 participant-block bootstrap resamples.
7. Rank by bootstrap-median contrast, then observed contrast, then lexical pair/orientation.
8. Greedily retain at most 10 relations, with each gene used by at most two selected relations.
9. At least five relations are required.

No SOURCE_REFERENCE observation can influence signature construction.

If fewer than five qualifying feature-disjoint relations exist, status is `INSUFFICIENT_SOURCE_SIGNATURE`. Relation-disjoint but feature-overlapping signatures may be reported only as sensitivity analyses and cannot rescue the V1 primary endpoint.

Signature relations receive equal weight. A specimen's validation score `v(x)` is the mean strict satisfaction across the complete frozen signature. The primary structural analysis does not silently drop unavailable signature relations.

## 8. Source-reference contrast

Execute the frozen assignment artifact on SOURCE_REFERENCE participants. Apply the frozen validation signature without refitting.

For each participant and each group separately, average `v(x)` over that participant's accepted specimens assigned to that group. Then define:

`C_source = mean(participant-level v in index profile) − mean(participant-level v in comparator)`.

This equal-participant construction prevents participants contributing two specimens to one group from receiving twice the weight.

Source-reference adequacy requires:

- specimen assignment coverage ≥0.70;
- participant assignment coverage (at least one accepted specimen) ≥0.70;
- at least 8 distinct participants contributing to the index group;
- at least 8 contributing to the comparator;
- `C_source ≥ 0.15`;
- one-sided 95% participant-bootstrap lower bound for `C_source` > 0.

Use 2,000 participant-block bootstrap resamples, seed 20260914.

If source-reference adequacy fails, authentic target structural status is `INSUFFICIENT_SOURCE_REFERENCE`. Target execution and diagnostics may still be reported, but they cannot be reclassified as a positive structural-retention result.

## 9. Target contrast and primary structural decision

For each authentic target, before evaluation labels are opened:

1. execute the unchanged artifact;
2. compute validation score using the unchanged signature;
3. aggregate within participant exactly as for source reference;
4. calculate `C_target`;
5. calculate `ΔC = C_target − C_source`;
6. if `C_source` passes its floor, report descriptive `R = C_target / C_source`.

Adequate target support requires:

- specimen assignment coverage ≥0.70;
- participant assignment coverage ≥0.70;
- at least 8 distinct index-profile participants;
- at least 8 distinct comparator participants;
- validation evaluability ≥0.90 among accepted specimens.

Define the non-inferiority-style retention quantity:

`D = C_target − 0.50*C_source`.

For uncertainty use 2,000 participant-block bootstrap resamples. Source-reference bootstrap draws are reused across targets/method comparisons; target participants are resampled independently. Seed 20260915.

### Decision

**PASS**
- support adequate;
- `C_target ≥ 0.10`; and
- one-sided 95% lower bootstrap bound of `D` > 0.

**FAIL**
- support adequate; and either
- one-sided 95% upper bound of `C_target` ≤ 0; or
- one-sided 95% upper bound of `D` < 0.

**INSUFFICIENT_SUPPORT**
- support is inadequate; or
- neither PASS nor FAIL criterion is established.

Absence of PASS is not automatically evidence of absence of retention.

Both authentic lung targets are reported. Neither may replace the other or be replaced by another dataset.

## 10. Mechanistic module B — order-preserving geometry change

On each Stage-2 target apply:

`x'_ij = exp(t*z_i) * x_ij + t*s_source*u_i`

for `t ∈ {0, 0.25, 0.5, 1, 2}`.

- `z_i` and `u_i` are fixed standard-normal draws generated in lexically sorted GSM order using NumPy `PCG64(20260916)`;
- `s_source` is the median within-specimen SD across the mapped common-gene source matrix and is computed source-only;
- `exp(t*z_i)>0`, therefore all finite non-tied within-sample order relations are mathematically preserved.

Required implementation check: the strict relation flip rate for all frozen assignment/signature relations must be exactly zero apart from explicitly reported numerical ties.

Report separately:

1. ARI between frozen RPP assignments and external evaluation labels after Stage 3;
2. ARI between frozen RPP assignments and deterministic target value-space clustering;
3. ARI between target value-space clustering before and after perturbation.

Target value-space clustering is a diagnostic only: deterministic Euclidean k-medoids, K fixed to source-selected K, no target tuning. Quantity 1 must remain unchanged if frozen assignments and reporting subset are unchanged. Quantities 2/3 may change. This module does **not** explain the historical label-ARI; it demonstrates separation of estimands.

A rank-space frozen comparator should be invariant under B as an explanatory control; this is not claimed as RPP-specific superiority.

## 11. Mechanistic module C — validation-specific damage

Assignment-core genes are untouched.

Doses: `f ∈ {0, .10, .25, .50, .75, 1}`.

Unique frozen validation-signature genes are ordered by `SHA256("C_DAMAGE|<gene_id>")`; the first `ceil(f*m)` are disrupted. For each disrupted gene, independently permute participant blocks within pair-status class and independently swap the two specimen positions inside paired blocks with probability 0.5. Use PCG64 seed 20260917 with deterministic gene/replicate substreams. This uses no phenotype labels.

Use 20 fixed perturbation replicates per nonzero dose.

Required invariant: assignment outputs for the unmodified assignment block must be byte-identical to intact-target assignments at every C dose. If not, the C experiment is invalid.

Report `C_target`, `D`, R where defined, and decision state versus dose. Monotonic decline is not required for every realization; all doses/replicates are retained.

## 12. Mechanistic module D — assignment-core feature loss

Unique frozen assignment-core genes are ordered by:

`SHA256("D_DROPOUT|<gene_id>")`.

Doses: `f ∈ {0, .10, .25, .50, .75, 1}`. Remove the first `ceil(f*m)` unique core genes at each nested dose. Do not search for substitutes.

Report:

- per-profile executable coverage;
- accepted/unassigned fractions and reasons;
- scores and margins;
- structural analysis only where support remains adequate.

Per-profile executable coverage must not increase under nested deletion. Score, margin and overall abstention need not be monotone at every intermediate dose. At `f=1`, every sample must be `UNASSIGNED` because all assignment-core genes are unavailable; failure is an implementation error.

## 13. Conditional false reassurance

For each intact Stage-2 target, use 999 full validation-block randomisations following the label-free C randomisation mechanism, seed 20260918. Assignment remains fixed.

For every replicate apply the already frozen structural decision rule. Define:

`false reassurance = number of null replicates classified PASS / 999`.

Report the exact binomial 95% interval as Monte Carlo uncertainty conditional on the observed target matrix and randomisation mechanism.

Endpoint calibration is considered acceptable for this preliminary if the **upper 95% binomial bound is ≤0.10 for each target**. If this criterion fails, an authentic PASS cannot be presented as strong preliminary support; the calibration failure must be reported.

This is conditional endpoint calibration, not unconditional type-I error control for the entire learning procedure.

## 14. Full-source null diagnostic

As a separate source-only diagnostic, generate 99 complete source-null replicates using the same participant-block / within-pair gene randomisation principle and rerun the assignment/signature/reference procedure with the fixed V1 rules. Seed 20260912.

Report the fraction of null replicates that produce an eligible source artifact, valid independent signature and adequate positive source-reference contrast. This diagnostic cannot rescue or redefine any authentic target result.

## 15. Explanatory frozen baselines

Secondary explanatory controls:

- frozen value centroid;
- frozen rank centroid.

They use source-only RR_DIRECT pseudo-labels from ASSIGNMENT_LEARNING, the same 60 discovery-gene budget and source-only calibration. They are frozen before target execution and cannot tune on target data. Their purpose is to separate benefits of relational execution from benefits of freezing a source-defined partition. They do not rescue the primary RPP endpoint.

If implementation time prevents a baseline, report it as missing rather than modifying the primary protocol.

## 16. Secondary labels and ARI/NMI

ARI/NMI are retained and reported; they are not redefined.

ARI can compare frozen assignments directly with existing target evaluation labels and does **not** require target reclustering. It answers a different question from structural retention.

After Stage 3 report, separately:

- ARI/NMI of frozen accepted/forced assignments versus external labels, with the reporting population explicit;
- value-space clustering agreement diagnostics from B;
- no label may alter an earlier frozen result.

The historical median direct-RPP label ARI is not evidence that the new structural endpoint will succeed and is not “explained away” by Module B.

## 17. Nuisance and construct-validity controls

Feature-disjoint validation is necessary for non-circularity but does not prove statistical independence. Where metadata are available after freeze, evaluate whether retained validation contrast is explainable by known technical nuisance, composition or a continuous source axis. Such analyses are secondary and cannot change the frozen structural result.

Structural retention supports persistence of the specified independent molecular contrast. It does not by itself prove a discrete universal biological subtype.

## 18. Required outputs

No empirical result file may exist before the corresponding access stage.

Minimum outputs:

`preliminary/results/`
- `source_fit_summary.csv`
- `source_signature.csv`
- `source_reference.csv`
- `authentic_transfer.csv`
- `perturbation_B.csv`
- `perturbation_C.csv`
- `perturbation_D.csv`
- `conditional_null.csv`
- `false_reassurance_summary.csv`
- `baseline_summary.csv`
- `pilot_summary.csv`

`preliminary/artifacts/`
- frozen assignment artifact JSON + SHA-256
- frozen validation signature JSON + SHA-256
- source mapping/config/software hashes
- prelabel target assignment CSVs + SHA-256

`preliminary/figures/`
- `pilot_states.svg`
- `ari_vs_structural_retention.svg`
- `feature_loss_abstention.svg`

`preliminary/docs/`
- `RESULTS_V1.md`
- `INTERPRETATION_V1.md`
- `LIMITATIONS_V1.md`
- `REPRODUCIBILITY.md`
- `GRANT_READY_SUMMARY.md`

Every figure and grant-ready number must be generated from frozen machine-readable outputs.

## 19. V1 integrity rules

V1 is immutable after target values are opened.

If V1 reveals a methodological defect:

1. retain all V1 outputs;
2. diagnose the defect;
3. define V2 prospectively;
4. commit/hash V2 before new result-generating execution;
5. never rewrite V1 history.

No V2 result may be presented as if it were the original V1.

## 20. Interpretation ceiling

A successful V1 may support:

- feasibility of a non-circular artifact-level structural-retention endpoint;
- empirical separation of assignment executability from independent structural retention;
- correct abstention under missing assignment features;
- controlled false-reassurance behaviour;
- development-exposed real-data evidence in the lung technical anchor.

It cannot establish:

- universal transportability;
- a general biological subtype identity;
- H2 family-level predictive validity;
- colorectal confirmation;
- clinical usefulness.

Those remain future project questions.
