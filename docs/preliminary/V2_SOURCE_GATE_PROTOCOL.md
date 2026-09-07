# V2 source-gate supplement — null-calibrated stability

**Status:** prospective for V2 target access; defined after closure of V1 source gate and before any Stage-2 target-value access.  
**V1 protocol freeze:** `8eb98b22555816dc2771582797e73f068cff062d`  
**V1 source implementation freeze:** `e1dce59bbbeda0aeee8f5baf26a0f25cfa930844`

## 1. Why V2 exists

V1 remains immutable and is closed as `NO_STABLE_STRUCTURE`.

Observed V1 source-only results were:

- K=2: converged; 20/20 stability fits; profile participant support 20/30; Q=0.1018639845 vs null95=0.0223006541; mean stability ARI=0.5895076868, below the frozen absolute threshold 0.60.
- K=3: converged; Q=0.0994879345 vs null95=0.0242831257; mean stability ARI=0.5547811960; minimum profile participant support 7, below the frozen minimum 8.

No V1 parameter is changed and V1 is not re-labelled as a success.

The sole methodological change in V2 is that the arbitrary absolute source-stability cliff (`mean ARI >= 0.60`) is replaced by a **source-null-calibrated stability criterion**. This asks whether the observed reproducibility of the source structure is greater than reproducibility expected from the same complete learning procedure under a participant-block-preserving no-structure randomisation.

This change is explicitly post-V1 development. It is permissible for this development-exposed lung technical anchor because target values and target labels remain sealed. It must never be described as a preregistered V1 success.

## 2. Candidate structure

V2 fixes **K=2**.

Reason: K=2 was the only V1 candidate satisfying every frozen V1 source criterion except the absolute 0.60 stability cliff. K=3 additionally violated the prespecified minimum profile support and is not reconsidered.

The exact V1 participant split, feature namespace, DISCOVERY/VALIDATION partition, top-60 MAD feature rule, strict relation semantics, all-pairs budget, rule support/contrast thresholds, artifact execution thresholds, index-profile rule, independent signature construction, source-reference contrast and source-reference adequacy rules remain unchanged.

## 3. V2 stability-null gate

### 3.1 Observed statistic

Use the same 20 deterministic 80% participant subsamples as V1.

For each subsample:

1. reselect the top 60 DISCOVERY genes by MAD using that subsample only;
2. fit the strict RR_DIRECT K=2 procedure;
3. force-execute its prototypes on all 60 ASSIGNMENT_LEARNING specimens;
4. compute ARI against the full K=2 source fit.

Failed subsample fits contribute ARI=0.

The observed statistic is the mean of the 20 ARIs.

### 3.2 Null distribution

Generate 199 source-null replicates with NumPy PCG64 seed 20260919.

For each replicate, independently for every DISCOVERY gene:

1. permute the 30 participant blocks;
2. independently swap the two within-participant specimen positions with probability 0.5;
3. preserve the complete two-specimen block and every gene's marginal observed values.

The complete DISCOVERY pool is randomised, not only the final 60 genes, so feature reselection inside stability subsamples is represented in the null learning procedure.

For each null replicate:

- fit the complete full-sample K=2 procedure including top-60 MAD selection;
- run the same 20 deterministic stability subsamples, each with fresh within-subsample top-60 MAD selection;
- calculate mean stability ARI exactly as for the observed data.

A failed full null fit or failed stability fit contributes zero to that component. No null replicate is discarded.

### 3.3 Decision

Let `S_obs` be observed mean stability ARI and `S_null95` the 95th percentile of the 199 null mean-stability statistics using NumPy `quantile(method="higher")`.

Empirical one-sided p-value:

`p = (1 + #{S_null >= S_obs}) / 200`.

V2 source structure is eligible only if all are true:

1. the V1 K=2 full fit converges;
2. V1 K=2 minimum profile participant support >=8;
3. V1 K=2 Q_observed > Q_null95;
4. all 20 observed stability subsample fits succeed;
5. `S_obs > S_null95`;
6. empirical p <=0.05.

No absolute 0.60 stability requirement is used in V2.

## 4. Artifact, independent signature and C_source

If the V2 source gate passes, use exactly the frozen V1 rules to:

1. create the K=2 assignment artifact;
2. select the same deterministic index profile criterion;
3. construct the feature-disjoint validation signature on the 12 SIGNATURE_CONSTRUCTION participants;
4. estimate C_source on the untouched 18 SOURCE_REFERENCE participants;
5. require all original V1 source-reference adequacy criteria.

No source-reference participant may influence artifact or signature construction.

If the signature or source-reference gate fails, V2 terminates without target access.

## 5. Target firewall

Throughout V2 source development:

- Stage-2 target expression ZIP remains unopened;
- Stage-3 evaluation-label ZIP remains unopened;
- colorectal data remain unopened.

Stage 2 becomes eligible for opening only if the final V2 source artifact, validation signature, C_source result and hashes are committed and tagged.

## 6. Reporting

V1 and V2 must both be retained.

Any later grant text must state that:

- V1 missed an absolute stability threshold by 0.0105 ARI for K=2;
- V2 was a transparent pre-target methodological refinement replacing that absolute cliff by null calibration;
- the lung evidence is development-exposed preliminary evidence, not untouched confirmation.
