# Structural-retention pilot V1 — design draft

**NOT FROZEN. NOT AUTHORISATION TO OPEN TARGET OUTCOMES. NO V1 RESULTS EXIST.**

This document specifies the proposed implementation direction after the repository audit. It is intentionally not named or hashed as a completed freeze. Input checksums, participant mapping, annotation releases, source calibration and numerical operating characteristics must be resolved before a committed V1 protocol. No target outcomes may be inspected to resolve these choices.

## Scientific question and scope

Test whether unchanged assignment and independent molecular contrast distinguish (A) authentic execution, (B) geometry change preserving within-sample order, (C) loss of validation contrast despite executable assignment, and (D) loss of executability. These controlled tests validate endpoint behaviour; they do not establish discrete biological subtypes or the generality of transport.

Use GSE19804 source (reported 60 participants/120 paired specimens), with GSE27262 and GSE32863 as prespecified preliminary targets after input audit and protocol commit. Counts must be verified, not inferred from matrix column ordering. Colorectal GSE39582/GSE14333/GSE33113 expression, labels and outcomes remain sealed.

## Historical transition — required wording in the eventual committed freeze

“The 2026-08-25 EXTERNAL_VALIDATION_FREEZE_v2.md remains an immutable historical record. From the new protocol's actual commit timestamp, GSE27262 and GSE32863 are intentionally designated development-exposed preliminary/technical-anchor targets. After this access they cannot be described as untouched confirmatory evidence. The colorectal family remains sealed and unused. The new structural-retention protocol is committed before generating or inspecting its results.”

Do not backdate this text. In this draft the transition has NOT occurred.

## Estimand: one fixed object, one independent source reference

A cross-fitted mean across different learned profiles/signatures need not estimate the source contrast of the final full-source refit. Label switching and changing K compound the problem. Do not divide a target contrast of one object by a source contrast of another without an explicit common estimand.

Recommended pilot: fixed participant-level source splits for assignment learning, independent signature construction and independent effect estimation. Freeze the assignment model before signature construction. Use the same frozen assignment model and signature for the source reference and all target evaluations; do NOT subsequently refit using reference participants.

Candidate allocation for metadata verification: 30/10/20 participants (assignment/signature/reference), all paired specimens kept together. Use a deterministic shuffled sorted participant list, proposed seed 20260907, without phenotype stratification. Small signature and reference sizes are limitations and may produce insufficient support. Do not try several splits and select the successful one.

Source stability/K diagnostics can use participant-level inner folds on the assignment subset. All feature screening, K selection, profile induction, index selection and tuning relevant to each diagnostic must occur inside its training folds. Such diagnostics describe a learning procedure, not unbiased exact-edge recovery of the final artifact.

This is a proposed resolution of the grant's full-refit/cross-fitting ambiguity. It must be stated explicitly in the final protocol/grant; this draft does not silently redefine the confirmatory endpoint.

## Input contract and preprocessing

Require, for each dataset:
- accession/platform, specimen IDs, participant IDs, file checksum and provenance;
- original expression scale and preprocessing record; no joint source-target normalisation;
- source specimen identity checked independently of outcome labels;
- annotation release/checksum, one-to-many identifier handling, per-sample aggregation and feature namespace;
- a target matrix access firewall, separate from metadata/platform access.

Relations must operate on a declared comparable measurement scale. Avoid feature-wise restandardising the RPP input and then claiming invariance of original abundance order. The single-sample claim applies to model execution after the stated preprocessing.

Freeze a common gene universe using platform annotation only. Do not match Affymetrix probe IDs directly to Illumina probe IDs. Ambiguous mappings must be handled by a deterministic policy established without target values. Within-sample median aggregation is a candidate policy, not a guarantee of platform invariance.

## Assignment learning and primary profile

Reuse repaired RR_DIRECT and the coverage-aware artifact executor. Candidate source-only K set {2,3}; choose by a prespecified null-calibrated source score subject to convergence, minimum support, stability and complexity. No supported K means NO_STABLE_STRUCTURE. A fixed K=2 feasibility sensitivity is allowed only if separately labelled; it cannot rescue the primary gate.

Candidate budgets before source calibration: 60 discovery features, all 1770 pairs among those features, at most 10 relations/profile. This avoids the existing lexicographic truncation of the first 1500 pairs. Reject a profile system if its final rules fail the required support/contrast rather than relying on the legacy fallback. Freeze score, winning margin and coverage thresholds on source-only calibration.

For K>2 select one index profile by a fixed source-only support/stability criterion and deterministic tie-break. No target profile replacement. When comparator is “other assigned profiles”, freeze whether the comparator uses source prevalence weights; otherwise changing mixtures can change C without loss of individual component identities.

## Validation signature

Reserve validation features outside **every feature used to discover the primary assignment structure**, not just genes surviving into the final core. Use a source-only deterministic discovery/validation feature-pool partition if needed. No shared genes, mapped aliases or duplicate probes across pools.

On the signature-construction participants, execute the already frozen assignment model and select oriented validation relations from the reserved pool. Proposed cap: 10 relations, deterministic contrast/support ranking, bounded feature reuse. No source-reference participant may influence this selection. If no supported feature-disjoint signature exists, primary status is INSUFFICIENT_SOURCE_SIGNATURE; relation-disjoint overlap is sensitivity only.

Feature disjointness prevents direct reuse; it is not statistical independence. Shared biological/technical latent factors must be examined with prespecified nuisance controls. A stable continuous latent axis can legitimately retain a contrast without demonstrating a discrete subtype.

## Contrast and uncertainty

For fixed signature h, compute per-specimen mean strict relation satisfaction v(x). C is mean v among accepted index-profile members minus mean v among the frozen comparator. Missing signature observations cannot silently change the averaged relation set: primary analysis requires the full fixed signature, with validation evaluability and excluded observations reported separately from assignment coverage.

Estimate C_source on the untouched source-reference participants; C_target after target assignments have been recorded and hashed. Report participant and specimen counts separately, including distinct participants in profile/comparator and overlap from paired specimens.

R=C_target/C_source is descriptive, defined only above a frozen source evidence floor. Report DeltaC and both contrasts with uncertainty. Bootstrap participant blocks; keep paired specimens together. Use the SAME source-bootstrap draw when comparing both targets or methods to avoid pretending their source-reference uncertainty is independent.

Candidate decision form (final numerical constants require pre-target operating-characteristic verification):
- PASS only with adequate assignment/validation support, positive target contrast and a one-sided lower confidence bound for C_target − rho*C_source above zero, plus any frozen absolute target-effect floor;
- FAIL when adequate precision supports degradation/nonpositive contrast under the frozen failure rule;
- INSUFFICIENT_SUPPORT when count/precision requirements fail or the interval does not establish either outcome.

Failure to prove retention is not automatically evidence that retention is absent. Report the three outcomes without excluding unsuccessful directions. Proposed bootstrap seed 20260908 and 2000 resamples; no seed search. Exact support thresholds, alpha, rho, source floor and handling of degenerate resamples must be frozen before target access.

## Conditions

| Module | Proposed perturbation | Required interpretation |
|---|---|---|
| A | Unchanged source artifact on both authentic targets | Report both, including failures; development-exposed evidence |
| B | x'_ij=exp(t*z_i)*x_ij+t*s_source*u_i, fixed z/u, t in {0,.25,.5,1,2} | Preserves order mathematically; must preserve RPP and rank-based execution up to declared numeric checks |
| C | Perturb only validation block at fractions {0,.10,.25,.50,.75,1} using prespecified participant-block randomisation | Assignment unchanged; validation contrast may degrade; no dose chosen after results |
| D | Nested deletion of assignment-core features at fractions {0,.10,.25,.50,.75,1} | Per-profile executable coverage cannot increase; margin and abstention need NOT be monotone at every intermediate dose |

For B, z/u and any source scale are generated/frozen without target outcomes. Compare at least:
1. ARI(frozen RPP assignments, external labels), where available;
2. ARI(frozen RPP assignments, target value-space clustering);
3. ARI(value clustering before perturbation, value clustering after perturbation).

Under unchanged RPP assignments and unchanged reference labels, quantity 1 cannot fall except through a changed reporting subset. Quantity 2 or 3 may fall. Failure to distinguish these ARIs would misrepresent the mechanism. A rank-space baseline should also be invariant in B; this is not RPP-specific superiority.

For C, simple permutation of paired participants while aligning tumour and normal positions can accidentally preserve the dominant tumour/normal contrast. Randomisation must break assignment–validation association without using target labels: include independently randomised within-pair orientations, or another mathematically specified zero-association construction. Record its assumptions and verify it with unit tests before real outcomes.

For D, removing contradictory rules can increase a remaining score or margin. Assert monotonic coverage, not monotonic abstention at every dose. Full core loss must cause abstention. No surrogate relation search.

## Matched nulls and false reassurance

Use conditional block randomisation of the fixed validation block to preserve relation count, feature reuse/degree, marginal signature distribution and assignment/profile prevalence while destroying its association with fixed assignment. Preserve participant structure under the declared randomisation. Report separately:
- conditional endpoint false reassurance, given the source artifact;
- full-pipeline false structure under source-null simulations, which requires refitting the complete learning/selection pipeline.

Do not call the former unconditional error control for the entire learning procedure. Likewise, a continuous distribution is not automatically a null for narrow operational retention.

Proposed 999 fixed-seed conditional null replicates per intact target, with a binomial uncertainty interval explicitly interpreted as Monte Carlo uncertainty conditional on the observed matrices/randomisation model. Calibration and evaluation null seeds must differ. Add a nuisance-preserving control where feasible; arbitrary gene permutation alone is too easy. The final protocol must define the controlled non-transport property, error denominator and prespecified error bound.

## Explanatory baselines

Minimum: frozen value centroid, frozen rank centroid, post-hoc RPP and RR_DIRECT, with source-only transformations/rejection calibration. Same discovery feature budget and independent validation information for representation-only comparisons. End-to-end comparison is secondary and may evaluate different discovered groups; report each C_source alongside retention. No target refit, target-selected complexity or favorable-coverage matching.

## Implementation sequence

1. Obtain inputs and verify metadata/mappings without target result inspection.
2. Finish the protocol/config/manifests; commit the genuine freeze and reclassification.
3. Add versioned independent-signature artifact, source splits, K/source gate and nested diagnostics.
4. Test participant isolation, entire discovery/signature separation and identical source/target artifact use.
5. Fit source, select signature, estimate reference contrast and freeze full artifact.
6. Execute A, record assignments; evaluate contrasts and then secondary labels.
7. Execute B/C/D with fixed seeds/doses; matched nulls; all baselines.
8. Write exact CSVs, schema/manifest checks and code-generated figures.
9. Publish all findings, including insufficient support; no V1 edits after results.
10. Only then write empirical grant replacement text and BEFORE/AFTER evidence assessment.

Outputs follow the user's requested preliminary/results, figures and docs structure. Before real-data execution, all such empirical result files remain absent. A preflight test log is not a pilot results table.
