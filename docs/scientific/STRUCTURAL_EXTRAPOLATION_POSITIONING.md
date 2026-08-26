# 25. EXTRAPOLATION POSITIONING AND RELATIONAL TRANSPORTABILITY RADIUS

## Why the extrapolation analogy is useful

The project does have an extrapolative component, but it is important to use the word carefully.

Classical **extrapolation** usually means predicting a numerical response outside the range/support represented in the training data. The proposed RPP framework instead asks whether a **scientific group definition** learned in one cohort can be executed beyond the cohort in which it was discovered.

A precise reviewer-facing formulation is therefore:

> **RPP transfer has a structural-extrapolation interpretation: we carry the semantics of a frozen group definition beyond its discovery cohort, rather than extrapolating feature values or a numerical response.**

For formal sections of the proposal, the preferred terms remain **transportability**, **cross-cohort transfer**, **domain generalization** and **robustness under distribution shift**. “Structural extrapolation” is useful as an intuitive explanation, not as the main technical label.

## Relational Transportability Radius (RTR)

The extrapolation analogy suggests a stronger formal object than a generic robustness score: a **source-derived radius within which a frozen RPP assignment is certified to remain unchanged for a declared perturbation class**.

For a frozen profile relation `r=(a,b,d)` and sample `x`, define the oriented margin

`m_r(x)=d(x_a-x_b)`.

Under a symmetric feature-wise additive perturbation `|delta_i|<=epsilon`, the relation is guaranteed satisfied when `m_r(x)>2epsilon` and guaranteed unsatisfied when `m_r(x)<-2epsilon`.

For profile `P_k`, this yields robust lower and upper score bounds:

`L_k(x,epsilon) = weighted fraction of relations with m_r(x)>2epsilon`

`U_k(x,epsilon) = weighted fraction of relations with m_r(x)>-2epsilon`.

For a sample currently assigned to profile `k`, the assignment is certified at perturbation `epsilon` if:

1. executable coverage remains above the frozen coverage threshold;
2. `L_k(x,epsilon)` remains above the frozen score threshold;
3. `L_k(x,epsilon) - max_{l != k} U_l(x,epsilon)` remains above the frozen assignment-margin threshold.

The **pointwise transportability radius** is

`rho(x,k)=sup{epsilon: all three conditions hold}`.

This is a direct robustness/extrapolation radius for the **full frozen assignment rule**, not merely for one gene pair.

## Source-population certificate without reusing fitting data

A statistical certificate must not be calibrated on the same observations used to select the profile. Therefore the confirmatory source workflow will separate:

`SOURCE-FIT -> learn RPP and thresholds`

from

`SOURCE-CALIBRATION -> assign with frozen RPP -> compute rho(x,k)`.

Within each profile, the calibration radii are then converted into a one-sided **distribution-free tolerance bound based on order statistics**. For desired population coverage `q` and confidence `1-alpha`, choose the largest order rank `r` satisfying

`P[Binomial(n,1-q) >= r] >= 1-alpha`.

The observed `r`-th smallest radius is a lower tolerance bound `R_k(q,1-alpha)`. Under the calibration assumptions, with confidence at least `1-alpha`, at least a fraction `q` of the source-profile population has a pointwise certified perturbation radius of at least `R_k`.

If the calibration sample is too small, the method returns **NOT_CERTIFIABLE**. This is scientifically important: a modest discovery cohort is not made reliable by rhetoric. Limited `n` directly limits the strength of the certificate.

For example, even the most conservative first-order 90%-coverage / 95%-confidence one-sided tolerance bound requires at least **29 independent calibration observations** within a profile. Higher coverage/confidence requires more.

## What RTR does and does not mean

RTR provides a mathematically explicit answer to:

> **How far can this frozen relational group definition be perturbed, under a declared shift model, before we lose a source-supported guarantee of preserving its assignment?**

It does **not** guarantee real cross-cohort success. A target cohort may differ through structural/mixture shift, feature-specific distortions outside the declared budget, unmodelled biology or missing platform features. That is why the project retains prospective external validation and the applicability map.

The proposed sequence is therefore stronger than empirical domain generalization alone:

**LEARN -> TRUST -> CERTIFY -> FREEZE -> TRANSFER -> TEST THE CERTIFICATE -> MAP FAILURE MODES.**

The four-WP grant structure does not change; CERTIFY is part of TRUST.

## Novelty boundary relative to certified robustness literature

Certified robustness under distribution shift already exists, particularly for supervised predictive models and neural networks. The project should not claim that “certification under shift” is itself new. The methodological contribution is the adaptation and development of certification for a different scientific object:

- an **unsupervised** patient-group definition;
- represented as sparse executable within-sample relations;
- frozen before target evaluation;
- calibrated without phenotype labels;
- capable of abstention;
- coupled to a real prospective multi-cohort transfer test.

Relevant conceptual neighbours include certified model accuracy under bounded distribution shifts and formal robustness verification for real-world shifts. The proposal should cite them to strengthen, rather than weaken, the novelty argument.


## Literature note

The interpolation/extrapolation distinction is not invented here: Rosenfeld, Ravikumar and Risteski (AISTATS 2022) formalised domain generalization as interpolation among versus extrapolation beyond source distributions. That work concerns predictive risk across domains; RPP uses the analogy for a different object - a frozen unsupervised group definition. This supports the intuition while keeping the novelty claim narrow.
