# 20. RELATIONAL TRANSPORTABILITY CERTIFICATE - FORMAL CORE v2

## Purpose

The **Relational Transportability Certificate (RTC)** is a source-derived, model-relative certificate for a frozen Relational Patient Profile (RPP). It does **not** claim universal invariance or guarantee success in an arbitrary target cohort. It certifies preservation only under an explicitly declared perturbation class and is then tested prospectively against real cross-cohort transfer.

The central formal quantity is the **Relational Transportability Radius (RTR)**: a lower bound on how large a declared perturbation can become while the full frozen RPP assignment remains guaranteed to satisfy its score and margin rules.

## 1. Relation-level robustness

For directed relation `r=(a,b,d)`, `d in {+1,-1}`, define the signed margin

`m_r(x)=d(x_a-x_b)`.

The relation is satisfied when `m_r(x)>0`.

### Exact order invariance

For a common strictly increasing transformation `T`,

`x_a>x_b <=> T(x_a)>T(x_b)`.

This is exact for that transformation class only.

### Bounded feature-specific additive perturbation

Let `x'_i=x_i+delta_i`, with `|delta_i|<=epsilon_i`. Then the signed relation margin can decrease by at most `epsilon_a+epsilon_b`. Under a symmetric bound `epsilon_i<=epsilon`, relation `r` is:

- guaranteed satisfied if `m_r(x)>2epsilon`;
- guaranteed unsatisfied if `m_r(x)<-2epsilon`;
- uncertain otherwise.

## 2. Robust bounds for the complete RPP score

For frozen profile `P_k` with relation weights `w_kr`, define, for one sample and perturbation budget `epsilon`,

`L_k(x,epsilon)=sum_r w_kr I[m_r(x)>2epsilon] / sum_r w_kr`

and

`U_k(x,epsilon)=sum_r w_kr I[m_r(x)>-2epsilon] / sum_r w_kr`.

`L_k` is a worst-case lower bound on the profile score and `U_k` a worst-case upper bound, for the declared additive perturbation class. If features are unavailable, they are excluded exactly as specified by the frozen executable-coverage rule; platform availability is handled separately below.

For a sample assigned to profile `k`, let `s_min` be the frozen minimum score and `Delta_min` the frozen minimum assignment margin. Its assignment is certified at `epsilon` when:

1. executable coverage >= frozen `c_min`;
2. `L_k(x,epsilon) >= s_min`;
3. `L_k(x,epsilon) - max_{l!=k} U_l(x,epsilon) >= Delta_min`.

The **pointwise transportability radius** is

`rho(x,k)=sup{epsilon >=0: conditions 1-3 hold}`.

Because the robust score bounds change only when `epsilon` crosses `|m_r(x)|/2`, `rho(x,k)` can be computed deterministically from finitely many breakpoints rather than from an arbitrary numerical grid.

## 3. Independent source calibration is mandatory

A statistical certificate cannot reuse the same observations that selected the profile and its relations. The confirmatory source workflow is therefore:

**SOURCE-FIT -> learn RPP, K, relations, weights and assignment thresholds -> FREEZE**

**SOURCE-CALIBRATION -> execute the frozen RPP -> compute pointwise radii for assigned samples**.

Only the independent source-calibration subset is used to convert pointwise radii into a population-level RTC. Cross-fitting may be explored for efficiency, but the primary confirmatory certificate uses a clean fit/calibration separation.

## 4. Distribution-free profile-level tolerance certificate

For one frozen profile `k`, suppose `n_k` independent source-calibration samples are assigned to it, yielding radii

`rho_1,...,rho_nk`.

For desired population coverage `q` and confidence `1-alpha`, choose the largest order rank `r` satisfying

`P[Binomial(n_k,1-q) >= r] >= 1-alpha`.

Let `rho_(r)` be the `r`-th smallest observed radius. Then, under the fixed-model/IID calibration assumptions, `rho_(r)` is a one-sided distribution-free lower tolerance bound: with confidence at least `1-alpha`, at least proportion `q` of the source-profile population has pointwise transportability radius at least `rho_(r)`.

Define

`RTR_k(q,1-alpha)=rho_(r)`.

If no order rank satisfies the confidence requirement, the formal outcome is **NOT_CERTIFIABLE**.

The first-order existence condition is

`q^n <= alpha`.

Thus a 90%-coverage / 95%-confidence lower bound requires at least 29 independent calibration observations even before asking for a less conservative order rank. This makes limited source sample size an explicit scientific limitation rather than a hidden weakness.

## 5. Platform executability certificate

Cross-platform transfer adds a distinct failure mode: a mathematically robust relation may be unavailable because one or both features cannot be mapped.

For a platform mapping `A` defined **before source fitting and without target expression or labels**, define

`E_k(A)=sum_r w_kr I[a_r and b_r are mappable under A]/sum_r w_kr`.

This is a metadata-only upper bound on weighted executable profile coverage. Actual per-sample executable coverage may be lower because of missing measurements.

The final RPP artifact therefore separates:

- **assignment robustness certificate / RTR** under a declared value-perturbation class;
- **platform executability certificate** under a frozen mapping;
- **empirical target execution**: actual score, margin, coverage and ASSIGNED/UNASSIGNED status.

## 6. Scientific hypotheses enabled by RTC/RTR

**H2a - calibration validity.** Under controlled in-class perturbations, observed assignment-preservation rates should respect the distribution-free RTC at the predeclared coverage/confidence levels.

**H2b - radius-transfer association.** Larger source-derived RTR values should predict better frozen target execution across controlled shifts and, descriptively, across independent cohorts.

**H2c - informative refusal.** Small source groups may be NOT_CERTIFIABLE; low platform executability may cause UNASSIGNED; neither condition is repaired using target outcomes.

**H2d - failure specificity.** Strong RTC followed by real-cohort failure should identify a shift mechanism outside the certified class, a violation of source-calibration assumptions or a limitation of the relational representation.

## 7. Evaluation

Controlled experiments will report:

- empirical violation rate for pointwise and population-level certificates;
- coverage calibration at predeclared `(q,1-alpha)` levels;
- distribution of pointwise `rho(x,k)` and profile `RTR_k`;
- association between RTR and target ARI/NMI, executable coverage, score margin and rejection;
- false-reassurance rate: a positive certificate followed by failure under an in-class perturbation;
- refusal rate due to insufficient calibration sample size.

Real-cohort evaluation treats RTC/RTR as a **pre-target predictor**, not proof. External targets remain frozen and cannot choose certificate definitions or thresholds.

## 8. Relation to extrapolation and certified robustness

Conceptually, RPP transfer has an extrapolative direction: a frozen group definition is executed beyond its discovery cohort. Domain-generalisation theory has itself distinguished interpolation among observed source distributions from extrapolation beyond them (Rosenfeld, Ravikumar & Risteski, 2022). Here, however, we extrapolate the **semantics of a discovered group**, not a numerical response. The formal proposal therefore uses transportability/domain-generalization terminology and treats *structural extrapolation* only as an intuition.

Certified robustness under distribution shift is an existing research area, including certificates for predictive-model accuracy and neural-network robustness. The novelty claim is therefore not “the first certificate under shift”. The proposed contribution is a certificate tied to an **unsupervised, executable relational group definition**, calibrated without phenotype labels and tested through prospectively frozen multi-cohort transfer.

## 9. Claim boundary

Safe claim:

> The project will develop source-calibrated sufficient conditions and distribution-free tolerance certificates for preserving frozen relational-profile assignments under explicitly declared perturbation classes, and test whether these certificates predict real cross-cohort transfer.

Do not claim:

> The certificate guarantees transfer to any independent cohort or platform.

The certificate is conditional on a fixed profile, independent source calibration and a declared perturbation model. Prospective target studies test how useful those assumptions are in practice.

## 10. Methodological neighbours to cite

- Kumar A, Levine A, Goldstein T, Feizi S. *Certifying Model Accuracy under Distribution Shifts*. 2022. arXiv:2201.12440.
- Wu H, Tagomori T, Robey A, et al. *Toward Certified Robustness Against Real-World Distribution Shifts*. IEEE SaTML 2023:537-553. doi:10.1109/SaTML54575.2023.00042.
- Rosenfeld E, Ravikumar P, Risteski A. *An Online Learning Approach to Interpolation and Extrapolation in Domain Generalization*. AISTATS / PMLR. 2022;151:2641-2657.
- Meeker WQ, Hahn GJ, Escobar LA. *Statistical Intervals: A Guide for Practitioners and Researchers*. 2nd ed. Wiley; 2017. doi:10.1002/9781118594841.
