# From Cohort-Specific Omics Patterns to Transferable Relational Patient Profiles

## 1. SCIENTIFIC GOAL OF THE PROJECT

### 1.1 Scientific problem: why cohort-specific discovery is not enough

High-dimensional omics studies frequently discover patient groups in one cohort and then struggle to reuse those groups in another. The difficulty is not only that clustering algorithms disagree. A deeper problem is that the result of clustering is usually tied to a dataset-specific geometry: distances between normalised values, centroids, medoids, latent embeddings or graph neighbourhoods. When the cohort, laboratory, measurement platform, feature coverage or population composition changes, investigators commonly repeat the analysis from the beginning. A second internally stable partition may then be obtained, but there is no guarantee that it represents the same scientific groups as in the original study.

This is particularly limiting in biomedical omics. Cohorts are often modest, independent replication is expensive, and many clinically interesting subgroups are rare. A useful discovery should therefore be more than a partition of the source matrix. It should be representable as a compact scientific object that can be frozen, published, independently tested and executed for a new sample without reconstructing the target cohort.

The project asks:

**Can patient-group structure discovered in one omics cohort be encoded as a compact relational patient profile that remains executable, interpretable and scientifically testable under cohort and platform shift without target-guided retraining - and can we determine when such transfer should be rejected?**

The proposed answer is a **Relational Patient Profile (RPP)**: a sparse executable definition of a discovered group based on within-sample feature-order relations such as `gene A > gene B`. The name denotes a group-level profile against which an individual patient sample is scored. An RPP is not merely a post-hoc explanation of a conventional cluster. The profile and group membership are learned jointly, so the same object defines the group, explains assignment and can be frozen for future single-sample execution.

This proposal deliberately separates three stages that are often mixed: **unsupervised discovery**, **frozen transfer**, and **external phenotype evaluation**. Clinical diagnosis, subtype or outcome labels will not be used to create RPPs. They are revealed only after the source profile and target assignments have been frozen and are used as external evidence of biological usefulness. Supervised phenotype classification is therefore not a second co-equal aim of the project.

### 1.2 Why within-sample relations may transfer - and why they may fail

Suppose a patient sample is represented by measurements x=(x1,...,xp). A simple relational feature is `r_ab(x)=I(x_a > x_b)`. If all measurements in a sample undergo the same strictly increasing transformation, the ordering is preserved. More generally, many changes in global scale or monotone sample-wise transformation alter absolute values while leaving a substantial part of the within-sample order unchanged. This motivates relational representations as candidates for cross-study reuse.

However, **invariance is not assumed universally**. Gene-specific offsets, platform-specific probe effects, feature loss, nonlinear distortions, low signal-to-noise ratio, biological mixture changes and disease heterogeneity can reverse or remove individual relations. The project therefore treats transportability as a falsifiable property, not as a built-in promise. The goal is to identify the regimes in which relational profiles are useful and the regimes in which they should abstain.

Two explicit negative outcomes are part of the method:

- `NO_STABLE_STRUCTURE`: the source cohort does not support a sufficiently stable and identifiable group structure to justify a reusable RPP;
- `UNASSIGNED`: a source RPP exists, but a target sample has insufficient executable relations or an insufficient assignment margin for a justified transfer.

### 1.3 Main objective and hypotheses

The **main objective** is to develop the theoretical and computational foundations for converting cohort-specific omics discoveries into reusable relational patient profiles that can be transferred to independent cohorts and individual samples without target-guided retraining.

These profiles will be learned directly as sparse executable sets of within-sample feature-order relations, rather than derived post hoc from conventional value-space clusters.

The **central hypothesis** is conditional:

**Stable within-sample feature-order relations can encode patient-group structure in a form that is less sensitive to selected classes of cohort and platform shift than conventional value-space descriptions. When a group is supported by stable relative orderings, a sparse executable profile can retain meaningful source-to-target agreement and single-sample executability without target-guided retraining. When this assumption is violated, source-only reliability criteria and target-side coverage or assignment-margin criteria should lead to abstention rather than forced transfer.**

Three operational hypotheses make the programme falsifiable:

**H1 - LEARN.** When latent group structure is generated or supported by an identifiable sparse partial order, direct relational-profile induction can recover accurate assignments and a compact executable definition under matched information budgets.

**H2 - TRUST AND CERTIFY.** Source-only resampling stability, relation recurrence, assignment consistency, null-calibrated structure strength and identifiability diagnostics can distinguish supported structure from forced fixed-K partitions and can justify `NO_STABLE_STRUCTURE`. For explicitly declared perturbation classes, a frozen RPP evaluated on an independent SOURCE-CALIBRATION subset can additionally yield a distribution-free lower tolerance bound on an assignment-preserving Relational Transportability Radius (RTR). Insufficient calibration evidence yields `NOT_CERTIFIABLE` rather than a weakened post-hoc guarantee.

**H3 - TRANSFER.** A complete source-fitted RPP artifact can retain useful target assignments without target-driven tuning under specific classes of cohort/platform shift; degradation in executable coverage, relation-flip rate and assignment margin will identify conditions in which transfer should be rejected.

The project follows a simple sequence: **LEARN -> TRUST -> TRANSFER -> MAP THE LIMITS**.

### 1.4 Scientific objectives

**O1. Formalise RPPs and direct induction.** Define sparse partial-order prototypes, equivalence/redundancy, identifiability and deterministic direct-induction objectives.

**O2. Establish trust, certification and abstention.** Develop source-only criteria for clusterability, profile stability and reliability; derive conditional Relational Transportability Certificates (RTC) for declared perturbation classes and frozen platform mappings; and retain explicit `NO_STABLE_STRUCTURE` and `UNASSIGNED` decisions.

**O3. Test frozen transfer prospectively.** Apply one unchanged source artifact to independent cohorts under a strict source/target firewall, including a predeclared two-target lung module.

**O4. Map applicability boundaries.** Determine which controlled distribution shifts preserve, degrade or destroy transferability and release a reusable profile specification and open benchmark.

The project is basic methodological research. Omics is the primary empirical testbed, but the mathematical object is more general: it applies wherever scientifically meaningful within-object order relations exist.

## 2. SIGNIFICANCE AND STATE OF THE ART

### 2.1 Reproducibility versus transportability

Cluster stability is an established concern in unsupervised learning [2-4]. Resampling-based stability and prediction strength test whether a partition persists under perturbation of a source dataset [2,3]. Such procedures are valuable but answer a different question from cross-cohort reuse. A group can be stable under resampling and still fail after platform change, feature loss or population shift. Conversely, a useful transferable group definition must be executable without redefining the group on the target cohort.

This distinction parallels a broader reproducibility problem in computational biology: many analyses can be repeated computationally, but scientific conclusions are difficult to reproduce across independent datasets. The proposed project focuses on **transportability of the learned scientific object**, not merely rerunning the same software.

### 2.2 Relative orderings in biomedical modelling

Within-sample relative expression orderings have a substantial biomedical history [6-10,16]. Top-Scoring Pair and related Relative Expression Analysis methods showed that a small number of comparisons can yield transparent supervised decisions [6-8]. REO-based methods later supported reference-based single-sample inference and robust differential analysis. System-level approaches use rank conservation or network-level ordering to describe biological states.

These methods establish that relative orderings can carry useful information and can reduce dependence on some forms of absolute scale. They do not solve the problem proposed here. Most foundational pair-based methods are supervised: outcome labels guide relation selection. Other REO methods compare a sample against reference distributions or transform the data before a conventional downstream analysis. In contrast, the proposed RPP must be **discovered without phenotype labels**, must define the group itself, and must remain executable after the target domain is sealed off from model selection.

The PI's 2026 taxonomic review of rank-based relational methods provides a direct conceptual foundation. It also identifies limitations relevant to the present project: combinatorial relation spaces, lack of a general statistical foundation, heterogeneous evaluation protocols and the need for reproducible benchmarks.

### 2.3 Interpretable clustering and the novelty boundary

Interpretable clustering is an established field [11-15]. Conceptual clustering, order-preserving biclustering, rule-based clustering and optimisation-based interpretable partitions all demonstrate that clusters can have human-readable descriptions. The proposal therefore does **not** claim the first rule-described cluster, the first interpretable clustering algorithm, or the first use of pairwise relations.

The specific methodological gap is narrower and testable. The proposed object combines four properties:

1. the cluster definition is a **sparse within-sample partial order**, not merely an axis-aligned value rule or a high-dimensional transformed representation;
2. membership and profile are **learned jointly**, making the rule set the operational cluster definition rather than a post-hoc explanation;
3. the resulting artifact is **frozen and executable on one unseen sample** without target reclustering or target-guided parameter adaptation;
4. the framework has **explicit rejection modes** and an empirical applicability map under distribution shift.

This combination is the scientific novelty. The project will test it against matched transform-then-cluster and post-hoc alternatives rather than relying on novelty by terminology.

A useful way to delimit the contribution is by the combination of properties. Supervised TSP/Relative Expression methods provide executable pair rules but do not discover unsupervised patient groups. Transform-then-cluster REO methods discover groups but the transformed representation is not itself a frozen group definition. Interpretable/rule-based clustering can provide readable cluster descriptions, but frozen single-sample execution, explicit source/target separation and abstention are not generally the defining object. Cluster-then-classify pipelines permit future assignment, but the classifier is a post-hoc supervised surrogate rather than the discovered unsupervised group definition. **RPP combines unsupervised group discovery, an executable within-sample partial-order definition, frozen single-sample execution, explicit abstention and a source-derived conditional transportability certificate.** No claim is made that any one of these ingredients is individually unprecedented.


### 2.4 Relation to domain generalisation and transfer learning

Distribution shift is extensively studied in supervised machine learning through domain adaptation, transfer learning and domain generalisation. Many domain-adaptation methods explicitly use target-domain observations to adapt a representation or model. That is useful in prediction, but it is not the scientific question here.

The project asks a stricter unsupervised question: **can a patient-group definition be learned in the source, sealed, and then executed in the target without using target outcomes or target-distribution information to redefine the group?** This resembles domain generalisation more than conventional target adaptation, but the output is an unsupervised group definition rather than a supervised predictor.

The project will therefore not market the method as "transfer learning". Instead it studies **frozen cross-cohort transfer under distribution shift**.

Domain-generalisation theory itself has been framed in terms of interpolation among observed source distributions versus extrapolation beyond them [28]. Conceptually, RPP reuse has the same extrapolative direction, but the object carried beyond the discovery cohort is a **frozen patient-group definition**, not a numerical response. I use the phrase *structural extrapolation* only as an intuition; the formal terminology throughout the proposal remains transportability, domain generalization and robustness under distribution shift.

Certified robustness under distribution shift is itself an established research direction, including guarantees for supervised predictive-model accuracy and neural-network behaviour under specified shifts [25,26]. The proposed novelty is not the generic idea of certification; it is the development of a source-calibrated certificate for an **unsupervised executable relational group definition** and its prospective test against frozen multi-cohort transfer.

### 2.5 Why reuse matters for small and heterogeneous cohorts

A small discovery cohort cannot be made reliable merely by choosing a different representation. The project does not claim otherwise. Its contribution is that a result from a modest cohort can be expressed in a form that is **explicitly falsifiable in later cohorts**. A subsequent study can execute the same profile and report success, failure or abstention rather than rediscovering a new set of clusters and comparing them only retrospectively.

This changes the unit of evidence from "a partition observed once" to "a versioned group definition that can accumulate external evidence". In the longer term, such profiles could be stored in a registry together with source provenance, required features, thresholds, applicable platforms and external-validation results. Building a clinical portal is outside the project; establishing the scientific specification required for such reuse is inside the project.

## 3. CONCEPT AND WORK PLAN

### WP1. LEARN - Learning Relational Patient Profiles (M1-M12)

**Aim.** Formalise RPPs and develop deterministic algorithms that learn group membership and sparse relational definitions jointly.

An RPP for group k will be represented as:

`P_k = {(r_kj, w_kj)} for j=1,...,m_k`

where each relation is an executable comparison such as `x_a > x_b`, and `w_kj` is an optional non-negative weight. A simple score for sample x is:

`s(x,P_k) = sum_j w_kj I[r_kj(x) is executable and satisfied] / sum_j w_kj I[r_kj(x) is executable]`.

The denominator makes missing target features explicit rather than silently converting missing relations into failures or ties.

The direct objective will balance: within-group relation agreement; between-group contrast; profile sparsity; redundancy; concentration on too few genes; and executable coverage. Deterministic assignment/profile updates and canonical tie-breaking will ensure that a fixed input and configuration produce byte-identical artifacts.

The formal work will distinguish a **profile** from a **specific edge list**. In some data-generating regimes, many relation sets are equivalent descriptions of the same latent group. Exact designated-edge recovery is then not identifiable. WP1 will define equivalence classes and partial-order reductions so that the method is not rewarded or penalised for arbitrary representation of equivalent structure.

**Matched baselines.**
- value-space PAM or another deterministic medoid baseline;
- relation-space PAM using pair-Hamming or related relational distance;
- rank-space distance baseline;
- `RR_POSTHOC`: relational clustering followed by profile extraction;
- direct RPP induction (`RR_DIRECT` family).

All primary comparisons will use matched source feature budgets and predeclared relation budgets.

**Milestone M12.** Formal RPP definition, deterministic direct-induction implementation, identifiability/equivalence tests and controlled synthetic benchmark.

### WP2. TRUST - Knowing When to Trust a Profile (M7-M24)

**Aim.** Determine whether a source cohort supports a reusable profile before any target outcome is considered.

WP2 adds a formal **Relational Transportability Certificate (RTC)** and **Relational Transportability Radius (RTR)** to the source-side trust layer. The certificate is calibrated **after** the profile has been learned: a SOURCE-FIT subset learns and freezes the RPP, K, relations, weights and assignment thresholds; an independent SOURCE-CALIBRATION subset is then executed by this frozen artifact. For each assigned calibration sample, robust lower/upper profile-score bounds under a declared bounded perturbation yield a pointwise radius `rho(x,k)` within which the full score and assignment-margin conditions remain guaranteed. An exact one-sided distribution-free tolerance bound on these independent radii then yields `RTR_k(q,1-alpha)`: a lower bound such that, under the stated calibration assumptions, with confidence `1-alpha` at least proportion `q` of the source-profile population has a certified radius at least `RTR_k`. If the calibration sample is too small, the result is explicitly `NOT_CERTIFIABLE` rather than a weaker post-hoc claim.

A separate platform-executability certificate will be calculated from frozen annotation/mapping metadata only: the weighted fraction of profile relations whose two features are present and mappable on a target platform. It is an upper bound on executable profile coverage and uses no target expression values or labels. The RTC/RTR is explicitly conditional. It does not certify arbitrary structural shift, emergence of new latent groups or unknown gene-specific distortions. These cases remain empirical applicability questions.


The project rejects the assumption that fixed-K clustering always corresponds to meaningful structure. Trust will be evaluated on source data only using complementary components:

- perturbation and resampling stability of assignments;
- bootstrap recurrence of profile relations or partial-order edges;
- minimum cluster size and non-degeneracy;
- score-margin separation between competing profiles;
- consistency across admissible initialisations/tie resolutions;
- matched NULL processes;
- identifiability diagnostics and structural equivalence.

`NO_STABLE_STRUCTURE` is an explicit scientific output. It will be triggered by a frozen rule combining null-calibrated structure strength and non-degeneracy rather than by subjective visual inspection.

Unknown K will be treated as a separate source-side problem. Candidate K values will be predeclared. Selection will use only source diagnostics; external phenotype agreement cannot select K.

**Milestone M24.** Source-only trust/abstention protocol, independently calibrated RTC/RTR formalism, profile-stability measures, unknown-K extension and neutral matched benchmark.

### WP3. TRANSFER - Transferring Profiles across Cohorts (M18-M38)

**Aim.** Test whether frozen RPPs remain executable and scientifically meaningful in independent cohorts.

The source/target boundary is central. Before a target outcome is opened, the project will freeze:
- feature namespace and mapping rules;
- source preprocessing;
- candidate/selected relations;
- K;
- profile relations and weights;
- score, coverage and margin thresholds;
- software/configuration hashes.

Target data may only **execute** this artifact. Target labels, target phenotype frequencies and target clustering results cannot alter it.

#### Primary prospective lung module

The primary confirmatory design is frozen as:

**Source:** GSE19804 (previously used in the pilot as source evidence)
**Target 1:** GSE27262, GPL570, untouched outcome target
**Target 2:** GSE32863, GPL6884, untouched outcome target; metadata previously audited but no RR_DIRECT outcome evaluation/tuning

One source artifact must be applied unchanged to both targets.

The primary gate is intentionally stringent:

- executable coverage >= 0.70 on Target 1 **and** Target 2;
- forced all-sample ARI >= 0.50 on Target 1 **and** Target 2.

A failed target cannot be replaced after evaluation unseal. Assigned-sample ARI/NMI, rejection rate, score margins and relation coverage are secondary outcomes. The forced all-sample measure prevents selective abstention from artificially inflating primary performance.

The common gene universe will be defined from platform annotation metadata before source fitting. Target expression values and labels cannot determine the universe. Identifier namespace, annotation release, multi-probe aggregation and mapping hashes will be frozen before source fitting. Low coverage produces `UNASSIGNED` or gate failure, not target-guided relation substitution.

#### Secondary colorectal module

The module is frozen at accession level:

**Source:** GSE39582
**Targets:** GSE14333 and GSE33113.

Its exact subtype/phenotype endpoint will be finalised after a sample-level metadata audit but **before target outcome labels are used for model evaluation**. The colorectal module cannot rescue failure of the primary lung gate.

**Milestone M38.** Completed prelabel-sealed source artifacts, target assignments and prospective two-target lung decision plus independent colorectal transfer evidence.

### WP4. MAP THE LIMITS - Mapping the Boundaries of Transferability (M31-M48)

**Aim.** Determine which shifts preserve, weaken or destroy RPP reuse.

Controlled experiments will separate shift mechanisms rather than treating "batch effect" as one undifferentiated nuisance:

**Class A - order-preserving or approximately preserving shifts.**
Shared strictly monotone sample-wise transformations preserve pair order exactly. Bounded feature-specific additive perturbations admit margin-based sufficient conditions and will provide the primary RTC calibration regime. Global scale changes are included as favourable controls.

**Class B - relation-flipping shifts.**
Gene-specific offsets/scales, local nonlinearities or feature-specific calibration changes. These directly test the limits of pairwise invariance.

**Class C - missing-feature/platform-overlap shifts.**
Feature loss, annotation mismatch and reduced platform overlap. These test executability and `UNASSIGNED`.

**Class D - structural/mixture shifts.**
Changes in cluster prevalence, emergence/disappearance of a group, nuisance structure and continuous rather than discrete structure. These test whether a source profile should be transferred at all.

The final applicability map will relate source trust metrics and target-side executable coverage/margins to transfer outcomes. Negative regimes are a planned scientific result.

**Milestone M48.** Applicability map, registry-ready machine-readable profile schema, open implementation, benchmark and final synthesis.

## 4. PRELIMINARY RESULTS AND FEASIBILITY

### 4.1 Frozen Representation Audit pilot

A deterministic source-only pilot compared VALUE, RELATIONAL and HYBRID representations across **630 controlled source-target pairs**. It identified the generating representation family in **93.3%** of signal replicates, with median target ARI regret **0.000**, NULL false-structure rate **6.7%**, no HYBRID selections in pure VALUE/RELATIONAL regimes, and source-audit/target-performance Spearman association **0.854**. This supported the feasibility of source-observable reliability diagnostics.

The pilot deliberately retained a negative external result. In the original bidirectional GSE10072/GSE19804 transfer, one direction matched the retrospective oracle, whereas the reverse exceeded the frozen regret allowance. **Gate C remained STOP.** This result was not repaired by changing thresholds. It showed that within-cohort adequacy is not sufficient for transportability.

### 4.2 Direct RPP pilot (RR_DIRECT)

A subsequent direct-prototype pilot tested whether the cluster definition itself could be learned as sparse relations. The implementation used deterministic candidate-pair generation, source-only feature screening, within-support/between-contrast rule selection, prototype-score reassignment, explicit missingness and frozen target scoring.

The prospective pilot achieved **4/5 criteria** and therefore remained formally **STOP**. Key results included:
- synthetic RELATIONAL regimes: median RR_DIRECT ARI **1.000**;
- best frozen transfer: ARI **0.960** with **92.5%** coverage;
- opposite frozen direction: ARI **0.771** with **96.7%** coverage;
- NULL false-structure criterion passed;
- real-data behaviour was heterogeneous rather than universally favourable.

The failed criterion was exact recovery of a designated six-pair ground truth. Inspection showed a non-identifiable generator: shifting two blocks of six genes created many equally discriminatory cross-block relations, so exact recovery of one arbitrary designated pair set was not mathematically justified.

### 4.3 Identifiability diagnostic

A preregistered-style diagnostic used separated order blocks with the **same hyperparameters**, not a tuned rescue. Across **120 replicates**, median source ARI, exact-pair recovery, target ARI and coverage were all **1.0**, and all replicates achieved exact recovery and ARI >=0.75. This does not convert the v2 STOP into a PASS; instead it demonstrates that identifiability is a genuine scientific variable requiring formal treatment in WP1/WP2.

### 4.4 Eleven real omics datasets

RR_DIRECT was also run descriptively on eleven real expression datasets. Some cohorts showed strong agreement, including GSE10072 (ARI 0.926), GSE19804 (0.839), Colon (0.446) and DLBCL (0.329), but the median RR_DIRECT ARI across all eleven datasets was only **0.033**. Relation/PAM median ARI was **0.022**. These results explicitly argue against a universal-superiority narrative and motivate the applicability-map objective.

### 4.5 Feasibility and reproducibility infrastructure

The current public repository contains deterministic Python code, frozen protocols, real-data adapters, prospective external-cohort freeze documentation and a test suite with **155 passing tests**. The infrastructure separates fit-time and evaluation labels, records source artifacts, and prevents target labels from modifying source preprocessing or model parameters. The SONATA BIS work will extend this infrastructure rather than rebuild already completed pilot components.

## 5. METHODOLOGY

### 5.1 Data model and leakage boundary

Each dataset will be represented by a measurement bundle and a separate evaluation-label layer. Fit and audit modules will not import target outcome labels. A primary run will proceed through:

1. metadata/feature audit;
2. source preprocessing and relation-candidate definition;
3. source-only profile learning and trust evaluation;
4. serialisation and hashing of the source artifact;
5. target execution without target-guided changes;
6. sealing of assignments, scores, coverage and hashes;
7. only then, outcome-label unsealing and evaluation.

A regression test will modify target labels and verify byte-identical source artifacts. A second test will compare individual versus batch execution to ensure that the assignment of one target sample does not depend on other target samples.

### 5.2 Formal profile model

Let candidate relation `r_j(x)` take values satisfied/not satisfied/non-executable. For group k, RPP `P_k` contains a sparse subset E_k of relations and optional weights. The executable-score form is:

`score_k(x) = sum_{j in E_k} w_kj I[r_j(x)=1] / sum_{j in E_k} w_kj I[r_j(x) is observed]`.

A target assignment requires:
- minimum executable coverage `c(x,P_k) >= c_min`;
- score above a source-calibrated threshold;
- margin over the second-best profile above `delta_min`.

Otherwise the result is `UNASSIGNED`.

Direct induction will optimise a transparent multi-component objective of the form:

`L = within-profile disagreement + lambda1*profile_size + lambda2*redundancy + lambda3*instability + penalties for degenerate groups/low coverage`.

The exact component scales and finite hyperparameter grids will be declared before confirmatory evaluation.

### 5.3 Relational Transportability Certificate and Radius

The RTC makes the theoretical prediction explicit before target outcomes are known. For relation `r=(a,b,d)` define signed margin `m_r(x)=d(x_a-x_b)`. Under a symmetric feature-wise additive perturbation `|delta_i|<=epsilon`, the relation is guaranteed satisfied when `m_r(x)>2epsilon` and guaranteed unsatisfied when `m_r(x)<-2epsilon`.

For frozen profile `P_k`, let `L_k(x,epsilon)` be the weighted fraction of executable relations guaranteed satisfied and `U_k(x,epsilon)` the weighted fraction not guaranteed unsatisfied. If `s_min`, `Delta_min` and `c_min` are the frozen score, margin and executable-coverage thresholds, the assignment of sample `x` to `k` is certified at `epsilon` when `coverage>=c_min`, `L_k>=s_min`, and `L_k-max_{l!=k}U_l>=Delta_min`. The **pointwise transportability radius** `rho(x,k)` is the largest perturbation budget satisfying all three conditions. Because the bounds change only at finite relation-margin breakpoints, the radius is deterministic and exactly computable for the declared perturbation class.

To avoid post-selection optimism, the same source observations will **not** both learn and certify the profile. The confirmatory source workflow will use a SOURCE-FIT subset to learn and freeze the RPP, K, relations, weights and thresholds, followed by an independent SOURCE-CALIBRATION subset on which the frozen artifact is executed. Calibration samples assigned to profile `k` provide pointwise radii `rho_1,...,rho_n`.

A one-sided **distribution-free tolerance bound** converts those radii into a population-level Relational Transportability Radius `RTR_k(q,1-alpha)`. We choose the largest order rank `r` satisfying `P[Binomial(n,1-q)>=r] >= 1-alpha`; the `r`-th smallest calibration radius is then a lower bound such that, under the fixed-model/IID calibration assumptions, with confidence at least `1-alpha`, at least proportion `q` of the source-profile population has certified radius at least `RTR_k`. If the calibration subset is too small, the formal outcome is `NOT_CERTIFIABLE`. For example, even a first-order 90%-coverage / 95%-confidence bound requires at least 29 independent calibration observations within a profile. Thus modest source cohorts explicitly limit certificate strength rather than being portrayed as reliable by assumption.

Under platform shift, a separate metadata-only certificate records the weighted fraction of profile relations executable after a mapping frozen before source fitting. The formal RTC does not certify arbitrary structural/mixture shift or unknown feature-specific distortions outside the declared perturbation class. Controlled experiments test certificate calibration; real targets test whether stronger source certificates predict better frozen transfer. Failure despite a strong certificate is not repaired: it identifies a real shift mechanism outside the certified model or a limitation of the source assumptions.

### 5.4 Candidate relations and computational control

The full pair space is O(p^2) and will not be materialised blindly at genome scale. Source-only screening will first define an admissible feature set using frozen variance/MAD and missingness criteria. Candidate pairs will then be bounded by a declared maximum or by source-only stability/contrast screening. No target values will nominate substitute features or relations.

Inverse relations are redundant (`A>B` versus `B<A`) and will be canonicalised. Excessive reuse of a single gene will be penalised or capped to avoid profiles that appear sparse in edge count but depend on one unstable anchor.

Primary computation is deterministic, CPU-first Python. Parallel execution may be used for independent resamples/cohorts but not in a manner that changes numerical tie-breaking.

### 5.5 Identifiability and structural equivalence

Exact edge recovery is meaningful only when the generating structure is identifiable. WP1 will distinguish:
- **assignment identifiability**: whether the patient partition is recoverable;
- **profile identifiability**: whether a unique sparse relation set is supported;
- **equivalence**: whether different relation sets imply the same partial order or nearly identical executable region.

Transitive reduction of directed acyclic partial-order representations will be explored to remove edges implied by other edges. Where several sparse descriptions are statistically indistinguishable, the output will be an equivalence class or stable core rather than an unjustified unique rule list.

### 5.6 Source stability and NO_STABLE_STRUCTURE

Source trust will combine multiple statistics rather than one aggregate score hiding failure:
- bootstrap/subsample assignment agreement;
- relation-edge Jaccard/recurrence;
- minimum group fraction;
- profile-score separation;
- null-calibrated structure strength;
- sensitivity to reasonable preprocessing perturbations;
- identifiability/equivalence diagnostics.

The `NO_STABLE_STRUCTURE` rule will be frozen on controlled data before use in prospective real-cohort transfer. The method is allowed to conclude that the source cohort does not support a defensible discrete partition.

### 5.7 Cross-cohort mapping and platform shift

Cross-platform mapping is a major source of hidden flexibility. For the prospective lung module, the feature universe will be established using platform annotation metadata only, before source fitting. The protocol will freeze:
- identifier namespace;
- annotation release;
- mapping table;
- rule for multiple probes per gene (working default: within-sample median, subject to exact pre-fit freeze);
- hashes of the mapping artifacts.

Target expression values and labels cannot influence feature-universe construction. Relations requiring unavailable genes become non-executable and reduce coverage.

### 5.8 Controlled shift experiments

Synthetic generators will independently manipulate:
- global sample-wise scale and monotone transformations;
- gene-specific additive/multiplicative distortions;
- relation-flip probability;
- missing-feature fraction and structured platform overlap;
- signal strength and dimensionality;
- cluster-number uncertainty;
- unequal prevalence and emergence/disappearance of groups;
- nuisance variables correlated or uncorrelated with group structure;
- continuous versus discrete latent structure.

Factorial or fractional-factorial designs will be used where computationally justified. The objective is not to generate one average performance number but a **mechanistic response surface** showing which shift components drive failure.

### 5.9 Baselines and role of supervised classification

Primary baselines will be unsupervised and matched:
- value-space PAM/medoid clustering;
- rank-space clustering;
- pair-relation clustering;
- post-hoc relational profile extraction;
- selected standard clustering alternatives used for sensitivity analysis.

Supervised phenotype classifiers may be included only as **secondary comparators** to answer a limited question: how much phenotype association is obtained when labels are explicitly optimised, compared with the association retained by unsupervised frozen RPPs? They cannot select the RPP, K, feature universe, threshold or target direction. Thus classification is contextual evidence, not a second aim.

### 5.10 Statistical evaluation

For controlled data, outcomes include assignment ARI/NMI, relation/profile recovery where identifiable, regret to a retrospective oracle, false-structure rate under NULL, profile size and stability, RTC/RTR violation rate, calibration of population-level assignment preservation at predeclared coverage/confidence, NOT_CERTIFIABLE rate due to limited calibration n, and false-reassurance rate under in-class perturbations.

For real transfer, outcomes include:
- executable coverage;
- forced all-sample ARI/NMI;
- assigned-sample ARI/NMI as secondary outcomes;
- `UNASSIGNED` rate;
- score/margin distribution;
- relation flip rate;
- profile-state retention;
- profile length/gene count;
- bootstrap profile Jaccard;
- post-freeze phenotype/pathway associations.

Uncertainty will be estimated by patient-level resampling. Paired samples and technical replicates will remain grouped. All preregistered cohorts and directions will be reported, including negative results. Multiple biological association tests will use appropriate false-discovery control and effect sizes with uncertainty intervals.

### 5.11 Biological interpretation

RPP relations will be mapped to genes and pathways only after the profile is frozen. Enrichment and clinical associations are interpretive/evaluation analyses. They cannot retroactively change profile membership or primary transfer decisions. This prevents biological plausibility from becoming an unrecorded source of model selection.

### 5.12 Reproducibility and open science

Every primary result will be linked to:
- dataset/accession manifest;
- source/target role;
- preprocessing and mapping hashes;
- software commit and dependency lock;
- configuration and seed schedule;
- prelabel/evaluation state;
- immutable compact evidence tables.

Public source datasets will be referenced by accession rather than unnecessarily redistributed. Code, synthetic generators, frozen protocols, machine-readable profile artifacts and derived evidence will be released openly. Major releases will be archived with persistent identifiers.

## 6. COMPOSITION AND QUALIFICATIONS OF THE RESEARCH TEAM

The financed scientific core is intentionally compact: **PI + one doctoral researcher + one postdoctoral researcher**. This architecture satisfies the scientific need for method construction, reliability analysis and independent transfer validation while supporting the SONATA BIS goal of creating a new team.

### Principal Investigator

The PI will lead RPP theory, algorithm design, preregistered experimental architecture, scientific integration and supervision. His research trajectory provides direct but non-duplicative preparation:

- PRELUDIUM 5 (2013/09/N/ST6/04083; 2014-2017) - independent interpretable algorithm design and model-tree induction;
- OPUS 17 (2019/33/B/ST6/02386; 2020-2024) - relative dependencies in integrated omics data;
- 2026 taxonomic review in *Artificial Intelligence Review* - synthesis of rank-based relational methods and their open methodological questions;
- international research experience including the 2026 University of Girona stay and the 2024 Broad Institute/CUNY visit.

The proposed project moves beyond the PI's earlier supervised/exploratory relational work to a distinct unsupervised programme centred on reusable group definitions under distribution shift.

### Doctoral researcher - at least 36 months

The doctoral researcher will develop source-side trust, identifiability and stability methodology, controlled shift experiments, unknown-K extensions and parts of direct profile induction. This forms a coherent dissertation trajectory in interpretable unsupervised learning rather than routine technical support.

### Postdoctoral researcher - 36 months

A postdoctoral researcher selected in an open competition will lead independent cohort curation, platform mapping, frozen-transfer execution, neutral benchmark design and reproducibility audits. Separating substantial external-transfer work from direct method construction reduces confirmation bias and provides the postdoc with an independent scientific line.

Short external consultations may support biological interpretation or statistical review but will not select methods after label unsealing.

## 7. REASONS FOR CREATING A NEW RESEARCH TEAM

The project is not a continuation of the personnel structure of OPUS 17. It creates a new team around **unsupervised relational patient profiling and transportability**, with two new researcher trajectories:

1. reliability, identifiability and abstention of relational profiles;
2. independent multi-cohort transfer and platform-shift benchmarking.

Former collaborators are not automatically included as funded team members. Recruitment will follow current SONATA BIS eligibility rules. The new team structure is scientifically necessary because a credible transfer study benefits from separation between algorithm construction and external-cohort assembly/evaluation.

The programme is also sufficiently broad for team development but sufficiently focused to remain coherent: one scientific object (RPP), one central source/target boundary and one four-stage logic.

## 8. RISK ANALYSIS AND ALTERNATIVE OUTCOMES

**R1. No discrete stable group structure exists.**
Response: `NO_STABLE_STRUCTURE` is a valid primary output. The project will not force K groups because the algorithmic framework expects them.

**R2. Direct induction does not outperform post-hoc profiles.**
Response: the result becomes a negative methodological conclusion defining when joint optimisation is unnecessary. `RR_POSTHOC` remains the supported comparator; no target-driven rescue is allowed.

**R3. Profiles are accurate but not uniquely identifiable.**
Response: report stable cores, partial orders or equivalence classes rather than arbitrary exact edges.

**R4. Cross-platform feature coverage is too low.**
Response: retain the source-defined mapping, produce `UNASSIGNED` or prospective gate failure, and include the event in the applicability map.

**R5. Profiles are stable but phenotype-discordant.**
Response: structural reproducibility and biological/clinical association are reported separately. Labels cannot rescue profile learning.

**R6. Relation invariance is weaker than expected under gene-specific shifts.**
Response: this is a central falsification pathway. The project will identify which shift types break the representation instead of claiming universal robustness.

**R7. Performance depends on one favourable cohort/direction.**
Response: the primary lung protocol requires both independent targets to pass; a failed target cannot be replaced.

**R8. Pair spaces are computationally expensive.**
Response: source-only feature/relation screening, bounded candidate sets, redundancy control and deterministic CPU implementation. No full p-by-p pair matrix is required for primary analyses.

**R9. Source calibration is too small for a formal RTC.**
Response: return `NOT_CERTIFIABLE` at the predeclared coverage/confidence rather than weaken the certificate after seeing the data. Empirical stability/transfer analyses remain separate and cannot be relabelled as formal certification.

## 9. EXPECTED RESULTS AND IMPACT

The central expected result is not simply a new clustering algorithm. It is a **scientific framework for turning a cohort-specific unsupervised discovery into a reusable and falsifiable group definition**.

Expected outputs are:

- a formal theory of sparse relational patient profiles, identifiability, equivalence and conditional transportability certification;
- deterministic direct-induction algorithms or a clear negative result on when direct induction adds value;
- source-only trust and abstention criteria;
- prospective frozen cross-cohort transfer evidence;
- a mechanistic applicability map for cohort/platform/distribution shift;
- an open benchmark and reproducible software implementation;
- a machine-readable, registry-ready RPP specification.

A frozen RPP can in principle be published and executed locally at another centre on one molecular sample. Only the compact profile needs to be shared; the target site's raw molecular matrix need not be centrally pooled for execution. This property is compatible with future distributed or federated architectures, but the present project does not claim federated learning or formal privacy guarantees.

The broader methodological contribution is a shift from asking only **"can we find a stable partition?"** to asking **"can the discovered group be defined in a form that another study can execute, falsify and reuse?"**

## 10. PROJECT LITERATURE

1. Kaufman L, Rousseeuw PJ. *Finding Groups in Data: An Introduction to Cluster Analysis*. Wiley; 1990.
2. Ben-Hur A, Elisseeff A, Guyon I. A stability based method for discovering structure in clustered data. *Pacific Symposium on Biocomputing*. 2002;7:6-17.
3. Tibshirani R, Walther G. Cluster validation by prediction strength. *Journal of Computational and Graphical Statistics*. 2005;14:511-528.
4. von Luxburg U. Clustering stability: an overview. *Foundations and Trends in Machine Learning*. 2010;2:235-274.
5. Hubert L, Arabie P. Comparing partitions. *Journal of Classification*. 1985;2:193-218.
6. Geman D, d'Avignon C, Naiman DQ, Winslow RL. Classifying gene expression profiles from pairwise mRNA comparisons. *Statistical Applications in Genetics and Molecular Biology*. 2004;3:Article 19.
7. Tan AC, Naiman DQ, Xu L, Winslow RL, Geman D. Simple decision rules for classifying human cancers from gene expression profiles. *Bioinformatics*. 2005;21:3896-3904.
8. Eddy JA, Sung J, Geman D, Price ND. Relative expression analysis for molecular cancer diagnosis and prognosis. *Technology in Cancer Research & Treatment*. 2010;9:149-159.
9. Czajkowski M, Jurczuk K, Kretowski M. Rank-based relational methods for interpretable biomedical modeling: a taxonomic review. *Artificial Intelligence Review*. 2026. doi:10.1007/s10462-026-11679-3.
10. Yan J, Zeng Q, Wang X. RankCompV3: a differential expression analysis algorithm based on relative expression orderings and applications in single-cell RNA transcriptomics. *BMC Bioinformatics*. 2024;25:259.
11. Michalski RS, Stepp RE. Learning from observation: conceptual clustering. In: Michalski RS, Carbonell JG, Mitchell TM, eds. *Machine Learning: An Artificial Intelligence Approach*. Tioga; 1983:331-363.
12. Ben-Dor A, Chor B, Karp R, Yakhini Z. Discovering local structure in gene expression data: the order-preserving submatrix problem. *Journal of Computational Biology*. 2003;10:373-384.
13. Bertsimas D, Orfanoudaki A, Wiberg H. Interpretable clustering: an optimization approach. *Machine Learning*. 2021;110:89-138.
14. Carrizosa E, Kurishchenko K, Marin A, Romero Morales D. On clustering and interpreting with rules by means of mathematical optimization. *Computers & Operations Research*. 2023;154:106180.
15. Hu L, Jiang M, Dong J, Liu X, He Z. Interpretable Clustering: A Survey. *ACM Computing Surveys*. 2026;58(8):Article 215.
16. Guan Q, et al. Differential expression analysis for individual cancer samples based on robust within-sample relative gene expression orderings across multiple profiling platforms. *Oncotarget*. 2016;7:68909-68920.
17. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*. 2007;8:118-127.
18. McCall MN, Bolstad BM, Irizarry RA. Frozen robust multiarray analysis (fRMA). *Biostatistics*. 2010;11:242-253.
19. Boulesteix AL, Lauer S, Eugster MJA. A plea for neutral comparison studies in computational sciences. *PLoS ONE*. 2013;8:e61562.
20. Ioannidis JPA, et al. Repeatability of published microarray gene expression analyses. *Nature Genetics*. 2009;41:149-155.
21. Landi MT, et al. Gene expression signature of cigarette smoking and its role in lung adenocarcinoma development and survival. *PLoS ONE*. 2008;3:e1651.
22. Lu TP, et al. Identification of a novel biomarker, SEMA5A, for non-small cell lung carcinoma in nonsmoking women. *Cancer Epidemiology, Biomarkers & Prevention*. 2010;19:2590-2597.
23. Marisa L, et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. *PLoS Medicine*. 2013;10:e1001453.
24. Guinney J, et al. The consensus molecular subtypes of colorectal cancer. *Nature Medicine*. 2015;21:1350-1356.
25. Kumar A, Levine A, Goldstein T, Feizi S. Certifying Model Accuracy under Distribution Shifts. arXiv:2201.12440; 2022.
26. Wu H, Tagomori T, Robey A, et al. Toward Certified Robustness Against Real-World Distribution Shifts. *IEEE Conference on Secure and Trustworthy Machine Learning (SaTML)*. 2023:537-553. doi:10.1109/SaTML54575.2023.00042.
27. Meeker WQ, Hahn GJ, Escobar LA. *Statistical Intervals: A Guide for Practitioners and Researchers*. 2nd ed. Wiley; 2017. doi:10.1002/9781118594841.

28. Rosenfeld E, Ravikumar P, Risteski A. An Online Learning Approach to Interpolation and Extrapolation in Domain Generalization. *Proceedings of AISTATS / PMLR*. 2022;151:2641-2657.
