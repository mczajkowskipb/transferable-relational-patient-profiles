# H2 metadata-only family feasibility audit

**Freeze date:** 2026-09-08  
**Purpose:** pre-outcome capacity audit for SONATA BIS H2 (TRUST)  
**Rows:** 30 candidate source-target blocks  
**Status of every row:** `CAPACITY_SCREENED_NOT_FINAL_ELIGIBLE`

## Result

A conservative desk screen identified **30 non-reused candidate blocks** outside the
development-exposed lung anchor and the principal colorectal H3 family:

- MetaGxOvarian: 9
- MetaGxBreast / curated breast compendium: 10
- MetaGxPancreas: 4
- curatedPCaData: 7

The screen requires nominal source N >= 100 and nominal target N >= 50. No named
dataset is reused in this matching. The candidate pool therefore provides headroom
for the grant's **minimum H2 feasibility floor of 20 final independent blocks**.

## Critical interpretation

This file is **not** evidence that all 30 rows are already statistically independent
or finally eligible. It is a metadata-only capacity screen performed before any H2
RPP outcome is inspected.

At M6, the outcome-blind dependency/eligibility audit must:

1. collapse participant overlap, publication derivatives, subset studies and platform
   splits into dependency components;
2. verify disease/material/layer compatibility and participant-level effective N;
3. verify identifier mapping, feature overlap and source-defined executability rules;
4. freeze one deterministic non-reuse matching;
5. exclude failing blocks without outcome-driven replacement.

Primary H2 proceeds only with >=20 independent blocks and >=8 `Y_ext=1` plus >=8
`Y_ext=0` blocks. Pre-outcome operating-characteristic simulation may raise these
minima, never lower them. Otherwise H2 is reported as inconclusive.

## Why the supplied 21-disease list was not copied verbatim

A different disease label does not by itself prove eligibility, adequate N, compatible
material or absence of hidden dependence. Several proposed pairs are plainly too small
for the stated planning screen. GEO reports only 47 individual tissue samples in
GSE8397 and 18 samples in GSE20141 for the proposed Parkinson block; GSE1919 contains
15 samples and GSE55235 79 individuals for the proposed rheumatoid-arthritis block;
GSE41762 contains 77 islet samples and GSE25724 only 13 for the proposed type-2-diabetes
block. Those rows were therefore not used to manufacture N.

The present census instead uses curated transcriptomic resources with explicit sample
counts and, for MetaGxBreast/MetaGxOvarian/MetaGxPancreas, built-in duplicate-removal
machinery.

## Evidence base

- MetaGxData: Gendoo DMA et al. *Scientific Reports* 2019;9:8770.
  https://doi.org/10.1038/s41598-019-45165-4
- MetaGxOvarian current vignette:
  https://bioconductor.posit.co/packages/3.24/data/experiment/vignettes/MetaGxOvarian/inst/doc/MetaGxOvarian.html
- MetaGxBreast duplicate handling:
  https://www.bioconductor.org/packages/release/data/experiment/manuals/MetaGxBreast/man/MetaGxBreast.pdf
- Curated breast sample-count table:
  https://aacrjournals.org/clincancerres/article/22/2/337/175104/Subtype-Specific-Metagene-Based-Prediction-of
- MetaGxPancreas current vignette:
  https://bioconductor.statistik.tu-dortmund.de/packages/3.23/data/experiment/vignettes/MetaGxPancreas/inst/doc/MetaGxPancreas.html
- curatedPCaData: Laajala TD et al. *Scientific Data* 2023;10:430.
  https://doi.org/10.1038/s41597-023-02335-4
- curatedPCaData overview:
  https://bioconductor.posit.co/packages/3.23/data/experiment/vignettes/curatedPCaData/inst/doc/overview.html

## Scope boundary

The lung pilot and the principal colorectal confirmatory family are deliberately absent.
No expression data were downloaded and no RPP outcomes were computed for this census.
