# SONATA BIS 16 - prospective external-cohort freeze v2

**Freeze timestamp:** 2026-08-25 19:50 CEST
**Status:** metadata-only cohort selection; no confirmatory RR_DIRECT outcome analysis has been run on the frozen targets.

## Primary confirmatory module - lung adenocarcinoma

### Source: GSE19804
Previously used in the pilot; frozen discovery/source cohort for future confirmatory transfer.

### Untouched target 1: GSE27262
- 25 stage-I lung adenocarcinoma tumour/adjacent-normal pairs (50 samples);
- Affymetrix HG-U133 Plus 2.0 (GPL570);
- same-platform confirmatory target;
- absent from the frozen 11-dataset pilot;
- no RR_DIRECT transfer result inspected before freeze.

### Untouched target 2: GSE32863
- 58 lung adenocarcinoma and 58 matched adjacent non-tumour tissues;
- Illumina HumanWG-6 v3.0 (GPL6884);
- deliberate cross-platform confirmatory target;
- metadata previously listed in proposal-stage audit, but no RR_DIRECT confirmatory outcome evaluation or target-guided tuning.

## Frozen G5 criterion
A single source artifact fitted on GSE19804 must be applied unchanged to both targets.

For each target:
1. executable coverage >= 0.70;
2. forced all-sample ARI from frozen prototype scores >= 0.50;
3. assigned-sample ARI/NMI, rejection rate, score, margin and executable-rule coverage are secondary diagnostics;
4. every target sample and failed direction is reported.

**G5 PASS:** both targets satisfy criteria 1-2.
**G5 FAIL:** either target fails. A failed target cannot be replaced after evaluation labels are opened.

## Mapping freeze
- eligible common gene universe defined before source fitting from platform annotation only;
- target expression values/labels cannot determine the universe;
- identifier namespace, annotation release, multi-probe aggregation rule and mapping hashes frozen before source fitting;
- practical default: within-sample median aggregation when multiple platform features map to one gene;
- low target coverage causes UNASSIGNED/FAIL rather than target-guided replacement.

## Secondary colorectal module
- source: GSE39582
- target 1: GSE14333
- target 2: GSE33113

This module cannot rescue a failed primary lung G5. Exact subtype/phenotype endpoint is finalised only after sample-level metadata audit and before labels are opened.

## Label firewall
1. create sample/feature manifest;
2. keep evaluation labels outside fit/assignment namespace;
3. fit source artifact and write target assignments, scores, coverage and hashes;
4. seal prelabel artifact;
5. only then open phenotype/subtype labels;
6. no rerun may modify the source artifact after unsealing.

## Prohibited
- joint source-target normalization;
- target-dependent batch correction;
- target-performance-guided mapping/feature/relation selection;
- target-driven K or threshold tuning;
- target reclustering;
- replacing failed target after unsealing;
- reporting only favourable directions.
