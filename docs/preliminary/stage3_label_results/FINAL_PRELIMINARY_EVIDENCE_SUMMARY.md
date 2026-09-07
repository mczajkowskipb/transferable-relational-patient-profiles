# Final development-exposed lung preliminary evidence

## Audit chronology

- V1 remained closed as `NO_STABLE_STRUCTURE`.
- V2 replaced the absolute stability cliff with a prospectively frozen pre-target null-calibrated stability gate; it passed that gate but stopped for insufficient signature-subset index support.
- V3 froze an assignment-stratified 12/18 holdout allocation before validation-gene analysis and then obtained `SOURCE_REFERENCE_PASS`.
- Stage 2 executed both authentic targets and froze all structural results before Stage-3 labels were opened.
- The evidence is therefore development-exposed preliminary evidence, not untouched confirmation.

## Source reference

- Frozen artifact SHA-256: `9f3eb50a5b72fcbb4997b180f11cccf51614abcb3edfb551e69ea4d313508a6b`
- Frozen validation signature: 10 relations; SHA-256 `64eb3631d6a25f0fdf0deede0e0c2da17c63fb25d82ee4612aaa8bd12470128d`
- `C_source = 0.605`; one-sided 95% lower bound `0.421`.

## Authentic cross-cohort structural retention

- **GSE27262: PASS** — C_target=0.892, D=0.589, R=1.473; assignment coverage=1.000 specimen / 1.000 participant; conditional false reassurance=0/999, exact 95% upper bound=0.0037.
- **GSE32863: PASS** — C_target=0.615, D=0.312, R=1.015; assignment coverage=0.810 specimen / 0.949 participant; conditional false reassurance=0/999, exact 95% upper bound=0.0037.

Both authentic targets therefore met the frozen structural PASS rule, and neither target produced a PASS in 999 conditional validation-block randomisations.

## Mechanistic controls

- **B / GSE27262:** frozen assignment/signature relation flips were zero at every dose; target value-space clustering agreement with the intact geometry fell as low as ARI=0.281.
- **C / GSE27262:** 0: 1/1 PASS; 0.1: 20/20 PASS; 0.25: 20/20 PASS; 0.5: 4/20 PASS; 0.75: 0/20 PASS; 1: 0/20 PASS.
- **D / GSE27262:** first non-PASS state occurred at core-gene deletion dose 0.25, with assigned fraction 0.000; complete core loss yielded no assigned samples.
- **B / GSE32863:** frozen assignment/signature relation flips were zero at every dose; target value-space clustering agreement with the intact geometry fell as low as ARI=-0.003.
- **C / GSE32863:** 0: 1/1 PASS; 0.1: 20/20 PASS; 0.25: 3/20 PASS; 0.5: 0/20 PASS; 0.75: 0/20 PASS; 1: 0/20 PASS.
- **D / GSE32863:** first non-PASS state occurred at core-gene deletion dose 0.25, with assigned fraction 0.000; complete core loss yielded no assigned samples.

## Post-freeze external-label agreement

- **GSE27262:** accepted-only ARI=0.920, NMI=0.878 on 50/50 specimens; forced ARI=0.920, NMI=0.878 on 50/50 specimens.
  Post-hoc descriptive binary mapping: accepted accuracy=0.980, balanced accuracy=0.980; forced accuracy=0.980, balanced accuracy=0.980.
- **GSE32863:** accepted-only ARI=0.518, NMI=0.520 on 94/116 specimens; forced ARI=0.471, NMI=0.500 on 116/116 specimens.
  Post-hoc descriptive binary mapping: accepted accuracy=0.862, balanced accuracy=0.841; forced accuracy=0.845, balanced accuracy=0.845.

## Interpretation ceiling

These results support feasibility of frozen source-defined relational execution, independent source-side structural contrast, and retention of that contrast in two development-exposed lung targets under explicit abstention and conditional false-reassurance controls. They do **not** establish a universal biological subtype, clinical utility, H2 family-level predictive validity, or independent colorectal confirmation.
