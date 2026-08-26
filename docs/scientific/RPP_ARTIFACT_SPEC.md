# RPPArtifact/v1 - executable scientific-object specification

The grant claim that an RPP is a reusable scientific object is represented in code by `RelationalPatientProfileArtifact/v1`.

Required provenance fields:

- artifact id and source dataset;
- feature namespace;
- mapping, preprocessing and configuration hashes;
- software commit;
- frozen profile relations and weights;
- frozen minimum score, assignment margin and executable-coverage threshold;
- explicit `created_prelabel=true` state;
- optional RTC/RTR payload.

Execution of one sample returns:

- `ASSIGNED` or `UNASSIGNED`;
- rejection reason (`LOW_COVERAGE`, `LOW_SCORE`, `LOW_MARGIN`, etc.);
- best profile id;
- score;
- margin;
- **executable coverage**;
- per-profile scores/coverage.

No target-cohort statistic or imputation is required for execution.

The transportability module adds:

- robust lower/upper profile-score bounds for a declared perturbation class;
- exact finite-breakpoint pointwise Relational Transportability Radius;
- exact distribution-free one-sided tolerance calibration on an independent `SOURCE-CALIBRATION` set;
- artifact-level RTC construction using only frozen assignments;
- metadata-only weighted platform-executability fraction derived from a frozen feature mapping.

The certificate is conditional and does not cover arbitrary structural/mixture shift. The implementation is additive and does not change the historical RR_DIRECT pilot implementation.

See `../../../examples/rpp_artifact_example.json` and `../../../scripts/31_demo_rpp_artifact.py` for a synthetic demonstration. No prospective target data are included in the example.
