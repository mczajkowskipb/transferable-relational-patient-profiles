from __future__ import annotations

import numpy as np

from relational_patient_profiles import (
    RPPRelation, RPPPrototype, RelationalPatientProfileArtifact, execute_rpp_artifact,
)


def _artifact(min_cov: float = 0.75) -> RelationalPatientProfileArtifact:
    p0 = RPPPrototype("P0", (
        RPPRelation("A", "B", 1, 1.0),
        RPPRelation("C", "D", 1, 1.0),
    ))
    p1 = RPPPrototype("P1", (
        RPPRelation("A", "B", -1, 1.0),
        RPPRelation("C", "D", -1, 1.0),
    ))
    return RelationalPatientProfileArtifact(
        artifact_id="a1", source_dataset="SOURCE", feature_namespace="gene",
        mapping_hash="m", preprocessing_hash="p", config_hash="c",
        software_commit="deadbeef", prototypes=(p0, p1), min_score=0.75,
        min_margin=0.25, min_executable_coverage=min_cov, created_prelabel=True,
    )


def test_execution_is_single_sample_and_reports_coverage() -> None:
    art = _artifact()
    X = np.array([[3.0, 1.0, 4.0, 2.0], [1.0, 3.0, 2.0, 4.0]])
    out = execute_rpp_artifact(X, ["A", "B", "C", "D"], art)
    assert [x.assignment for x in out] == ["ASSIGNED", "ASSIGNED"]
    assert [x.best_profile_id for x in out] == ["P0", "P1"]
    assert out[0].executable_coverage == 1.0
    assert out[1].executable_coverage == 1.0


def test_missing_features_trigger_unassigned_without_imputation() -> None:
    art = _artifact(min_cov=0.75)
    X = np.array([[3.0, 1.0, np.nan, np.nan]])
    out = execute_rpp_artifact(X, ["A", "B", "C", "D"], art)[0]
    assert out.assignment == "UNASSIGNED"
    assert out.reason == "LOW_COVERAGE"
    assert out.executable_coverage == 0.5


def test_artifact_serialisation_is_deterministic() -> None:
    art = _artifact()
    assert art.canonical_bytes() == art.canonical_bytes()
    assert art.sha256() == art.sha256()
    assert len(art.sha256()) == 64
