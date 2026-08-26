from __future__ import annotations

import numpy as np

from relational_patient_profiles import (
    RPPRelation, RPPPrototype, RelationalPatientProfileArtifact,
    minimum_n_for_first_order_bound, distribution_free_lower_tolerance_bound,
    assignment_is_certified_at_epsilon, pointwise_transportability_radius,
    platform_executability_fraction, calibrate_profile_transportability,
    calibrate_artifact_transportability,
)


def _artifact() -> RelationalPatientProfileArtifact:
    p0 = RPPPrototype("P0", (RPPRelation("A", "B", 1, 1.0),))
    p1 = RPPPrototype("P1", (RPPRelation("A", "B", -1, 1.0),))
    return RelationalPatientProfileArtifact(
        artifact_id="a1", source_dataset="S", feature_namespace="gene",
        mapping_hash="m", preprocessing_hash="p", config_hash="c",
        software_commit="x", prototypes=(p0, p1), min_score=0.75,
        min_margin=0.25, min_executable_coverage=1.0, created_prelabel=True,
    )


def test_exact_first_order_sample_size() -> None:
    assert minimum_n_for_first_order_bound(0.90, 0.95) == 29
    no = distribution_free_lower_tolerance_bound([1.0] * 28, coverage_q=0.90, confidence=0.95)
    yes = distribution_free_lower_tolerance_bound([1.0] * 29, coverage_q=0.90, confidence=0.95)
    assert not no.certifiable
    assert yes.certifiable
    assert yes.order_rank == 1
    assert yes.lower_bound == 1.0


def test_pointwise_radius_matches_pair_margin() -> None:
    art = _artifact()
    row = np.array([3.0, 1.0])  # signed margin = 2, breakpoint epsilon=1
    assert assignment_is_certified_at_epsilon(row, ["A", "B"], art, "P0", 0.9)
    assert not assignment_is_certified_at_epsilon(row, ["A", "B"], art, "P0", 1.0)
    radius = pointwise_transportability_radius(row, ["A", "B"], art, "P0")
    assert np.isclose(radius, 1.0)


def test_profile_calibration_can_refuse_certificate() -> None:
    no = calibrate_profile_transportability([0.5] * 10, profile_id="P0", coverage_q=0.90, confidence=0.95)
    yes = calibrate_profile_transportability([0.5] * 29, profile_id="P0", coverage_q=0.90, confidence=0.95)
    assert no.transportability_radius is None
    assert yes.transportability_radius == 0.5


def test_platform_executability_is_metadata_only_fraction() -> None:
    p = RPPPrototype("P", (
        RPPRelation("A", "B", 1, 1.0),
        RPPRelation("C", "D", 1, 3.0),
    ))
    assert platform_executability_fraction(p, ["A", "B", "C"]) == 0.25
    assert platform_executability_fraction(p, ["A", "B", "C", "D"]) == 1.0


def test_artifact_level_calibration_uses_frozen_assignments() -> None:
    art = _artifact()
    X = np.repeat(np.array([[3.0, 1.0]]), 29, axis=0)
    rtc = calibrate_artifact_transportability(
        X, ["A", "B"], art, coverage_q=0.90, confidence=0.95
    )
    assert rtc.artifact_sha256 == art.sha256()
    by_id = {p.profile_id: p for p in rtc.profiles}
    assert by_id["P0"].transportability_radius == 1.0
    assert by_id["P0"].calibration_n == 29
    assert by_id["P1"].transportability_radius is None
    assert by_id["P1"].calibration_n == 0
