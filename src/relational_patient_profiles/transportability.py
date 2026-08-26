"""Source-only conditional transportability certificates for frozen RPPs.

The certificate is deliberately model-relative. It guarantees assignment
preservation only for an explicitly declared bounded perturbation class around
an independently calibrated source-profile population. It is not a guarantee
for arbitrary real-world target cohorts.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence
import math
import numpy as np
from scipy.stats import binom

from .artifact import RelationalPatientProfileArtifact, RPPPrototype


RTC_SCHEMA = "RelationalTransportabilityCertificate/v1"


@dataclass(frozen=True, slots=True)
class ToleranceBound:
    coverage_q: float
    confidence: float
    n: int
    order_rank: int | None
    lower_bound: float | None

    @property
    def certifiable(self) -> bool:
        return self.order_rank is not None and self.lower_bound is not None


@dataclass(frozen=True, slots=True)
class ProfileTransportabilityCertificate:
    profile_id: str
    perturbation_class: str
    coverage_q: float
    confidence: float
    calibration_n: int
    transportability_radius: float | None
    order_rank: int | None
    assumptions: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "perturbation_class": self.perturbation_class,
            "coverage_q": self.coverage_q,
            "confidence": self.confidence,
            "calibration_n": self.calibration_n,
            "transportability_radius": self.transportability_radius,
            "order_rank": self.order_rank,
            "assumptions": list(self.assumptions),
            "status": "CERTIFIED" if self.transportability_radius is not None else "NOT_CERTIFIABLE",
        }


@dataclass(frozen=True, slots=True)
class RelationalTransportabilityCertificate:
    source_dataset: str
    artifact_sha256: str
    calibration_policy: str
    profiles: tuple[ProfileTransportabilityCertificate, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": RTC_SCHEMA,
            "source_dataset": self.source_dataset,
            "artifact_sha256": self.artifact_sha256,
            "calibration_policy": self.calibration_policy,
            "profiles": [p.to_dict() for p in self.profiles],
        }


def minimum_n_for_first_order_bound(coverage_q: float, confidence: float) -> int:
    """Minimum n for any nonparametric one-sided lower tolerance bound.

    For the first order statistic, the exact condition is q**n <= 1-confidence.
    """
    if not 0 < coverage_q < 1 or not 0 < confidence < 1:
        raise ValueError("coverage_q and confidence must be in (0,1)")
    alpha = 1.0 - confidence
    return int(math.ceil(math.log(alpha) / math.log(coverage_q)))


def distribution_free_lower_tolerance_bound(
    values: Sequence[float],
    *,
    coverage_q: float,
    confidence: float,
) -> ToleranceBound:
    """Exact one-sided nonparametric lower tolerance bound via order statistics.

    Returns X_(r), choosing the largest r such that
    P[Binomial(n, 1-q) >= r] >= confidence.
    Then, under IID calibration sampling from a fixed population and continuous
    distribution (ties make the result conservative), with at least `confidence`
    probability, at least `coverage_q` of that population lies at or above X_(r).
    """
    if not 0 < coverage_q < 1 or not 0 < confidence < 1:
        raise ValueError("coverage_q and confidence must be in (0,1)")
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n == 0:
        return ToleranceBound(coverage_q, confidence, 0, None, None)
    p = 1.0 - coverage_q
    valid: list[int] = []
    for r in range(1, n + 1):
        if float(binom.sf(r - 1, n, p)) >= confidence:
            valid.append(r)
    if not valid:
        return ToleranceBound(coverage_q, confidence, n, None, None)
    r = max(valid)
    xs = np.sort(x)
    return ToleranceBound(coverage_q, confidence, n, r, float(xs[r - 1]))


def _oriented_margin(row: np.ndarray, index: dict[str, int], rel: Any) -> float | None:
    if rel.feature_a not in index or rel.feature_b not in index:
        return None
    xa = row[index[rel.feature_a]]
    xb = row[index[rel.feature_b]]
    if not (np.isfinite(xa) and np.isfinite(xb)):
        return None
    return float(rel.direction) * float(xa - xb)


def robust_score_bounds(
    row: np.ndarray,
    feature_ids: Sequence[str],
    profile: RPPPrototype,
    epsilon: float,
) -> tuple[float | None, float | None, float]:
    """Worst-case score bounds under |delta_i|<=epsilon for observed features.

    Lower score counts only relations guaranteed satisfied (margin > 2epsilon).
    Upper score counts relations that are not guaranteed unsatisfied
    (margin > -2epsilon). Missing/unmappable relations reduce coverage and are
    excluded from the executable denominator, matching frozen execution.
    """
    if epsilon < 0:
        raise ValueError("epsilon must be >=0")
    row = np.asarray(row, dtype=float)
    index = {str(f): j for j, f in enumerate(feature_ids)}
    total_weight = sum(r.weight for r in profile.relations)
    executable_weight = 0.0
    lower_num = 0.0
    upper_num = 0.0
    for rel in profile.relations:
        m = _oriented_margin(row, index, rel)
        if m is None:
            continue
        executable_weight += rel.weight
        if m > 2.0 * epsilon:
            lower_num += rel.weight
        if m > -2.0 * epsilon:
            upper_num += rel.weight
    coverage = executable_weight / total_weight if total_weight > 0 else 0.0
    if executable_weight <= 0:
        return None, None, coverage
    return lower_num / executable_weight, upper_num / executable_weight, coverage


def assignment_is_certified_at_epsilon(
    row: np.ndarray,
    feature_ids: Sequence[str],
    artifact: RelationalPatientProfileArtifact,
    assigned_profile_id: str,
    epsilon: float,
) -> bool:
    bounds = [robust_score_bounds(row, feature_ids, p, epsilon) for p in artifact.prototypes]
    ids = [p.profile_id for p in artifact.prototypes]
    if assigned_profile_id not in ids:
        raise ValueError("assigned_profile_id not in artifact")
    k = ids.index(assigned_profile_id)
    lower_k, _, cov_k = bounds[k]
    if lower_k is None or cov_k < artifact.min_executable_coverage:
        return False
    if lower_k < artifact.min_score:
        return False
    competitor_uppers = [b[1] for j, b in enumerate(bounds) if j != k and b[1] is not None and b[2] >= artifact.min_executable_coverage]
    if not competitor_uppers:
        return False
    robust_margin = lower_k - max(float(x) for x in competitor_uppers)
    return robust_margin >= artifact.min_margin


def pointwise_transportability_radius(
    row: np.ndarray,
    feature_ids: Sequence[str],
    artifact: RelationalPatientProfileArtifact,
    assigned_profile_id: str,
) -> float:
    """Largest certified symmetric additive perturbation radius for one sample.

    The exact robust conditions change only when epsilon crosses |margin|/2, so
    the search is over finite breakpoints rather than an arbitrary numeric grid.
    """
    row = np.asarray(row, dtype=float)
    index = {str(f): j for j, f in enumerate(feature_ids)}
    breaks = {0.0}
    for profile in artifact.prototypes:
        for rel in profile.relations:
            m = _oriented_margin(row, index, rel)
            if m is not None:
                breaks.add(abs(m) / 2.0)
    ordered = sorted(breaks)
    best = 0.0
    # Conditions use strict margin inequalities; evaluate just below each break.
    for b in ordered[1:]:
        eps = np.nextafter(float(b), 0.0)
        if assignment_is_certified_at_epsilon(row, feature_ids, artifact, assigned_profile_id, eps):
            best = float(b)
        else:
            break
    return best



def platform_executability_fraction(
    profile: RPPPrototype,
    available_feature_ids: Sequence[str],
) -> float:
    """Metadata-only weighted executability upper bound for one platform.

    `available_feature_ids` must be derived from a frozen platform annotation /
    mapping table, not from target expression values or outcome labels.
    """
    available = {str(x) for x in available_feature_ids}
    total = sum(r.weight for r in profile.relations)
    if total <= 0:
        return 0.0
    executable = sum(
        r.weight for r in profile.relations
        if r.feature_a in available and r.feature_b in available
    )
    return float(executable / total)


def calibrate_artifact_transportability(
    X_calibration: np.ndarray,
    feature_ids: Sequence[str],
    artifact: RelationalPatientProfileArtifact,
    *,
    coverage_q: float = 0.80,
    confidence: float = 0.95,
) -> RelationalTransportabilityCertificate:
    """Build a profile-level RTC from an independent SOURCE-CALIBRATION set.

    The artifact MUST already be frozen. Calibration samples are executed
    independently by that artifact. For each frozen profile, the certificate
    population is the source-calibration population conditionally assigned to
    that profile by the frozen assignment rule. UNASSIGNED samples do not create
    a profile certificate and are reported separately by the caller if needed.
    """
    # Local import avoids a circular import at module import time.
    from .artifact import execute_rpp_artifact

    X = np.asarray(X_calibration, dtype=float)
    executions = execute_rpp_artifact(X, feature_ids, artifact)
    radii_by_profile: dict[str, list[float]] = {
        p.profile_id: [] for p in artifact.prototypes
    }
    for row, execution in zip(X, executions):
        if execution.assignment != "ASSIGNED" or execution.best_profile_id is None:
            continue
        pid = execution.best_profile_id
        radii_by_profile[pid].append(
            pointwise_transportability_radius(row, feature_ids, artifact, pid)
        )

    profiles = tuple(
        calibrate_profile_transportability(
            radii_by_profile[p.profile_id],
            profile_id=p.profile_id,
            coverage_q=coverage_q,
            confidence=confidence,
        )
        for p in artifact.prototypes
    )
    return RelationalTransportabilityCertificate(
        source_dataset=artifact.source_dataset,
        artifact_sha256=artifact.sha256(),
        calibration_policy=(
            "independent SOURCE-CALIBRATION; frozen RPP execution; "
            "profile population conditional on ASSIGNED status"
        ),
        profiles=profiles,
    )

def calibrate_profile_transportability(
    radii: Sequence[float],
    *,
    profile_id: str,
    coverage_q: float = 0.80,
    confidence: float = 0.95,
) -> ProfileTransportabilityCertificate:
    """Calibrate an exact distribution-free profile-level radius.

    `radii` MUST come from an independent source-calibration subset after the
    profile artifact and thresholds are fixed. This avoids post-selection use of
    the same samples that created the relational profile.
    """
    tb = distribution_free_lower_tolerance_bound(
        radii, coverage_q=coverage_q, confidence=confidence
    )
    radius = None if tb.lower_bound is None else max(0.0, float(tb.lower_bound))
    return ProfileTransportabilityCertificate(
        profile_id=profile_id,
        perturbation_class="symmetric_featurewise_additive_Linf",
        coverage_q=coverage_q,
        confidence=confidence,
        calibration_n=tb.n,
        transportability_radius=radius,
        order_rank=tb.order_rank,
        assumptions=(
            "profile and thresholds fixed before calibration",
            "calibration samples independent of profile fitting",
            "IID calibration sampling within the source-profile population",
            "feature-wise additive perturbation satisfies |delta_i|<=epsilon",
            "certificate does not cover structural/mixture shift or unavailable features",
        ),
    )
