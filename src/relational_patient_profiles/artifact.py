"""Versioned executable Relational Patient Profile artifacts.

This module is additive: it does not modify the frozen RR_DIRECT pilot implementation.
It turns learned sparse relational prototypes into a stable, provenance-rich object
that can be executed on one unseen sample without target-cohort fitting.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence
import hashlib
import json
import math
import numpy as np


RPP_SCHEMA = "RelationalPatientProfileArtifact/v1"
RPP_EXECUTION_SCHEMA = "RPPExecution/v1"


@dataclass(frozen=True, slots=True)
class RPPRelation:
    feature_a: str
    feature_b: str
    direction: int
    weight: float
    within_support: float | None = None
    contrast: float | None = None

    def __post_init__(self) -> None:
        if self.direction not in (-1, 1):
            raise ValueError("direction must be -1 or +1")
        if not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("weight must be finite and >0")
        if self.feature_a == self.feature_b:
            raise ValueError("relation requires two different features")

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "feature_a": self.feature_a,
            "feature_b": self.feature_b,
            "direction": ">" if self.direction == 1 else "<",
            "weight": self.weight,
        }
        if self.within_support is not None:
            out["within_support"] = self.within_support
        if self.contrast is not None:
            out["contrast"] = self.contrast
        return out


@dataclass(frozen=True, slots=True)
class RPPPrototype:
    profile_id: str
    relations: tuple[RPPRelation, ...]

    def __post_init__(self) -> None:
        if not self.profile_id:
            raise ValueError("profile_id is required")
        if not self.relations:
            raise ValueError("profile must contain at least one relation")

    def to_dict(self) -> dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "relations": [r.to_dict() for r in self.relations],
        }


@dataclass(frozen=True, slots=True)
class RelationalPatientProfileArtifact:
    artifact_id: str
    source_dataset: str
    feature_namespace: str
    mapping_hash: str
    preprocessing_hash: str
    config_hash: str
    software_commit: str
    prototypes: tuple[RPPPrototype, ...]
    min_score: float
    min_margin: float
    min_executable_coverage: float
    created_prelabel: bool = True
    transportability_certificate: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if not self.artifact_id or not self.source_dataset:
            raise ValueError("artifact_id and source_dataset are required")
        if len(self.prototypes) < 2:
            raise ValueError("at least two profiles are required")
        for name, value in (
            ("min_score", self.min_score),
            ("min_margin", self.min_margin),
            ("min_executable_coverage", self.min_executable_coverage),
        ):
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if not 0 <= self.min_score <= 1:
            raise ValueError("min_score must be in [0,1]")
        if not 0 <= self.min_executable_coverage <= 1:
            raise ValueError("min_executable_coverage must be in [0,1]")
        if self.min_margin < 0:
            raise ValueError("min_margin must be >=0")
        if not self.created_prelabel:
            raise ValueError("confirmatory RPP artifacts must be created prelabel")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": RPP_SCHEMA,
            "artifact_id": self.artifact_id,
            "source_dataset": self.source_dataset,
            "feature_namespace": self.feature_namespace,
            "mapping_hash": self.mapping_hash,
            "preprocessing_hash": self.preprocessing_hash,
            "config_hash": self.config_hash,
            "software_commit": self.software_commit,
            "min_score": self.min_score,
            "min_margin": self.min_margin,
            "min_executable_coverage": self.min_executable_coverage,
            "created_prelabel": self.created_prelabel,
            "prototypes": [p.to_dict() for p in self.prototypes],
            "transportability_certificate": self.transportability_certificate,
        }

    def canonical_bytes(self) -> bytes:
        return (json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_bytes()).hexdigest()


@dataclass(frozen=True, slots=True)
class RPPExecution:
    assignment: str
    reason: str
    best_profile_id: str | None
    best_score: float | None
    margin: float | None
    executable_coverage: float | None
    profile_scores: tuple[float | None, ...]
    profile_coverages: tuple[float, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": RPP_EXECUTION_SCHEMA,
            "assignment": self.assignment,
            "reason": self.reason,
            "best_profile_id": self.best_profile_id,
            "best_score": self.best_score,
            "margin": self.margin,
            "executable_coverage": self.executable_coverage,
            "profile_scores": self.profile_scores,
            "profile_coverages": self.profile_coverages,
        }


def _profile_score_coverage(
    row: np.ndarray,
    index: dict[str, int],
    profile: RPPPrototype,
) -> tuple[float | None, float]:
    numerator = 0.0
    executable_weight = 0.0
    total_weight = sum(r.weight for r in profile.relations)
    for rel in profile.relations:
        if rel.feature_a not in index or rel.feature_b not in index:
            continue
        xa = row[index[rel.feature_a]]
        xb = row[index[rel.feature_b]]
        if not (np.isfinite(xa) and np.isfinite(xb)):
            continue
        executable_weight += rel.weight
        observed = xa > xb
        satisfied = observed if rel.direction == 1 else (not observed and xa != xb)
        if satisfied:
            numerator += rel.weight
    coverage = executable_weight / total_weight if total_weight > 0 else 0.0
    if executable_weight <= 0:
        return None, coverage
    return numerator / executable_weight, coverage


def execute_rpp_artifact(
    X: np.ndarray,
    feature_ids: Sequence[str],
    artifact: RelationalPatientProfileArtifact,
) -> tuple[RPPExecution, ...]:
    """Execute a frozen RPP artifact independently for every target sample.

    No target-cohort statistic, imputation or refitting is used. Missing features
    reduce executable coverage. A sample can remain UNASSIGNED.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != len(feature_ids):
        raise ValueError("feature_ids must match X columns")
    index = {str(f): j for j, f in enumerate(feature_ids)}
    if len(index) != len(feature_ids):
        raise ValueError("feature_ids must be unique")

    outputs: list[RPPExecution] = []
    for row in X:
        scores: list[float | None] = []
        coverages: list[float] = []
        eligible: list[int] = []
        for j, proto in enumerate(artifact.prototypes):
            score, cov = _profile_score_coverage(row, index, proto)
            scores.append(score)
            coverages.append(cov)
            if score is not None and cov >= artifact.min_executable_coverage:
                eligible.append(j)

        if not eligible:
            outputs.append(RPPExecution(
                assignment="UNASSIGNED", reason="LOW_COVERAGE",
                best_profile_id=None, best_score=None, margin=None,
                executable_coverage=max(coverages, default=0.0),
                profile_scores=tuple(scores), profile_coverages=tuple(coverages),
            ))
            continue

        ranked = sorted(eligible, key=lambda j: (-float(scores[j]), artifact.prototypes[j].profile_id))
        best = ranked[0]
        best_score = float(scores[best])
        best_cov = float(coverages[best])
        if best_score < artifact.min_score:
            outputs.append(RPPExecution(
                assignment="UNASSIGNED", reason="LOW_SCORE",
                best_profile_id=artifact.prototypes[best].profile_id,
                best_score=best_score, margin=None, executable_coverage=best_cov,
                profile_scores=tuple(scores), profile_coverages=tuple(coverages),
            ))
            continue

        if len(ranked) < 2:
            outputs.append(RPPExecution(
                assignment="UNASSIGNED", reason="INSUFFICIENT_COMPARABLE_PROFILES",
                best_profile_id=artifact.prototypes[best].profile_id,
                best_score=best_score, margin=None, executable_coverage=best_cov,
                profile_scores=tuple(scores), profile_coverages=tuple(coverages),
            ))
            continue

        second = ranked[1]
        margin = best_score - float(scores[second])
        if margin < artifact.min_margin:
            outputs.append(RPPExecution(
                assignment="UNASSIGNED", reason="LOW_MARGIN",
                best_profile_id=artifact.prototypes[best].profile_id,
                best_score=best_score, margin=margin, executable_coverage=best_cov,
                profile_scores=tuple(scores), profile_coverages=tuple(coverages),
            ))
            continue

        outputs.append(RPPExecution(
            assignment="ASSIGNED", reason="OK",
            best_profile_id=artifact.prototypes[best].profile_id,
            best_score=best_score, margin=margin, executable_coverage=best_cov,
            profile_scores=tuple(scores), profile_coverages=tuple(coverages),
        ))
    return tuple(outputs)


def artifact_from_rr_direct(
    rr_result: Any,
    *,
    artifact_id: str,
    source_dataset: str,
    feature_namespace: str,
    mapping_hash: str,
    preprocessing_hash: str,
    config_hash: str,
    software_commit: str,
    min_score: float,
    min_margin: float,
    min_executable_coverage: float,
) -> RelationalPatientProfileArtifact:
    """Convert a frozen RR_DIRECT result without importing or changing pilot code."""
    protos: list[RPPPrototype] = []
    for p in rr_result.prototypes:
        relations: list[RPPRelation] = []
        for a, b, direction, support, contrast in p.rules:
            weight = max(1e-12, float(support) * max(float(contrast), 1e-6))
            relations.append(RPPRelation(
                feature_a=str(a), feature_b=str(b), direction=int(direction),
                weight=weight, within_support=float(support), contrast=float(contrast),
            ))
        protos.append(RPPPrototype(profile_id=f"P{int(p.cluster)}", relations=tuple(relations)))
    return RelationalPatientProfileArtifact(
        artifact_id=artifact_id,
        source_dataset=source_dataset,
        feature_namespace=feature_namespace,
        mapping_hash=mapping_hash,
        preprocessing_hash=preprocessing_hash,
        config_hash=config_hash,
        software_commit=software_commit,
        prototypes=tuple(protos),
        min_score=min_score,
        min_margin=min_margin,
        min_executable_coverage=min_executable_coverage,
        created_prelabel=True,
    )
