#!/usr/bin/env python3
import numpy as np
from relational_patient_profiles import (
    RPPRelation, RPPPrototype, RelationalPatientProfileArtifact,
    execute_rpp_artifact, calibrate_artifact_transportability,
)

p0 = RPPPrototype("P0", (RPPRelation("A", "B", 1, 1.0),))
p1 = RPPPrototype("P1", (RPPRelation("A", "B", -1, 1.0),))
artifact = RelationalPatientProfileArtifact(
    artifact_id="demo", source_dataset="SYNTHETIC_SOURCE", feature_namespace="gene",
    mapping_hash="demo", preprocessing_hash="demo", config_hash="demo", software_commit="demo",
    prototypes=(p0,p1), min_score=.75, min_margin=.25, min_executable_coverage=1.0,
)
X = np.array([[3.,1.],[1.,3.]])
for out in execute_rpp_artifact(X,["A","B"],artifact):
    print(out.to_dict())
Xcal = np.repeat(np.array([[3.,1.]]),29,axis=0)
print(calibrate_artifact_transportability(Xcal,["A","B"],artifact,coverage_q=.90,confidence=.95).to_dict())
