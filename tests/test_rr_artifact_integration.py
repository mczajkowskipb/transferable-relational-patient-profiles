"""Pre-pilot integration regressions; all observations here are synthetic."""
import unittest
from dataclasses import replace
import numpy as np
from relational_patient_profiles import fit_rr_direct, execute_rpp_artifact
from relational_patient_profiles.artifact import artifact_from_rr_direct
from relational_patient_profiles.rr_direct import (
    SparseRelationalPrototype, _score_from_prototypes, _build_prototypes,
)


class RRArtifactIntegration(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[4.,1.], [5.,2.], [1.,4.], [2.,5.]])
        self.features = ['A', 'B']
        self.model = fit_rr_direct(self.x, self.features, feature_budget=2,
                                  max_pairs=1, max_rules=1)
        self.kwargs = dict(artifact_id='regression', source_dataset='SYNTHETIC',
            feature_namespace='gene', mapping_hash='m', preprocessing_hash='p',
            config_hash='c', software_commit='test', min_score=.6,
            min_margin=.05, min_executable_coverage=.7)
        self.artifact = artifact_from_rr_direct(self.model, **self.kwargs)

    def test_binary_direction_is_converted_to_signed(self):
        self.assertEqual({r.direction for p in self.artifact.prototypes for r in p.relations}, {-1, 1})

    def test_fitted_model_and_artifact_agree(self):
        scores = _score_from_prototypes(self.x, self.features, self.model.prototypes)
        out = execute_rpp_artifact(self.x, self.features, self.artifact)
        np.testing.assert_allclose([o.profile_scores for o in out], scores)
        self.assertTrue(all(o.assignment == 'ASSIGNED' for o in out))

    def test_invalid_direction_is_rejected(self):
        rule = ('A','B',2,1.,.5)
        model = replace(self.model, prototypes=(SparseRelationalPrototype(0,(rule,)),)*2)
        with self.assertRaises(ValueError):
            artifact_from_rr_direct(model, **self.kwargs)

    def test_ties_satisfy_neither_orientation(self):
        x = np.array([[3.,3.]])
        scores = _score_from_prototypes(x, self.features, self.model.prototypes)
        np.testing.assert_array_equal(scores, [[0.,0.]])
        out = execute_rpp_artifact(x, self.features, self.artifact)[0]
        self.assertEqual(out.assignment, 'UNASSIGNED')
        self.assertEqual(out.profile_scores, (0.,0.))

    def test_ties_do_not_inflate_induction_support(self):
        # Group 0: one reverse relation and one tie. Reverse support must be .5.
        b = np.array([[0], [0], [1], [1]], dtype=np.uint8)
        rev = np.array([[1], [0], [0], [0]], dtype=np.uint8)
        p = _build_prototypes(b, np.array([0,0,1,1]), [(0,1)], ['A','B'],
                             2, 1, .8, .1, reverse_B=rev)
        self.assertEqual(p[0].rules[0][2], 0)
        self.assertEqual(p[0].rules[0][3], .5)

    def test_order_preserving_execution_invariance(self):
        shifted = self.x * np.array([.5,2.,4.,.25])[:,None] + np.array([1.,-10.,30.,5.])[:,None]
        self.assertEqual(execute_rpp_artifact(self.x,self.features,self.artifact),
                         execute_rpp_artifact(shifted,self.features,self.artifact))

    def test_target_cohort_composition_cannot_change_execution(self):
        first = execute_rpp_artifact(self.x[:1],self.features,self.artifact)[0]
        batch = np.vstack([self.x[:1], [[1e9,-1e9],[-1e9,1e9]]])
        self.assertEqual(first,execute_rpp_artifact(batch,self.features,self.artifact)[0])

    def test_feature_loss_causes_nonexecution(self):
        out = execute_rpp_artifact(self.x[:,:1], ['A'], self.artifact)
        self.assertTrue(all(o.assignment == 'UNASSIGNED' and o.reason == 'LOW_COVERAGE' for o in out))

    def test_execution_does_not_mutate_artifact(self):
        before = self.artifact.canonical_bytes()
        execute_rpp_artifact(self.x, self.features, self.artifact)
        self.assertEqual(before, self.artifact.canonical_bytes())

    def test_conversion_is_deterministic(self):
        self.assertEqual(self.artifact.sha256(), artifact_from_rr_direct(self.model, **self.kwargs).sha256())


if __name__ == '__main__':
    unittest.main()
