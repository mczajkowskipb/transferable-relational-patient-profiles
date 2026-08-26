import numpy as np
from relational_patient_profiles import fit_rr_direct, assign_frozen_prototypes


def test_rr_direct_is_deterministic_and_frozen_assignable():
    X = np.array([
        [4.,1.,5.,2.], [5.,1.,4.,2.], [4.,2.,6.,1.], [5.,2.,5.,1.],
        [1.,4.,2.,5.], [1.,5.,2.,4.], [2.,4.,1.,6.], [2.,5.,1.,5.],
    ])
    fids = ["A","B","C","D"]
    a = fit_rr_direct(X, fids, k=2, feature_budget=4, max_pairs=6, max_rules=4, min_support=.75, min_contrast=.05)
    b = fit_rr_direct(X, fids, k=2, feature_budget=4, max_pairs=6, max_rules=4, min_support=.75, min_contrast=.05)
    assert a == b
    assert len(set(a.labels)) == 2
    assigned, score, margin = assign_frozen_prototypes(X, fids, a.prototypes, min_score=.5, min_margin=.0)
    assert np.all(assigned >= 0)
    assert np.isfinite(score).all()
    assert np.isfinite(margin).all()
