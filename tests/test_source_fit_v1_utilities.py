import importlib.util
from pathlib import Path
import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_source_fit_v1.py"
spec=importlib.util.spec_from_file_location("sourcefitv1",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_ari_permutation_invariant():
    assert abs(m.ari([0,0,1,1],[1,1,0,0])-1.0)<1e-12

def test_ari_unrelated_not_one():
    assert m.ari([0,0,1,1],[0,1,0,1]) < 1.0

def test_strict_ties_satisfy_neither():
    X=np.array([[1.,1.],[2.,1.],[1.,2.]])
    B,R=m.binary_forward_reverse(X,[(0,1)])
    assert B[:,0].tolist()==[0,1,0]
    assert R[:,0].tolist()==[0,0,1]

def test_stability_subsets_are_deterministic_and_24():
    p=[str(i) for i in range(30)]
    a=m.deterministic_subsets(p); b=m.deterministic_subsets(p)
    assert a==b and len(a)==20 and all(len(x)==24 for x in a)

def test_participant_bootstrap_degenerate_is_zero():
    vals={"INDEX":{"1":1.0},"COMPARATOR":{"2":0.0}}
    z=m.bootstrap_contrast(vals,["1"],5,123)
    assert np.all(z==0.0)
