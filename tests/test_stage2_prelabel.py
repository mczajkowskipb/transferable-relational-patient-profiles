import importlib.util
from pathlib import Path
import numpy as np

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_stage2_prelabel.py"
spec=importlib.util.spec_from_file_location("s2",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_target_bootstrap_is_deterministic_and_dataset_specific():
    a=m.bootstrap_indices(10,27262,50); b=m.bootstrap_indices(10,27262,50); c=m.bootstrap_indices(10,32863,50)
    assert np.array_equal(a,b)
    assert not np.array_equal(a,c)

def test_draw_contrast_handles_dual_group_participants():
    ai=np.array([1.0,np.nan,0.8]); ac=np.array([0.2,0.3,np.nan])
    inds=np.array([[0,1,2],[0,0,1]])
    z=m.draw_contrast(ai,ac,inds)
    assert np.isfinite(z).all()

def test_order_preserving_affine_relation_states():
    X=np.array([[1.,3.,2.],[5.,2.,4.]])
    genes=["1","2","3"]
    rel=[{"feature_a":"1","feature_b":"2","direction":"<"},
         {"feature_a":"3","feature_b":"2","direction":"<"}]
    s=m.relation_states(X,genes,rel)
    scale=np.array([2.,0.5])[:,None]; shift=np.array([7.,-3.])[:,None]
    Y=scale*X+shift
    assert np.array_equal(s,m.relation_states(Y,genes,rel))

def test_ari_perfect_under_label_permutation():
    assert abs(m.ari([0,0,1,1],[1,1,0,0])-1.0)<1e-12

def test_kmedoids_deterministic():
    Z=np.array([[0.,0.],[0.1,0.],[10.,10.],[10.1,10.]])
    a,ma=m.kmedoids(Z,["a","b","c","d"],2)
    b,mb=m.kmedoids(Z,["a","b","c","d"],2)
    assert np.array_equal(a,b) and ma==mb
