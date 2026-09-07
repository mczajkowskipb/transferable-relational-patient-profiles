import importlib.util
from pathlib import Path
import numpy as np

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_source_v2.py"
spec=importlib.util.spec_from_file_location("v2src",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_null_permutation_preserves_each_gene_values():
    X=np.arange(24,dtype=float).reshape(6,4)
    blocks=np.arange(6).reshape(3,2)
    Y=m.null_complete_discovery(X,blocks,0)
    for j in range(X.shape[1]):
        assert sorted(Y[:,j].tolist())==sorted(X[:,j].tolist())

def test_null_permutation_deterministic():
    X=np.arange(24,dtype=float).reshape(6,4)
    blocks=np.arange(6).reshape(3,2)
    assert np.array_equal(m.null_complete_discovery(X,blocks,7),
                          m.null_complete_discovery(X,blocks,7))
