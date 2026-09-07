import importlib.util
from pathlib import Path
import numpy as np

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_stage2_prelabel.py"
spec=importlib.util.spec_from_file_location("s2repair",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_participant_sort_key_keeps_numeric_source_order():
    x=["102","3","91","18"]
    assert sorted(x,key=m.participant_sort_key)==["3","18","91","102"]

def test_participant_sort_key_accepts_opaque_alphanumeric_ids():
    x=["06L54","06L2","A10","9"]
    y=sorted(x,key=m.participant_sort_key)
    assert "06L54" in y and "A10" in y and y[0]=="9"

def test_group_arrays_accepts_alphanumeric_participants():
    class E:
        def __init__(self,profile):
            self.assignment="ASSIGNED"; self.best_profile_id=profile
    maprows=[
        {"gsm":"G1","participant_id":"06L54"},
        {"gsm":"G2","participant_id":"06L2"},
    ]
    pids,ai,ac,sup=m.group_arrays(
        [E("P0"),E("P1")],np.array([1.0,0.0]),maprows,["G1","G2"],"P0"
    )
    assert len(pids)==2
    assert sup["participant_assignment_coverage"]==1.0

def test_B_seed_is_repeatable_per_target():
    a=np.random.default_rng(20260916).standard_normal(8)
    b=np.random.default_rng(20260916).standard_normal(8)
    assert np.array_equal(a,b)
