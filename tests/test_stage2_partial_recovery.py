import importlib.util
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_stage2_prelabel.py"
spec=importlib.util.spec_from_file_location("s2recovery",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_recovery_runner_uses_opaque_id_sorting():
    ids=["06L54","06L2","9","A1"]
    out=sorted(ids,key=m.participant_sort_key)
    assert out[0]=="9"
    assert set(out)==set(ids)

def test_numeric_source_order_unchanged():
    ids=["102","3","18","91"]
    assert sorted(ids,key=m.participant_sort_key)==["3","18","91","102"]
