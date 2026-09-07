import importlib.util
from pathlib import Path
import math

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_stage3_label_evaluation.py"
spec=importlib.util.spec_from_file_location("s3",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_ari_perfect_under_profile_name_swap():
    assert abs(m.ari(["P0","P0","P1","P1"],["B","B","A","A"])-1.0)<1e-12

def test_nmi_perfect_under_name_swap():
    assert abs(m.nmi(["P0","P0","P1","P1"],["B","B","A","A"])-1.0)<1e-12

def test_nmi_independent_checkerboard_zero():
    z=m.nmi(["P0","P0","P1","P1"],["A","B","A","B"])
    assert abs(z)<1e-12

def test_optimal_binary_mapping():
    mp,acc,bal=m.optimal_binary_mapping(["P0","P0","P1","P1"],["N","N","T","T"])
    assert acc==1.0 and bal==1.0 and mp["P0"]=="N"

def test_forced_profile_uses_frozen_scores_and_lexical_tie():
    row={"profile_scores":"[0.7, 0.7]"}
    assert m.forced_profile(row,["P0","P1"])=="P0"

def test_detect_column_is_case_insensitive():
    assert m.detect_column(["GSM","Evaluation_Label"],m.SAMPLE_ALIASES,"sample")=="GSM"
    assert m.detect_column(["GSM","Evaluation_Label"],m.LABEL_ALIASES,"label")=="Evaluation_Label"
