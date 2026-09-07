import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "public_replay_lung_pilot.py"
spec = importlib.util.spec_from_file_location("public_replay", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_gse27262_title_parser():
    pid, lab = m.participant_and_label("GSE27262", "Lung cancer tissue sample 134T")
    assert pid == "134"
    assert lab == "tumor"

def test_gse19804_title_parser():
    pid, lab = m.participant_and_label("GSE19804", "Lung Cancer 2T")
    assert pid == "2"
    assert lab == "tumor"

def test_gse32863_title_parser():
    pid, lab = m.participant_and_label(
        "GSE32863", "Adjacent non-tumor lung tissue 05L12_N (expression)"
    )
    assert pid == "05L12"
    assert lab == "adjacent_normal"

def test_gse32863_numeric_title_parser():
    pid, lab = m.participant_and_label(
        "GSE32863", "Lung adenocarcinoma 3036_T (expression)"
    )
    assert pid == "3036"
    assert lab == "tumor"
