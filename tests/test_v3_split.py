import importlib.util
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"prepare_v3_split.py"
spec=importlib.util.spec_from_file_location("v3prep",SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_split_rule_deterministic():
    rec=[]
    for i in range(30):
        rec.append({"participant_id":str(i+1),"contributes_index":1 if i<20 else 0,
                    "contributes_comparator":1 if i>=5 else 0,
                    "assigned_participant":1,"n_assigned_specimens":2})
    a,_,_=m.choose_split([dict(x) for x in rec])
    b,_,_=m.choose_split([dict(x) for x in rec])
    assert a==b

def test_split_rule_preserves_sizes_and_constraints_when_feasible():
    rec=[]
    for i in range(30):
        rec.append({"participant_id":str(i+1),"contributes_index":1 if i<20 else 0,
                    "contributes_comparator":1 if i>=5 else 0,
                    "assigned_participant":1,"n_assigned_specimens":2})
    z,_,tot=m.choose_split(rec)
    assert z is not None
    cost,tup,st,ref=z
    assert len(tup)==12
    n,si,sc,sap,sas=st; ri,rc,rap,ras=ref
    assert si>=6 and sc>=6 and ri>=8 and rc>=8 and rap>=13 and ras>=26
