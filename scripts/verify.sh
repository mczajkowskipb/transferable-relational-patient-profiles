#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN="${PYTHON_BIN:-python}"
"$PYTHON_BIN" -m pytest -q
"$PYTHON_BIN" scripts/demo.py >/dev/null
python -m json.tool schemas/relational_patient_profile_v1.schema.json >/dev/null
python -m json.tool schemas/relational_transportability_certificate_v1.schema.json >/dev/null
echo "RPP clean repository verification: PASS"
