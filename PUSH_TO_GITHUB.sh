#!/usr/bin/env bash
set -euo pipefail
OWNER="${OWNER:-mczajkowskipb}"
NEW_REPO="${NEW_REPO:-transferable-relational-patient-profiles}"
VISIBILITY="${VISIBILITY:-public}"

gh auth status >/dev/null

if [[ -d .git ]]; then
  echo "ERROR: this clean package already contains .git; inspect it before pushing." >&2
  exit 2
fi
if [[ -z "$(git config --get user.name || true)" || -z "$(git config --get user.email || true)" ]]; then
  echo "ERROR: configure your Git identity first, for example:" >&2
  echo "  git config --global user.name 'Marcin Czajkowski'" >&2
  echo "  git config --global user.email 'YOUR_GITHUB_EMAIL'" >&2
  exit 3
fi
if gh repo view "$OWNER/$NEW_REPO" >/dev/null 2>&1; then
  echo "ERROR: https://github.com/$OWNER/$NEW_REPO already exists. Refusing to overwrite it." >&2
  exit 4
fi

PYTHONPATH=src python -m pytest -q
PYTHONPATH=src python scripts/demo.py >/dev/null
python -m json.tool schemas/relational_patient_profile_v1.schema.json >/dev/null
python -m json.tool schemas/relational_transportability_certificate_v1.schema.json >/dev/null

git init -b main
git add -A
git diff --cached --check
git commit -m "Initial clean RPP/RTC research repository"

gh repo create "$OWNER/$NEW_REPO" --"$VISIBILITY" --source=. --remote=origin --push

echo
echo "Created and pushed: https://github.com/$OWNER/$NEW_REPO"
echo "Verify README/CI/links before archiving the historical pilot repository."
