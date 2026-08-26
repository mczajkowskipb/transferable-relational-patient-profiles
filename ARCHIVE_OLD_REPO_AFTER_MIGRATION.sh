#!/usr/bin/env bash
set -euo pipefail
OWNER="${OWNER:-mczajkowskipb}"
NEW_REPO="${NEW_REPO:-transferable-relational-patient-profiles}"
OLD_REPO="${OLD_REPO:-omics-representation-audit-pilot}"
CONFIRM="${CONFIRM_ARCHIVE_OLD:-NO}"

if [[ "$CONFIRM" != "YES" ]]; then
  echo "Refusing to archive without CONFIRM_ARCHIVE_OLD=YES" >&2
  echo "Example: CONFIRM_ARCHIVE_OLD=YES bash $0" >&2
  exit 2
fi

gh auth status >/dev/null
# Require the new repository to exist first.
gh repo view "$OWNER/$NEW_REPO" >/dev/null

notice="> **Historical pilot archive.** The current SONATA BIS programme and clean RPP/RTC implementation are maintained at https://github.com/$OWNER/$NEW_REPO. This repository remains public for pilot provenance; historical GO/STOP decisions are not rewritten."

api="repos/$OWNER/$OLD_REPO/contents/README.md"
sha="$(gh api "$api" --jq .sha)"
old="$(gh api "$api" --jq .content | tr -d '\n' | base64 -d)"
if ! grep -q 'Historical pilot archive' <<<"$old"; then
  tmp="$(mktemp)"
  { printf '%s\n\n' "$notice"; printf '%s' "$old"; } > "$tmp"
  encoded="$(base64 -w0 "$tmp")"
  gh api -X PUT "$api" \
    -f message='Mark repository as historical pilot archive' \
    -f content="$encoded" \
    -f sha="$sha" >/dev/null
  rm -f "$tmp"
fi

gh repo edit "$OWNER/$OLD_REPO" \
  --description "Historical frozen evidence archive for the omics representation/RR_DIRECT pilots; current RPP programme moved to $NEW_REPO."

gh api -X PATCH "repos/$OWNER/$OLD_REPO" -F archived=true >/dev/null

echo "Archived public evidence repository: https://github.com/$OWNER/$OLD_REPO"
