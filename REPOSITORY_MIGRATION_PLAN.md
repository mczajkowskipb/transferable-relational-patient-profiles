# Repository migration plan

## Recommended strategy

Use a **new clean public repository** as the grant-facing primary repository:
`mczajkowskipb/transferable-relational-patient-profiles`.

Keep `mczajkowskipb/omics-representation-audit-pilot` as an **archived public evidence repository**, not as the primary project page. Do not make it private unless there is a separate reason: public archival history strengthens provenance and lets a reviewer inspect the frozen pilot decisions.

## Why this is preferable

The old repository mixes three historical layers: Representation Audit v1, later RR_DIRECT pilot v2/v2.1, and several obsolete grant narratives. Cleaning it in place risks confusing historical statements such as v1 `direct regions: NOT TESTED` with later RR_DIRECT evidence. The new repository presents one current object hierarchy and links immutable historical evidence by exact commit/tag.

## Migration order

1. Verify this clean repository locally (`scripts/verify.sh`).
2. Create/push the new GitHub repository.
3. Confirm README, CI and file links render correctly on GitHub.
4. Change the repository URL in OSF/grant documents if a URL is explicitly entered.
5. Add a short redirect notice to the old README if desired.
6. Archive the old repository **only after** the new one is verified.

## Optional old README redirect

At the top of the historical repository README:

> **Historical pilot archive.** The current SONATA BIS programme and clean RPP/RTC implementation are maintained at `mczajkowskipb/transferable-relational-patient-profiles`. This repository remains public and immutable for pilot provenance; historical GO/STOP decisions are not rewritten.

## Do not

- delete the old repository;
- rewrite or force-push its historical tag/commits;
- change frozen pilot STOP decisions;
- move or rename untouched target outcome data into the clean repository;
- run GSE27262/GSE32863 confirmatory evaluation before the planned unseal.
