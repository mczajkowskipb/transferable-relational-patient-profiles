# Stage-2 runner repair after target-value access, before target endpoint computation

**Date:** 2026-09-07  
**Original pre-access implementation commit:** `8a6a7d789c178a88c35d1b0808d538ccbadf512b`

The exact pre-hashed Stage-2 target-value payload was opened by the original
runner. The run then terminated at the first target inside `group_arrays`,
before `structural_decision` and before any authentic target assignment/result
CSV was written, because the runner attempted `int(participant_id)` and
encountered the valid opaque identifier `06L54`.

Stage-3 evaluation labels remain sealed.

This repair is implementation-only and is committed/tagged before resuming
target endpoint computation. It makes no change to the frozen artifact,
signature, C_source, structural thresholds, perturbation doses or target list.

Repairs:

1. Participant identifiers are treated as opaque strings. Purely numeric IDs
   retain the historical numeric ordering needed to reproduce the frozen source
   bootstrap. Non-numeric IDs use deterministic lexical ordering.
2. An already-opened Stage-2 extraction is not deleted or silently overwritten.
   Every extracted file is verified byte-for-byte against the exact pre-hashed
   ZIP before resumption, and the original access record is preserved.
3. The B-module RNG is corrected to the already frozen specification:
   initialise PCG64 with seed 20260916 independently for each target in lexical
   GSM order. The earlier code contained a contradictory dataset-specific
   substream line; no B result had been computed before the crash.
4. The status JSON no longer stores the hash of an intermediate version of its
   own manifest. The final manifest simply hashes every result file except the
   manifest itself.

No target outcome was used to choose any of these repairs.
