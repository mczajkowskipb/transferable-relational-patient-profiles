# Stage-2 partial-output recovery note

**Date:** 2026-09-07  
**Original Stage-2 pre-access implementation commit:** `8a6a7d789c178a88c35d1b0808d538ccbadf512b`  
**Opaque-ID repair commit:** `191f4447b5b71ab998a6dbfbda9384e163058da3`

## Correction to the previous repair note

The previous repair note stated that the original Stage-2 run terminated before
any authentic target assignment/result CSV was written.

A subsequent resume guard detected an existing:

`GSE27262_prelabel_assignments.csv`

Therefore that stronger statement was incorrect. The historical note is not
rewritten. Instead, this correction is appended transparently.

The correct chronology is:

1. the exact pre-hashed Stage-2 target-value payload was opened;
2. the original runner generated at least the GSE27262 prelabel assignment CSV;
3. the run later encountered the valid opaque participant identifier `06L54`
   and terminated at `int(participant_id)`;
4. Stage-3 evaluation labels remained sealed;
5. the first opaque-ID repair was committed before the next attempted resume;
6. that resume intentionally stopped when it detected the pre-existing partial
   output rather than overwriting it.

Before any further endpoint computation, this recovery workflow archives every
existing file in the local `STAGE2_PRELABEL_RESULTS` directory, writes SHA-256
hashes, commits/tags that archive, and only then removes the local derived
partial result files for deterministic recomputation.

The original Stage-2 access record is retained.

No target label is used. No artifact, validation signature, C_source, target,
threshold, perturbation dose, bootstrap count or decision rule is changed.

The rerun is therefore a deterministic completion of the already-opened
development-exposed Stage-2 analysis, with the pre-existing partial state
preserved as an immutable audit record.
