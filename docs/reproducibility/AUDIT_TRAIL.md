# Immutable lung-pilot audit trail

This page indexes the principal tags in chronological order. Earlier tags are never replaced by later outcomes.

| Stage | Tag | Purpose / terminal state |
|---|---|---|
| V1 protocol | `preliminary-v1-protocol-freeze-2026-09-07` | Structural-retention V1 frozen before target access |
| V1 source implementation | `preliminary-v1-source-implementation-2026-09-07` | Source implementation frozen |
| V2 source gate | `preliminary-v2-source-gate-freeze-2026-09-07` | Null-calibrated stability rule frozen |
| V2 outcome | `preliminary-v2-source-outcome-2026-09-07` | `INSUFFICIENT_SOURCE_SIGNATURE` retained |
| V3 split rule | `preliminary-v3-stratified-split-rule-freeze-2026-09-07` | Assignment-only stratified split rule frozen |
| V3 runner repair | `preliminary-v3-runner-repair-2026-09-07` | Implementation-only precondition repair |
| V3 exact split | `preliminary-v3-source-split-freeze-2026-09-07` | Exact 12/18 split frozen before validation-gene analysis |
| V3 source package | `preliminary-v3-source-package-freeze-2026-09-07` | `SOURCE_REFERENCE_PASS` |
| Stage-2 implementation | `preliminary-stage2-implementation-freeze-2026-09-07` | Target-side implementation frozen before value access |
| Stage-2 runner repair | `preliminary-stage2-runner-repair-2026-09-07` | Opaque participant-ID repair after Stage-2 access; labels still sealed |
| Stage-2 partial recovery | `preliminary-stage2-partial-recovery-2026-09-07` | Partial output archived before deterministic rerun |
| Stage-2 prelabel freeze | `preliminary-stage2-prelabel-freeze-2026-09-07` | Complete target assignments and structural endpoints frozen before labels |
| Stage-3 implementation | `preliminary-stage3-label-evaluation-implementation-freeze-2026-09-07` | Descriptive label evaluation frozen before label access |
| Final lung pilot | `preliminary-final-lung-pilot-freeze-2026-09-07` | Final post-label development-exposed preliminary package |

## Why the negative and failed stages remain

The repository intentionally preserves:

- V1 missing its absolute stability threshold;
- V2 stopping for signature-subset support;
- implementation errors encountered during V3 and Stage 2;
- the partial Stage-2 assignment file written before one failure;
- correction notes rather than rewritten history.

This is part of the evidence, not clutter. It makes clear which changes occurred before and after each access event.

## Access boundary

At the complete Stage-2 prelabel freeze:

- target values had been opened;
- target labels had not been opened;
- both structural target decisions were already PASS;
- all target assignments and structural outputs had committed hashes.

Only after that freeze was the Stage-3 label ZIP opened.

## Confirmatory boundary

The lung family is development-exposed preliminary evidence.

The colorectal family is intentionally not opened by this pilot and remains the planned independent confirmatory family.
