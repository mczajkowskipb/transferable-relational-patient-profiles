# V1 source outcome — immutable closure

V1 source-side execution completed before target-value access.

Formal V1 outcome: **NO_STABLE_STRUCTURE**.

Recorded source-gate results:

| K | converged | stability successes | mean ARI | profile support min | Q observed | Q null95 | V1 pass |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | yes | 20/20 | 0.5895076868 | 20 | 0.1018639845 | 0.0223006541 | no |
| 3 | yes | 20/20 | 0.5547811960 | 7 | 0.0994879345 | 0.0242831257 | no |

K=2 failed only the frozen absolute stability requirement `mean ARI >= 0.60`.
K=3 additionally failed minimum profile participant support.

This V1 outcome is not changed by V2. V2 is a separately versioned, post-V1,
pre-target methodological development.
