# Frozen Pilot v2 / v2.1 evidence snapshot

This clean repository does **not** rewrite the historical pilot. The authoritative reproducible archive remains:

- repository: `mczajkowskipb/omics-representation-audit-pilot`
- current evidence-alignment commit: `95ddbbdb962bc8519d59b148c27e88baf84d45e5`
- historical closeout tag: `pilot-closeout-2026-08-17`

## Pilot v2 formal status

**STOP (4/5 prospective criteria passed).** The failed criterion was exact designated-pair recovery in a non-identifiable generator. Thresholds were not relaxed after observing the result.

Key frozen observations:
- synthetic RR_DIRECT median ARI: **1.0000**;
- designated-pair recovery in the non-identifiable generator: **0.1667**;
- real RR_DIRECT median ARI across 11 datasets: **0.0332**;
- real RELATION/PAM median ARI: **0.0223**;
- best frozen transfer ARI: **0.9596**, coverage **0.9252**;
- opposite frozen direction: ARI **0.771**, coverage **0.967**.

## Pilot v2.1 diagnostic

The identifiability addendum is diagnostic only and does **not** convert v2 to GO. With unchanged hyperparameters and an identifiable separated-order generator, 120 replicates gave median source ARI, exact-pair recovery, target ARI and coverage all equal to **1.0**.

## Interpretation boundary

The low median real-data ARI is evidence against broad superiority. The strong frozen transfer directions are feasibility observations, not definitive validation of the future RPP/RTC framework. The SONATA BIS project therefore asks **when** a reusable relational profile is supported and transportable, with `NO_STABLE_STRUCTURE`, `UNASSIGNED`, and `NOT_CERTIFIABLE` as valid negative outcomes.
