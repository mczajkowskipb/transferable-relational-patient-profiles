# Lung pilot reproducibility

## Two layers

The repository distinguishes:

### 1. Historical prospective audit trail

This layer records what was frozen before each stage of data access.

It includes:

- V1/V2/V3 protocol/configuration files;
- exact source participant splits;
- frozen assignment artifact;
- frozen validation signature;
- Stage-2 target assignments/results frozen before labels;
- SHA-256 manifests;
- the post-freeze Stage-3 descriptive label analysis;
- archived failed/partial runs.

This layer answers: **what was known and fixed at each point in the original development-exposed pilot?**

### 2. Public replay

This layer downloads the public GEO Series Matrix files and re-executes the final frozen artifact/signature.

It answers: **can an external user reconstruct the public data representation and reproduce the reported final scientific quantities?**

The public replay cannot recreate the historical state of ignorance before unsealing. That is why both layers are retained.

## Environment

Python 3.11+ is recommended.

```bash
python -m venv .venv

# Windows
.venv\Scripts\python -m pip install -U pip
.venv\Scripts\python -m pip install -e ".[test]"

# Linux/macOS
.venv/bin/python -m pip install -U pip
.venv/bin/python -m pip install -e ".[test]"
```

## Tests

```bash
python -m pytest -q
```

At the Stage-3 implementation freeze, 45 tests passed before label access. Additional public-replay helper tests are included in the repository refresh.

## Download public inputs

```bash
python scripts/download_public_geo_inputs.py --output-dir data/public_geo
```

Expected files:

```text
data/public_geo/
  GSE19804_series_matrix.txt.gz
  GSE27262_series_matrix.txt.gz
  GSE32863_series_matrix.txt.gz
  GPL570.annot.gz
  GPL6884.annot.gz
  DOWNLOAD_MANIFEST_SHA256.csv
```

## Replay the final frozen scientific objects

```bash
python scripts/public_replay_lung_pilot.py \
  --data-dir data/public_geo \
  --output-dir public_replay_outputs
```

The script must verify:

- common eligible gene universe = 17,062;
- frozen artifact SHA matches the committed artifact;
- frozen validation-signature SHA matches the committed signature;
- reconstructed source-reference C_source agrees with the committed value within numerical tolerance;
- both target structural decisions reproduce as PASS;
- target C_target/R values reproduce within numerical tolerance;
- accepted and forced label ARI/NMI reproduce from public GEO sample titles.

The script writes:

```text
public_replay_outputs/
  public_replay_summary.csv
  public_replay_checks.json
```

## Regenerate figures

```bash
python -m pip install matplotlib
python scripts/generate_preliminary_figures.py
```

Figures are written to:

`docs/preliminary/figures/`

## Historical execution scripts

The original stage-aware scripts remain in `scripts/`:

- `run_source_fit_v1.py`
- `run_source_v2.py`
- `prepare_v3_split.py`
- `run_v3_signature_reference.py`
- `run_stage2_prelabel.py`
- `run_stage3_label_evaluation.py`

They are retained as part of the audit trail. Some contain stage/tag guards that deliberately prevent them from being used as generic public replay commands.

## Key final frozen values

Source:

- C_source = 0.6053475936
- one-sided lower 95% = 0.4207264957

GSE27262:

- C_target = 0.8916666667
- D = 0.5889928699
- R = 1.4729829211
- structural decision = PASS
- conditional false PASS = 0/999
- accepted ARI = 0.9199693344
- accepted NMI = 0.8782063703

GSE32863:

- C_target = 0.6146825397
- D = 0.3120087429
- R = 1.0154208032
- structural decision = PASS
- conditional false PASS = 0/999
- accepted ARI = 0.5177621225
- accepted NMI = 0.5202361015

## Interpretation

These are development-exposed preliminary results. The source pipeline underwent transparently preserved V1→V2→V3 refinement before final target execution. The colorectal family remains reserved for independent confirmatory work.
