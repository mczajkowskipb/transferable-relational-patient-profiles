# Public GEO data access

No raw or processed GEO expression matrix is redistributed in this repository.

The lung pilot uses three public NCBI GEO Series and two GEO platform annotations.

## Cohorts

| Role | GEO Series | Platform | Public accession |
|---|---|---|---|
| Source | GSE19804 | GPL570 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804 |
| External target 1 | GSE27262 | GPL570 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27262 |
| External target 2 | GSE32863 | GPL6884 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32863 |

GSE27262 contains 50 samples from 25 paired lung adenocarcinoma / adjacent-normal participants. GSE32863 contains 116 expression samples on GPL6884. GSE19804 contains 120 samples and is used as the source cohort.

## Direct files used by the public replay

Series Matrix:

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27262/matrix/GSE27262_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32863/matrix/GSE32863_series_matrix.txt.gz`

Platform annotations:

- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/annot/GPL6884.annot.gz`

Optional historical supplementary file, not required by the public final-artifact replay:

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32863/suppl/GSE32863_non-normalized.txt.gz`

## Download automatically

```bash
python scripts/download_public_geo_inputs.py --output-dir data/public_geo
```

The downloader records file sizes and SHA-256 values after download.

The two platform annotations are checked against the hashes preserved in the original Stage-1 manifest:

- GPL570 annotation: `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`
- GPL6884 annotation: `8bef75ebe5b7e28bf61cf398a31e28d63d3a9884405a3510f60ac534a0b59abb`

If NCBI replaces one of those exact annotation files, an exact historical replay should stop and report the mismatch rather than silently continuing.

## Feature mapping

The frozen preprocessing rule is:

1. read the GEO platform annotation;
2. retain a probe only when `Gene ID` contains exactly one numeric Entrez Gene ID;
3. discard missing/non-numeric/multi-ID mappings;
4. if multiple probes map to the same Entrez gene, aggregate by the within-specimen median;
5. define the cross-platform universe as the intersection of eligible GPL570 and GPL6884 Entrez IDs;
6. sort the common universe numerically.

For the frozen pilot this produced:

- GPL570 eligible genes: 20,848
- GPL6884 eligible genes: 19,405
- common cross-platform universe: **17,062 genes**

The public replay verifies that the common universe is exactly 17,062 before executing the frozen artifact.

## Participant identifiers and labels

For replay purposes, participant IDs are reconstructed from the public GEO sample titles.

The final public replay derives `tumor` / `adjacent_normal` only after the structural endpoint has been computed. This does not recreate the original prospective blinding event; that event is preserved by the historical commits/tags and Stage-2 prelabel hashes.

## What is and is not reproducible from public data

Two complementary forms of reproducibility are provided:

1. **Historical audit reproducibility** — frozen protocols, scripts, intermediate outcomes, hashes, commits and tags show what was fixed before each access stage.
2. **Public scientific replay** — public GEO files are downloaded again and the final frozen artifact/signature are re-executed to verify the final source/target metrics.

The public replay is intentionally not described as a new independent confirmation.
