# Public replay verification - 2026-09-07

The complete final-artifact public replay was executed from freshly downloaded NCBI GEO Series Matrix files after repository refresh commit:

`cb87d0dc0edd4288f56960e58e18bb008571ddb7`

The replay returned:

**ALL CHECKS PASSED: True**

## Reproduced source quantities

- common Entrez universe: **17062**
- C_source: **0.6053475935828878**
- one-sided source lower 95%: **0.4207264957264959**

## Reproduced external targets

- **GSE27262** - decision=PASS, C_target=0.8916666667, R=1.4729829211, accepted ARI=0.9199693344, accepted NMI=0.8782063703
- **GSE32863** - decision=PASS, C_target=0.6146825397, R=1.0154208032, accepted ARI=0.5177621225, accepted NMI=0.5202361015

## Public input files used

- `GSE19804_series_matrix.txt.gz` - 21530998 bytes - SHA-256 `e96b0974243acd941824aeae9db13cbfc30663b45cbd9342bb14dad0d1372e91`
- `GSE27262_series_matrix.txt.gz` - 11944565 bytes - SHA-256 `df0a42a13be9959e4f47c85fa8bb9f291becf4f5e8cbf7e4130f7bbedba22bda`
- `GSE32863_series_matrix.txt.gz` - 29776898 bytes - SHA-256 `7f33b83dc8d2b0876e610a265f947893753b2a70854d360c247311fa3474f63f`
- `GPL570.annot.gz` - 8471521 bytes - SHA-256 `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`
- `GPL6884.annot.gz` - 6942892 bytes - SHA-256 `8bef75ebe5b7e28bf61cf398a31e28d63d3a9884405a3510f60ac534a0b59abb`

The GEO matrices themselves remain ignored and are not committed to Git.

The replay verified the final frozen artifact/signature against the committed source and target quantities. This is a reproducibility check of the development-exposed lung pilot, not a new independent confirmation.

Colorectal data remained unopened.
