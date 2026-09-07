#!/usr/bin/env python3
"""Download the public GEO inputs used by the lung-pilot public replay."""

from __future__ import annotations
import argparse, csv, hashlib, shutil, sys, urllib.request
from pathlib import Path

FILES = {
    "GSE19804_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19804/matrix/GSE19804_series_matrix.txt.gz",
    "GSE27262_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27262/matrix/GSE27262_series_matrix.txt.gz",
    "GSE32863_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32863/matrix/GSE32863_series_matrix.txt.gz",
    "GPL570.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GPL6884.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/annot/GPL6884.annot.gz",
}

OPTIONAL = {
    "GSE32863_non-normalized.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32863/suppl/GSE32863_non-normalized.txt.gz",
}

FROZEN_ANNOTATION_SHA256 = {
    "GPL570.annot.gz": "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
    "GPL6884.annot.gz": "8bef75ebe5b7e28bf61cf398a31e28d63d3a9884405a3510f60ac534a0b59abb",
}

HISTORICAL_COMPRESSED_SIZE = {
    "GSE19804_series_matrix.txt.gz": 21530998,
    "GSE27262_series_matrix.txt.gz": 11944565,
    "GSE32863_series_matrix.txt.gz": 29776898,
    "GPL570.annot.gz": 8471521,
    "GPL6884.annot.gz": 6942892,
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    req = urllib.request.Request(url, headers={"User-Agent": "RPP-public-replay/1.0"})
    with urllib.request.urlopen(req, timeout=120) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    tmp.replace(dest)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="data/public_geo")
    ap.add_argument("--include-gse32863-nonnormalized", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    requested = dict(FILES)
    if args.include_gse32863_nonnormalized:
        requested.update(OPTIONAL)

    rows = []
    for name, url in requested.items():
        dest = out / name
        if dest.exists() and not args.force:
            print(f"[keep] {dest}")
        else:
            print(f"[download] {url}")
            download(url, dest)

        digest = sha256(dest)
        size = dest.stat().st_size

        if name in FROZEN_ANNOTATION_SHA256:
            expected = FROZEN_ANNOTATION_SHA256[name]
            if digest != expected:
                raise SystemExit(
                    f"ERROR: {name} SHA-256 differs from the frozen Stage-1 annotation hash.\n"
                    f"expected={expected}\nobserved={digest}\n"
                    "Stop rather than silently using a changed annotation."
                )

        hist = HISTORICAL_COMPRESSED_SIZE.get(name)
        if hist is not None and size != hist:
            print(
                f"[warning] {name} size is {size}, while the historical download snapshot was {hist}. "
                "The SHA-256 is recorded below for audit.",
                file=sys.stderr,
            )

        rows.append({
            "filename": name,
            "url": url,
            "bytes": size,
            "sha256": digest,
        })

    manifest = out / "DOWNLOAD_MANIFEST_SHA256.csv"
    with manifest.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["filename", "url", "bytes", "sha256"])
        w.writeheader()
        w.writerows(rows)

    print(f"\nManifest: {manifest}")
    for row in rows:
        print(f"{row['filename']}: {row['bytes']} bytes  {row['sha256']}")

if __name__ == "__main__":
    main()
