#!/usr/bin/env python3
"""Turn the published skeletal ATAC peak workbook into per-element BED files.

The open-chromatin measure behind Figure 4B counts peaks in human embryonic
limb cartilage across eight skeletal elements. Those peaks are published as
one Excel workbook, one sheet per element, not as the per-element BED files
the analysis reads:

    Richard, D., Muthuirulan, P., Young, M., Yengo, L., Vedantam, S.,
    Marouli, E., Bartell, E., GIANT Consortium, Hirschhorn, J., and
    Capellini, T.D. (2025). Functional genomics of human skeletal development
    and the patterning of height heritability. Cell 188, 15-32.e24.
    GEO: GSE252289 (ATAC-seq), peaks called on hg19.

This script downloads the workbook and writes

    <reference>/peaks/GSE252289_ATAC_<Element>_E54_hg19.bed.gz

for the eight limb elements at embryonic day 54, which is the stage at which
all eight were profiled in the same specimens. The file name is what
`04_analysis/07_pathways_and_axis.py` expects: it globs `*hg19.bed.gz` and
takes the element from the third-from-last underscore-separated field.

    python 00_setup/04_prepare_atac_peaks.py
    python 00_setup/04_prepare_atac_peaks.py --stage 67     # the later stage
    python 00_setup/04_prepare_atac_peaks.py --keep-xlsx

Why this exists at all: the accession previously hard-coded in
`02_download_references.sh` was GSE170199, a two-sample HepG2 ChIP-seq series
from ENCODE. It is the wrong assay, the wrong tissue and the wrong species
context; see docs/KNOWN_ISSUES.md, issue 12.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

GSE = "GSE252289"
URL = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:-3]}nnn/{GSE}/suppl/"
       f"{GSE}_ATAC_peaks.xlsx")

# The eight limb elements: the proximal and distal end of each long bone.
# The workbook also holds the acetabulum, the femoral head and neck, and the
# vertebrae; the paper's measure is the limb, so those are not written.
ELEMENTS = ["DistFemur", "ProxFemur", "DistTibia", "ProxTibia",
            "DistHumerus", "ProxHumerus", "DistRadius", "ProxRadius"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", choices=["54", "67"], default="54",
                    help="embryonic day; both stages carry all eight elements")
    ap.add_argument("--keep-xlsx", action="store_true",
                    help="keep the downloaded workbook instead of deleting it")
    args = ap.parse_args()

    try:
        import pandas as pd
    except ImportError:
        print("pandas is required; activate the skeletal-breadth environment.")
        return 1

    out_dir = Path(cfg.REFERENCE) / "peaks"
    out_dir.mkdir(parents=True, exist_ok=True)
    xlsx = out_dir / f"{GSE}_ATAC_peaks.xlsx"

    if xlsx.exists():
        print(f"  workbook already present: {xlsx}")
    else:
        print(f"  downloading {URL}")
        req = urllib.request.Request(URL, headers={"User-Agent": "skeletal-breadth"})
        with urllib.request.urlopen(req, timeout=600) as r, xlsx.open("wb") as fh:
            fh.write(r.read())
        print(f"  wrote {xlsx} ({xlsx.stat().st_size / 1e6:.1f} MB)")

    written = 0
    for element in ELEMENTS:
        sheet = f"{element}_{args.stage}"
        # The sheets are headerless: the first peak sits in the header row.
        df = pd.read_excel(xlsx, sheet_name=sheet, header=None,
                           names=["chrom", "start", "end"], usecols=[0, 1, 2])
        df = df.dropna()
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        dest = out_dir / f"{GSE}_ATAC_{element}_E{args.stage}_hg19.bed.gz"
        with gzip.open(dest, "wt", newline="\n") as fh:
            for chrom, start, end in df.itertuples(index=False):
                fh.write(f"{chrom}\t{start}\t{end}\n")
        print(f"  {dest.name:48s} {len(df):>7,d} peaks")
        written += 1

    if not args.keep_xlsx:
        xlsx.unlink()

    print(f"\n  {written} element files written to {out_dir}")
    print("  Now run:  python 04_analysis/07_pathways_and_axis.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
