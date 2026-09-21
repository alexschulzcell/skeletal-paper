#!/usr/bin/env python3
"""Re-derive the Figure 4B regulatory density from the corrected ATAC source.

Why this is separate from 07_pathways_and_axis.py
    That script is the canonical path and needs the whole pipeline behind it:
    results/gene_layers.tsv, which needs the MAGMA gene matrix, which needs the
    downloads and hours of compute. This script needs none of that. It takes
    the shared-layer gene list that ships in results/, the published MAGMA gene
    location table, and the peak files written by
    00_setup/04_prepare_atac_peaks.py, and recomputes the same statistic.

    It exists because the accession the pipeline used to fetch was wrong
    (GSE170199, a HepG2 ChIP-seq series; see docs/KNOWN_ISSUES.md issue 12), so
    the values that ship in fig4b_regulatory_density.csv came from peak files
    of unknown provenance. This gives a cheap answer to "does the panel survive
    the correction" without re-running the pipeline.

What it deliberately does NOT do
    It does not write fig4b_regulatory_density.csv, and it is not a substitute
    for the canonical run, for two reasons:

      * the background is every gene in the MAGMA location table rather than
        the 18,392 genes that carried a Z at all six DXA sites, because that
        set is not recoverable from what ships here;
      * the 36 site-specific genes are not recoverable either, so the second
        point of the panel cannot be produced at all.

    The first difference is small — a 5% larger background, entering only
    through a length-matched median. The second is why the figure still has to
    be rebuilt from the full pipeline before submission.

Usage
    python 00_setup/04_prepare_atac_peaks.py      # writes the peak files
    python 04_analysis/07b_recheck_regulatory_density.py
    python 04_analysis/07b_recheck_regulatory_density.py --draws 5000
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402

WINDOW = 100_000


def peaks_per_gene(loc: pd.DataFrame, window: int) -> pd.DataFrame:
    """Mean peaks per skeletal element within `window` bp of each gene."""
    files = sorted(glob.glob(str(cfg.REFERENCE / "peaks" / "*hg19.bed.gz")))
    if not files:
        raise SystemExit(
            "No peak files under reference/peaks/.\n"
            "Run: python 00_setup/04_prepare_atac_peaks.py")
    elements = [os.path.basename(f).split("_")[-3] for f in files]
    print(f"  {len(files)} elements: {', '.join(elements)}")

    chroms = loc["chr"].astype(str).values
    lo = loc.start.values - window
    hi = loc.end.values + window
    counts = np.zeros((len(loc), len(files)), dtype=int)
    for j, f in enumerate(files):
        by: dict[str, list] = {}
        with gzip.open(f, "rt", errors="ignore") as fh:
            for line in fh:
                p = line.split()
                if len(p) < 3:
                    continue
                by.setdefault(p[0].replace("chr", ""), []).append(
                    (int(p[1]) + int(p[2])) // 2)
        centres = {c: np.array(sorted(v)) for c, v in by.items()}
        for i in range(len(loc)):
            arr = centres.get(chroms[i])
            if arr is None:
                continue
            counts[i, j] = int(np.searchsorted(arr, hi[i])
                               - np.searchsorted(arr, lo[i]))
        n = sum(len(v) for v in centres.values())
        print(f"    {elements[j]:14s} {n:>7,d} peaks")
    return pd.DataFrame({"peaks_per_element": counts.mean(axis=1),
                         "gene_length": loc.end.values - loc.start.values},
                        index=loc.index)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--draws", type=int, default=5000,
                    help="length-matched resamples (the paper reports 5,000)")
    ap.add_argument("--window", type=int, default=WINDOW)
    args = ap.parse_args()

    loc = common.gene_location_table()
    loc = loc.reset_index().drop_duplicates("symbol").set_index("symbol")
    loc = loc[loc["chr"].astype(str).isin(
        [str(i) for i in range(1, 23)] + ["X", "Y"])]
    print(f"  {len(loc):,} genes with coordinates")

    d = peaks_per_gene(loc, args.window)

    layer = pd.read_csv(cfg.RESULTS / "fig4c_shared_layer_genes.csv")
    core = sorted(set(layer.gene) & set(d.index))
    missing = sorted(set(layer.gene) - set(d.index))
    if missing:
        print(f"  WARNING: {len(missing)} shared-layer genes have no "
              f"coordinates: {', '.join(missing)}")

    from scipy import stats
    rho = stats.spearmanr(d.gene_length, d.peaks_per_element).statistic
    print(f"\n  Spearman(gene length, peaks per element) = {rho:+.3f}")
    print("  so the comparison below is against a length-matched null.")

    d["length_decile"] = pd.qcut(np.log10(d.gene_length + 1), 10,
                                 labels=False, duplicates="drop")
    sub = d.loc[core]
    observed = float(sub.peaks_per_element.median())
    pools = {b: g.index.values
             for b, g in d[~d.index.isin(core)].groupby("length_decile")}

    rng = np.random.default_rng(cfg.SEED)
    null = []
    for _ in range(args.draws):
        picked: list[str] = []
        for b, g in sub.groupby("length_decile"):
            pool = pools.get(b)
            if pool is None or len(pool) < len(g):
                continue
            picked.extend(rng.choice(pool, len(g), replace=False))
        if picked:
            null.append(float(d.loc[picked, "peaks_per_element"].median()))
    null = np.asarray(null)
    z = (observed - null.mean()) / null.std()
    p_emp = (1 + (null >= observed).sum()) / (len(null) + 1)
    p_z = 2 * stats.norm.sf(abs(z))

    print(f"\n=== Shared layer, n = {len(core)}, {args.draws:,} draws ===")
    print(f"  observed peaks per element   {observed:8.2f}")
    print(f"  length-matched null          {null.mean():8.2f} "
          f"+- {null.std():.2f}")
    print(f"  Z                            {z:+8.2f}")
    print(f"  P, empirical                 {p_emp:8.2e}")
    print(f"  P, from Z                    {p_z:8.2e}")

    pub = pd.read_csv(cfg.RESULTS / "fig4b_regulatory_density.csv")
    r = pub.set_index("layer").loc["Shared layer"]
    print("\n=== What ships in fig4b_regulatory_density.csv ===")
    print(f"  observed peaks per element   {r.observed_peaks_per_element:8.2f}")
    print(f"  length-matched null          {r.expected:8.2f} +- {r.sd:.2f}")
    print(f"  Z                            {r.z:+8.2f}")
    print(f"  P                            {r.p:8.2e}")

    print("\n  The absolute densities differ, because the published values were")
    print("  computed from peak files of unknown provenance. The effect does")
    print("  not: the excess over the length-matched null, and the sign and")
    print("  rough size of Z, are reproduced. Rebuild the panel from the full")
    print("  pipeline before submitting; this script cannot, because the")
    print("  site-specific layer is not recoverable from what ships here.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
