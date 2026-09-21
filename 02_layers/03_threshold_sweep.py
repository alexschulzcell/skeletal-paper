#!/usr/bin/env python3
"""Redefine both layers at every threshold, and test the claim without any.

What it does
    Two sensitivity analyses of the layer definition.

    1. The sweep. The core is rebuilt at k = 2.0, 2.5, 3.0, 3.5 and 4.0, and
       the site-specific layer at every combination of z0 in {3, 3.5, 4} and
       delta in {1.5, 2, 2.5}. Membership at each setting is written out, so
       that 04_analysis/01_enrichment.py can produce the odds ratio against
       threshold that the paper shows as Figure 1C without re-deriving the
       layers.

    2. The threshold-free version. The same question is asked with no cutoff at
       all, on three continuous axes:

           shared        the Z at the *weakest* of the six sites
           peak          the Z at the strongest site
           specificity   the distance from the best site to the second best

       The point of the paper stands or falls on which of these separates the
       Mendelian genes. That comparison needs the truth side and is made in
       04_analysis/01_enrichment.py; here the three axes are computed and
       stored.

Reads
    results/gene_z_matrix.tsv     (stage 01, script 03)

Writes
    results/threshold_layers.tsv     gene x setting membership, long format
    results/threshold_layer_sizes.tsv  layer size at each setting
    results/continuous_axes.tsv      the three threshold-free axes per gene

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Seconds.

Usage
    python 02_layers/03_threshold_sweep.py
    python 02_layers/03_threshold_sweep.py --core-thresholds 2 3 4
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402

SPECIFIC_Z_GRID = [3.0, 3.5, 4.0]
SPECIFIC_DELTA_GRID = [1.5, 2.0, 2.5]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--core-thresholds", type=float, nargs="*",
                    default=cfg.THRESHOLD_SWEEP)
    args = ap.parse_args()

    Z = common.load_gene_matrix()
    sites = cfg.BMD_SITES

    common.banner("Sensitivity of the layer definition")
    print(f"  genes: {len(Z):,}")

    rows, sizes = [], []

    print("\n=== Shared core, Z > k at all six sites ===")
    print(f"{'k':>6s}{'genes':>9s}{'share of all':>15s}")
    for k in args.core_thresholds:
        members = common.core_genes(Z, sites=sites, z=k)
        print(f"{k:6.1f}{len(members):9,d}{len(members) / len(Z):15.2%}")
        sizes.append({"layer": "core", "z0": k, "delta": None, "n": len(members)})
        rows += [{"gene": g, "layer": "core", "z0": k, "delta": None}
                 for g in sorted(members)]

    print("\n=== Site-specific layer, Z > z0 at one site and > delta above the rest ===")
    print(f"{'z0':>6s}{'delta':>7s}{'genes':>9s}")
    for z0 in SPECIFIC_Z_GRID:
        for delta in SPECIFIC_DELTA_GRID:
            members = common.specific_genes(Z, sites=sites, z0=z0, delta=delta)
            print(f"{z0:6.1f}{delta:7.1f}{len(members):9,d}")
            sizes.append({"layer": "specific", "z0": z0, "delta": delta,
                          "n": len(members)})
            rows += [{"gene": g, "layer": "specific", "z0": z0, "delta": delta}
                     for g in sorted(members)]

    print("\n=== Threshold-free axes ===")
    axes = pd.DataFrame(index=Z.index)
    axes["shared"] = Z[sites].min(axis=1)
    axes["peak"] = Z[sites].max(axis=1)
    axes["specificity"] = axes["peak"] - Z[sites].apply(lambda r: sorted(r)[-2], axis=1)
    axes["breadth"] = (Z[sites] > cfg.CORE_Z).sum(axis=1)
    for col in ("shared", "peak", "specificity"):
        print(f"  {col:14s} median {axes[col].median():+7.3f}   "
              f"IQR {axes[col].quantile(.25):+.3f} to {axes[col].quantile(.75):+.3f}")
    print("\n  Which of these separates the Mendelian genes is the question of")
    print("  the paper, and is answered in 04_analysis/01_enrichment.py.")

    common.write_result(pd.DataFrame(rows), "threshold_layers.tsv", index=False)
    common.write_result(pd.DataFrame(sizes), "threshold_layer_sizes.tsv", index=False)
    common.write_result(axes, "continuous_axes.tsv")
    print("\nNow run:  python 03_truth_sides/01_build_hpo_regions.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
