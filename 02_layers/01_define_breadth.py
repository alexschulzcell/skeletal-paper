#!/usr/bin/env python3
"""Define breadth, peak strength, the shared core and the site-specific layer.

What it does
    Turns the gene x site matrix into the one number the paper is about. For
    every gene:

        breadth   how many of the six DXA sites reach Z > 2          (0-6)
        peak      the strongest single-site Z
        layer     'core'      Z > 2 at *all six* sites
                  'specific'  Z > 4 at one site and more than 2 above
                              the best of the other five
                  'rest'      everything else

    Breadth is the paper's unit of claim and is reported as a continuous
    property. The two named layers exist because a contrast needs sides, and
    because the site-specific layer is what makes the result a dissociation
    rather than a correlation; but the gene *list* is not the claim. Redefining
    the structure in the replication cohort recovers only 19 % of the same
    genes while the property replicates at t +8.3 (see
    04_analysis/05_replication_gefos.py and docs/KNOWN_ISSUES.md).

Reads
    results/gene_z_matrix.tsv     (stage 01, script 03)

Writes
    results/gene_layers.tsv       gene, six site Z, breadth, peak, min Z, layer
    results/layer_summary.tsv     size and composition of each layer
    results/breadth_distribution.tsv   genes per breadth value

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
    python 02_layers/01_define_breadth.py
    python 02_layers/01_define_breadth.py --core-z 3.0 --out gene_layers_k3.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402

#: Textbook osteoporosis-genetics genes, used as a positive control. If the
#: core does not contain these, something upstream is wrong and the printed
#: check will say so.
POSITIVE_CONTROL = ["ESR1", "LRP5", "WNT16", "WNT4", "WLS", "TNFRSF11A",
                    "TNFRSF11B", "SOX6", "CPED1", "FGFRL1", "IDUA", "MDK", "NOTUM"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--core-z", type=float, default=cfg.CORE_Z,
                    help=f"threshold defining 'acts here' (default {cfg.CORE_Z})")
    ap.add_argument("--specific-z", type=float, default=cfg.SPECIFIC_Z,
                    help=f"site-specific threshold (default {cfg.SPECIFIC_Z})")
    ap.add_argument("--specific-delta", type=float, default=cfg.SPECIFIC_DELTA,
                    help=f"margin over the next site (default {cfg.SPECIFIC_DELTA})")
    ap.add_argument("--out", default="gene_layers.tsv")
    args = ap.parse_args()

    Z = common.load_gene_matrix()
    sites = cfg.BMD_SITES

    common.banner("Two layers of skeletal genetic architecture")
    print(f"  genes: {len(Z):,}")
    print(f"  core      : Z > {args.core_z} at all six sites")
    print(f"  specific  : Z > {args.specific_z} at one site, "
          f"> {args.specific_delta} above the best other")

    out = Z[sites].copy()
    out["breadth"] = (out[sites] > args.core_z).sum(axis=1)
    out["peak"] = out[sites].max(axis=1)
    out["min_z"] = out[sites].min(axis=1)
    # Distance from the best site to the second best: the continuous form of
    # site-specificity, used in the threshold-free comparison below.
    out["specificity"] = out["peak"] - out[sites].apply(
        lambda r: sorted(r)[-2], axis=1)

    core = common.core_genes(Z, sites=sites, z=args.core_z)
    specific = common.specific_genes(Z, sites=sites, z0=args.specific_z,
                                     delta=args.specific_delta)
    out["layer"] = ["core" if g in core else "specific" if g in specific else "rest"
                    for g in out.index]
    out["chr"] = Z["chr"]
    out["gene_start"] = Z["gene_start"]
    out["gene_end"] = Z["gene_end"]

    print(f"\n{'layer':16s}{'n':>7s}{'median breadth':>16s}{'median peak':>13s}"
          f"{'median min Z':>14s}")
    summary = []
    for name in ("core", "specific", "rest"):
        s = out[out.layer == name]
        print(f"{name:16s}{len(s):7,d}{s.breadth.median():16.1f}"
              f"{s.peak.median():13.2f}{s.min_z.median():14.2f}")
        summary.append({"layer": name, "n": len(s),
                        "median_breadth": round(float(s.breadth.median()), 3),
                        "median_peak": round(float(s.peak.median()), 3),
                        "median_min_z": round(float(s.min_z.median()), 3)})

    print("\n=== Breadth distribution ===")
    dist = (out.breadth.value_counts().sort_index()
            .rename("n_genes").rename_axis("breadth").reset_index())
    dist["share"] = (dist.n_genes / len(out)).round(4)
    for _, r in dist.iterrows():
        print(f"  breadth {int(r.breadth)}: {int(r.n_genes):6,d} genes  ({r.share:5.1%})")

    print("\n=== Positive control: the textbook of osteoporosis genetics ===")
    found = [g for g in POSITIVE_CONTROL if g in core]
    absent = [g for g in POSITIVE_CONTROL if g in out.index and g not in core]
    unknown = [g for g in POSITIVE_CONTROL if g not in out.index]
    print(f"  in the core     : {len(found)}/{len(POSITIVE_CONTROL)}  {', '.join(found)}")
    if absent:
        print(f"  tested, not core: {', '.join(absent)}")
    if unknown:
        print(f"  not in matrix   : {', '.join(unknown)}")
    if len(found) < len(POSITIVE_CONTROL) // 2:
        print("  WARNING: the core does not recover the known BMD genes. "
              "Check the harmonisation step before reading anything below.")

    print("\n=== Site-specific layer, by site of the peak ===")
    spec = out[out.layer == "specific"]
    by_site = spec[sites].idxmax(axis=1).value_counts()
    for site in sites:
        print(f"  {cfg.BMD_SITE_LABELS[site]:16s}{int(by_site.get(site, 0)):4d}")

    common.write_result(out, args.out)
    common.write_result(pd.DataFrame(summary), "layer_summary.tsv", index=False)
    common.write_result(dist, "breadth_distribution.tsv", index=False)
    print("\nNow run:  python 02_layers/02_locus_clumping.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
