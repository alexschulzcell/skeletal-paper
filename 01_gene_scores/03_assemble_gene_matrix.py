#!/usr/bin/env python3
"""Assemble the gene x trait matrix of MAGMA Z statistics.

What it does
    Collects MAGMA's per-trait gene output into one table indexed by gene
    symbol, with the six DXA sites first, then the clinical endpoints, then the
    replication cohort and the brain regions. A gene is kept only if it has a
    value at all six DXA sites, because breadth is undefined otherwise;
    endpoint columns may be missing and are carried as NaN so that a trait
    nobody could download does not shrink the analysis set of everything else.

    It also records how many genes hit MAGMA's gene-p floor per trait. MAGMA
    does not report a gene p below 5e-10, which truncates the top of the
    distribution: for height, 6.8 % of genes sit on the floor. Every analysis
    downstream is therefore rank-based where the top of the distribution could
    matter. See docs/KNOWN_ISSUES.md.

Reads
    $SKELBREADTH_DATA/magma/<trait>.genes.out    (stage 01, script 02)
    $SKELBREADTH_DATA/reference/NCBI37.3.gene.loc

Writes
    results/gene_z_matrix.tsv     gene x trait Z, plus chromosome and gene
                                  coordinates, which 02_layers/02 needs
    results/gene_coverage.tsv     per trait: genes, missing, on the p floor

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Under a minute.

Usage
    python 01_gene_scores/03_assemble_gene_matrix.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="gene_z_matrix.tsv",
                    help="filename under results/ (default: gene_z_matrix.tsv)")
    args = ap.parse_args()

    optional = ([t for t in cfg.ENDPOINTS] + list(cfg.GEFOS_TRAITS)
                + cfg.BRAIN_REGIONS)
    present = [t for t in optional if cfg.magma_out(t).exists()]
    absent = [t for t in optional if t not in present]

    common.banner("Assembling the gene x trait matrix")
    print(f"  required   : {', '.join(cfg.BMD_SITES)}")
    print(f"  also found : {', '.join(present) if present else '(none)'}")
    if absent:
        print(f"  not found  : {', '.join(absent)}")
        print("               the analyses that need them will skip and say so")

    Z = common.magma_matrix(cfg.BMD_SITES, extra=present)
    print(f"\n  {len(Z):,} genes with a value at all six DXA sites")

    # Coverage and the p floor, per trait.
    rows = []
    for trait in cfg.BMD_SITES + present:
        p = pd.read_csv(cfg.magma_out(trait), sep=r"\s+", index_col=0)["P"]
        on_floor = float((p <= cfg.MAGMA_P_FLOOR * 1.0001).mean())
        rows.append({
            "trait": trait,
            "label": cfg.BMD_SITE_LABELS.get(
                trait, cfg.ENDPOINTS.get(trait, cfg.GEFOS_TRAITS.get(trait, trait))),
            "N": cfg.TRAIT_N.get(trait, np.nan),
            "genes_tested": int(p.notna().sum()),
            "genes_in_matrix": int(Z[trait].notna().sum()),
            "frac_at_p_floor": round(on_floor, 4),
        })
    coverage = pd.DataFrame(rows)

    print(f"\n{'trait':12s}{'N':>12s}{'genes':>9s}{'in matrix':>11s}{'at p floor':>12s}")
    for _, r in coverage.iterrows():
        n = f"{int(r.N):,}" if pd.notna(r.N) else "-"
        print(f"{r.trait:12s}{n:>12s}{r.genes_tested:9,d}"
              f"{r.genes_in_matrix:11,d}{r.frac_at_p_floor:12.1%}")

    worst = coverage.loc[coverage.frac_at_p_floor.idxmax()]
    print(f"\n  Largest share on MAGMA's p floor of {cfg.MAGMA_P_FLOOR:g}: "
          f"{worst.trait} at {worst.frac_at_p_floor:.1%}.")
    print("  Analyses that could be affected by the truncated top use ranks.")

    common.write_result(Z, args.out)
    common.write_result(coverage, "gene_coverage.tsv", index=False)
    print("\nNow run:  python 02_layers/01_define_breadth.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
