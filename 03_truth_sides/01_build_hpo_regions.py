#!/usr/bin/env python3
"""Build the disease side: six anatomical body regions from the HPO.

What it does
    Mirrors the anatomy of the GWAS side on the disease side. Six body regions
    - skull, spine, thorax, pelvis, long bones, hands and feet - are defined as
    the full descendant sets of named HPO terms (the roots are listed in
    00_setup/config.py). For every gene with a Mendelian disease annotation,
    the script counts how many of the six regions its phenotype touches, and
    records the total number of HPO terms annotated to it.

    That total is the ascertainment control, and it is why the target is a
    *proportion* - regions out of six - rather than a count of terms. A
    well-studied gene accumulates terms; it does not thereby accumulate body
    regions. The correlation between breadth and total term count comes out at
    rho = -0.007, which is the number the paper puts in the figure.

    One bug is worth naming, because it inverted a conclusion once: the term
    count must be computed for *every* gene with any HPO disease, not only for
    the genes of the organ under test. Restricting it silently shrank the
    skeletal regression from 2,813 to 1,860 genes and moved p from 0.0012 to
    0.117. This script computes it over all genes, and prints both counts so
    the difference is visible.

Reads
    $SKELBREADTH_DATA/reference/hp.obo
    $SKELBREADTH_DATA/reference/phenotype.hpoa
    $SKELBREADTH_DATA/reference/genes_to_disease.txt

Writes
    results/truth_hpo.tsv         gene, n_hpo_terms, regions_skeletal,
                                  has_skeletal_disease, and one column per region
    results/hpo_region_sizes.tsv  how many terms each region contains

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    One to two minutes; dominated by parsing hp.obo.

Usage
    python 03_truth_sides/01_build_hpo_regions.py
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
    ap.add_argument("--out", default="truth_hpo.tsv")
    args = ap.parse_args()

    common.banner("Truth side 1 of 4: HPO body regions")

    regions = common.hpo_descendants(cfg.HPO_SKELETAL_REGIONS)
    print("\n=== Six skeletal body regions ===")
    sizes = []
    for name, terms in regions.items():
        print(f"  {name:14s}{len(terms):6,d} terms")
        sizes.append({"region": name, "n_terms": len(terms),
                      "roots": ";".join(cfg.HPO_SKELETAL_REGIONS[name])})

    overlap = sum(len(a & b) for i, a in enumerate(regions.values())
                  for b in list(regions.values())[i + 1:])
    print(f"  pairwise term overlap between regions: {overlap} "
          f"(the six roots are near-disjoint by construction)")

    terms_by_gene = common.gene_to_hpo_terms()
    print(f"\n  genes with any HPO disease annotation: {len(terms_by_gene):,}")

    rows = []
    for gene, terms in terms_by_gene.items():
        if not terms:
            continue
        hit = {name: int(bool(terms & members)) for name, members in regions.items()}
        rows.append({"gene": gene, "n_hpo_terms": len(terms),
                     "regions_skeletal": sum(hit.values()), **hit})
    d = pd.DataFrame(rows).drop_duplicates("gene").set_index("gene")
    d["has_skeletal_disease"] = (d.regions_skeletal > 0).astype(int)

    print(f"  genes with a phenotype term:           {len(d):,}")
    print(f"  genes with a skeletal region affected: "
          f"{int(d.has_skeletal_disease.sum()):,}")
    print("\n  The second number is the one the symmetry law runs on. The first")
    print("  is the one the term count must be computed over; using the second")
    print("  for both is the bug that once inverted the cross-organ control.")

    print("\n=== Distribution of affected regions, among skeletal disease genes ===")
    skeletal = d[d.has_skeletal_disease == 1]
    for k in range(1, 7):
        n = int((skeletal.regions_skeletal == k).sum())
        print(f"  {k} region(s): {n:5,d}  ({n / len(skeletal):5.1%})")
    print(f"  mean regions affected: {skeletal.regions_skeletal.mean():.2f}")
    print(f"  'generalised' (>= 4 regions): "
          f"{(skeletal.regions_skeletal >= 4).mean():.1%}")

    print("\n=== Ascertainment: how much of this is just annotation depth? ===")
    print(f"  median HPO terms, skeletal disease genes: "
          f"{skeletal.n_hpo_terms.median():.0f}")
    print(f"  median HPO terms, all disease genes:      {d.n_hpo_terms.median():.0f}")
    print("  The correlation with breadth is computed in "
          "04_analysis/02_symmetry_law.py.")

    common.write_result(d, args.out)
    common.write_result(pd.DataFrame(sizes), "hpo_region_sizes.tsv", index=False)
    print("\nNow run:  python 03_truth_sides/02_build_panelapp.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
