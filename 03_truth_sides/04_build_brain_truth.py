#!/usr/bin/env python3
"""Build the brain disease side, for the cross-organ control only.

What it does
    The same construction as 01_build_hpo_regions.py, applied to the brain:
    eight neuroanatomical HPO roots - cortex, basal ganglia, cerebellum,
    brainstem, corpus callosum, ventricles, white matter, hippocampus - and,
    per gene, how many of them its Mendelian phenotype touches.

    This exists for one purpose. The most dangerous objection to the paper is
    that breadth measures nothing but general gene importance, in which case
    breadth measured in the brain should predict skeletal phenotype extent just
    as well as breadth measured in the skeleton. Answering that needs a brain
    phenotype axis built exactly like the skeletal one, on the same ontology,
    with the same term-count covariate. That is what this table is. It carries
    no claim of its own, and the paper reports no brain result other than the
    control.

    The term count is again taken over *all* genes with an HPO disease, not
    only over brain-disease genes, for the reason set out in
    01_build_hpo_regions.py.

Reads
    $SKELBREADTH_DATA/reference/hp.obo
    $SKELBREADTH_DATA/reference/phenotype.hpoa
    $SKELBREADTH_DATA/reference/genes_to_disease.txt

Writes
    results/truth_brain.tsv       gene, n_hpo_terms, regions_brain,
                                  has_brain_disease, one column per region

Runtime
    One to two minutes.

Usage
    python 03_truth_sides/04_build_brain_truth.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="truth_brain.tsv")
    args = ap.parse_args()

    common.banner("Truth side 4 of 4: brain regions, for the cross-organ control")

    regions = common.hpo_descendants(cfg.HPO_BRAIN_REGIONS)
    print("\n=== Eight neuroanatomical regions ===")
    for name, terms in regions.items():
        print(f"  {name:18s}{len(terms):6,d} terms")

    terms_by_gene = common.gene_to_hpo_terms()
    rows = []
    for gene, terms in terms_by_gene.items():
        if not terms:
            continue
        hit = {name: int(bool(terms & members)) for name, members in regions.items()}
        rows.append({"gene": gene, "n_hpo_terms": len(terms),
                     "regions_brain": sum(hit.values()), **hit})
    d = pd.DataFrame(rows).drop_duplicates("gene").set_index("gene")
    d["has_brain_disease"] = (d.regions_brain > 0).astype(int)

    print(f"\n  genes with a phenotype term:        {len(d):,}")
    print(f"  genes with a brain region affected: "
          f"{int(d.has_brain_disease.sum()):,}")

    try:
        skeletal = common.load_truth("hpo")
        both = d.index.intersection(
            skeletal[skeletal.has_skeletal_disease == 1].index)
        both = [g for g in both if d.loc[g, "has_brain_disease"] == 1]
        print(f"  genes with both organ axes:         {len(both):,}")
        print("\n  That last number is the sample the 2x2 control runs on. It is")
        print("  smaller than either organ alone, and the control is powered by")
        print("  the full regression, not by this intersection - see")
        print("  04_analysis/03_cross_organ_control.py.")
    except SystemExit:
        print("  (run 03_truth_sides/01_build_hpo_regions.py to see the overlap)")

    print("\n=== Distribution of affected brain regions ===")
    brain = d[d.has_brain_disease == 1]
    for k in range(1, 9):
        n = int((brain.regions_brain == k).sum())
        if n:
            print(f"  {k} region(s): {n:5,d}  ({n / len(brain):5.1%})")
    print(f"  mean regions affected: {brain.regions_brain.mean():.2f}")

    common.write_result(d, args.out)
    print("\nNow run:  python 04_analysis/01_enrichment.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
