#!/usr/bin/env python3
"""What the shared layer is made of: Wnt and ossification, and its regulation.

What it does
    Three mechanistic readings of the core, in decreasing order of how much
    they depend on curation.

    1. Pathway enrichment. Eight signalling and developmental gene sets, taken
       from the Gene Ontology with full descendant propagation, tested in the
       core against background and reported next to the site-specific layer.
       The result is a Wnt/ossification axis, and the informative part is the
       contrast: the shared layer is the Wnt layer, the site-specific one is
       not. Figure 4A.

    2. Regulatory equipment, which is immune to publication bias because it is
       measured rather than curated. Open-chromatin peaks from human fetal limb
       ATAC are counted in a window around every gene, across eight skeletal
       elements. Peak count tracks gene length, so the comparison is made
       against a length-decile-matched null rather than against the raw
       background.

    3. The maturation axis, with external validation on it. A growth-plate
       pseudobulk matrix ordered resting -> proliferating -> prehypertrophic ->
       hypertrophic gives each gene a position on the axis, and the markers
       published independently in the skeletal-dysplasia literature are placed
       on that axis: GAS1, SFRP5 and PTHLH at the resting pole, SGMS2 at the
       maturation pole. Figure 4B. This part needs the pseudobulk matrix, which
       is not a public download; when it is absent the script says so and skips
       that section, and the first two sections still run.

Reads
    results/gene_layers.tsv                              (stage 02, script 01)
    $SKELBREADTH_DATA/reference/goa_human.gaf.gz, go-basic.obo
    $SKELBREADTH_DATA/reference/peaks/*hg19.bed.gz       (optional)
    $SKELBREADTH_DATA/reference/growth_plate_pseudobulk.tsv  (optional)

Writes
    results/pathway_enrichment.tsv     eight gene sets x three layers
    results/pathway_core_members.tsv   which core genes are in which pathway
    results/regulatory_density.tsv     peaks per element per layer, with the
                                       length-matched null
    results/maturation_axis.tsv        per-gene axis position (when available)
    results/maturation_axis_markers.tsv  the published markers on that axis

Numbers in the paper
    This script feeds Figure 4. The panel tables it ends up in are
    fig4a_pathways.csv, fig4b_regulatory_density.csv,
    fig4c_shared_layer_genes.csv.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Two to five minutes, most of it counting ATAC peaks.

Usage
    python 04_analysis/07_pathways_and_axis.py
    python 04_analysis/07_pathways_and_axis.py --skip-atac
"""

from __future__ import annotations

import argparse
import collections
import glob
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402

#: Gene Ontology roots for the eight sets, with descendants taken.
PATHWAYS = {
    "FGF":                  ["GO:0008543"],
    "Ossification":         ["GO:0001503"],
    "Wnt, canonical":       ["GO:0060070"],
    "Wnt, all":             ["GO:0016055"],
    "Cartilage development": ["GO:0051216"],
    "BMP":                  ["GO:0030509"],
    "Hedgehog":             ["GO:0007224"],
    "TGF-beta":             ["GO:0007179"],
}

#: Growth-plate zone order, resting to hypertrophic.
ZONE_ORDER = ["RestingChon", "ProlifChon", "PrehyperChon", "HyperChon"]

#: Markers published independently of this work, placed on our axis as external
#: validation. Sign is the expected pole: negative = resting, positive = mature.
EXTERNAL_MARKERS = {"GAS1": "resting", "SFRP5": "resting", "PTHLH": "resting",
                    "SGMS2": "mature"}

ATAC_WINDOW = 100_000


def go_descendants(roots: dict) -> dict:
    """Full descendant set of each GO root, following is_a and part_of."""
    path = cfg.require(cfg.REFERENCE / "go-basic.obo", "Gene Ontology (go-basic.obo)")
    children = collections.defaultdict(list)
    current = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current = None
            elif line.startswith("id: GO:"):
                current = line[4:]
            elif current and line.startswith("is_a:"):
                children[line[6:16]].append(current)
            elif current and line.startswith("relationship: part_of GO:"):
                children[line[23:33]].append(current)

    def descend(root):
        seen, stack = set(), [root]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(children.get(node, []))
        return seen

    return {name: set().union(*[descend(r) for r in terms])
            for name, terms in roots.items()}


def go_gene_sets(term_sets: dict) -> dict:
    """Map each pathway to its gene symbols, from the human GAF."""
    path = cfg.require(cfg.REFERENCE / "goa_human.gaf.gz", "GO annotations (GAF)")
    by_term = collections.defaultdict(set)
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) < 10 or f[6] in {"IEA", "ND"}:
                continue      # curated evidence only; IEA is machine-inferred
            by_term[f[4]].add(f[2])
    return {name: set().union(*[by_term.get(t, set()) for t in terms])
            for name, terms in term_sets.items()}


def count_peaks(d: pd.DataFrame, window: int) -> pd.DataFrame | None:
    """Peaks per skeletal element within `window` bp of each gene."""
    files = sorted(glob.glob(str(cfg.REFERENCE / "peaks" / "*hg19.bed.gz")))
    if not files:
        return None
    elements = [os.path.basename(f).split("_")[-3] if "_" in f else os.path.basename(f)
                for f in files]
    print(f"  {len(files)} peak file(s): {', '.join(elements)}")
    counts = np.zeros((len(d), len(files)), dtype=int)
    chroms = d["chr"].astype(str).values
    lo = d.gene_start.values - window
    hi = d.gene_end.values + window
    for j, f in enumerate(files):
        by_chrom: dict[str, list] = {}
        with gzip.open(f, "rt", errors="ignore") as fh:
            for line in fh:
                p = line.split()
                if len(p) < 3:
                    continue
                by_chrom.setdefault(p[0].replace("chr", ""), []).append(
                    (int(p[1]) + int(p[2])) // 2)
        sorted_centres = {c: np.array(sorted(v)) for c, v in by_chrom.items()}
        for i in range(len(d)):
            arr = sorted_centres.get(chroms[i])
            if arr is None:
                continue
            counts[i, j] = int(np.searchsorted(arr, hi[i])
                               - np.searchsorted(arr, lo[i]))
    return pd.DataFrame(counts, index=d.index, columns=elements)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--skip-atac", action="store_true")
    ap.add_argument("--window", type=int, default=ATAC_WINDOW)
    ap.add_argument("--draws", type=int, default=cfg.N_PERMUTATIONS)
    args = ap.parse_args()

    d = common.load_layers()
    core = set(d.index[d.layer == "core"])
    specific = set(d.index[d.layer == "specific"])
    rest = set(d.index[d.layer == "rest"])

    common.banner("Mechanism: what the shared layer is made of")
    print(f"  core {len(core)} | site-specific {len(specific)} | rest {len(rest):,}")

    # ---------------------------------------------------------------- 1
    print("\n=== Pathway enrichment (Figure 4A) ===")
    gene_sets = go_gene_sets(go_descendants(PATHWAYS))
    print(f"{'pathway':22s}{'core':>10s}{'specific':>11s}{'background':>12s}"
          f"{'OR':>8s}{'p':>11s}")
    rows, members = [], []
    for name, genes in gene_sets.items():
        S = genes & set(d.index)
        if len(S) < 20:
            print(f"{name:22s}   too few annotated genes ({len(S)})")
            continue
        a, b, c = len(core & S), len(specific & S), len(rest & S)
        orr, p = stats.fisher_exact([[a, len(core) - a], [c, len(rest) - c]])
        print(f"{name:22s}{f'{a}/{len(core)}':>10s}{f'{b}/{len(specific)}':>11s}"
              f"{c / len(rest):12.1%}{orr:8.2f}{p:11.2e}")
        rows.append({"pathway": name, "n_annotated": len(S),
                     "core_hits": a, "core_n": len(core),
                     "specific_hits": b, "specific_n": len(specific),
                     "background_share": round(c / len(rest), 5),
                     "odds_ratio": orr, "p": p})
        members += [{"pathway": name, "gene": g, "layer": "core"}
                    for g in sorted(core & S)]
        members += [{"pathway": name, "gene": g, "layer": "specific"}
                    for g in sorted(specific & S)]
    common.write_result(pd.DataFrame(rows), "pathway_enrichment.tsv", index=False)
    common.write_result(pd.DataFrame(members), "pathway_core_members.tsv",
                        index=False)

    wnt = gene_sets.get("Wnt, all", set())
    print(f"\n  Wnt genes in the core           : {', '.join(sorted(core & wnt))}")
    print(f"  Wnt genes in the specific layer : "
          f"{', '.join(sorted(specific & wnt)) or '(none)'}")
    print("\n  Wnt as a bone-density axis is long established. What is new is")
    print("  that the *shared* layer coincides with it and the site-specific")
    print("  one does not.")

    # ---------------------------------------------------------------- 2
    if args.skip_atac:
        print("\n=== Regulatory equipment: skipped (--skip-atac) ===")
    else:
        print(f"\n=== Regulatory equipment, ATAC peaks within "
              f"+/- {args.window // 1000} kb ===")
        counts = count_peaks(d, args.window)
        if counts is None:
            print("  No ATAC peak files found under reference/peaks/.")
            print("  Run 00_setup/02_download_references.sh to include them;")
            print("  the rest of this script does not depend on them.")
        else:
            d["peaks_per_element"] = counts.mean(axis=1)
            d["regulatory_breadth"] = (counts > 0).sum(axis=1)
            d["gene_length"] = d.gene_end - d.gene_start

            print(f"\n  {'layer':16s}{'n':>6s}{'peaks/element':>16s}"
                  f"{'elements with a peak':>22s}{'gene length kb':>16s}")
            for name, S in [("core", core), ("site-specific", specific),
                            ("rest", rest)]:
                g = d.loc[sorted(S)]
                print(f"  {name:16s}{len(g):6,d}{g.peaks_per_element.median():16.1f}"
                      f"{g.regulatory_breadth.median():22.1f}"
                      f"{g.gene_length.median() / 1000:16.0f}")

            rho, _ = stats.spearmanr(d.gene_length, d.peaks_per_element)
            print(f"\n  The confounder, stated: Spearman(gene length, peaks) "
                  f"= {rho:+.3f}")
            print("  so the comparison below is against a length-matched null,")
            print("  not against the raw background.")

            d["length_decile"] = pd.qcut(np.log10(d.gene_length + 1), 10,
                                         labels=False, duplicates="drop")
            rng = np.random.default_rng(cfg.SEED)
            reg_rows = []
            for name, S in [("core", core), ("site-specific", specific)]:
                sub = d.loc[sorted(S)]
                observed = float(sub.peaks_per_element.median())
                null = []
                pools = {b: g.index.values for b, g in
                         d[~d.index.isin(S)].groupby("length_decile")}
                for _ in range(args.draws):
                    picked = []
                    for b, g in sub.groupby("length_decile"):
                        pool = pools.get(b)
                        if pool is None or len(pool) < len(g):
                            continue
                        picked.extend(rng.choice(pool, len(g), replace=False))
                    if picked:
                        null.append(float(d.loc[picked, "peaks_per_element"].median()))
                null = np.asarray(null)
                z = (observed - null.mean()) / (null.std() or np.nan)
                p = (1 + (null >= observed).sum()) / (len(null) + 1)
                print(f"  {name:16s} peaks/element {observed:5.1f}   null "
                      f"{null.mean():5.1f} +- {null.std():.2f}   "
                      f"z {z:+.2f}   p {p:.4f}")
                reg_rows.append({"layer": name, "n": len(sub),
                                 "observed_peaks_per_element": round(observed, 3),
                                 "null_mean": round(float(null.mean()), 3),
                                 "null_sd": round(float(null.std()), 3),
                                 "z": round(float(z), 3), "p_empirical": float(p),
                                 "draws": len(null)})

            X = d[["breadth", "peak"]].copy()
            X["log_length"] = np.log10(d.gene_length + 1)
            X = (X - X.mean()) / X.std()
            model = sm.OLS(np.log1p(d.peaks_per_element), sm.add_constant(X)).fit()
            table = model.summary2().tables[1]
            print("\n  Continuous, with length in the model:")
            for term in ("breadth", "peak", "log_length"):
                print(f"    {term:12s} coef {table.loc[term, 'Coef.']:+.4f}   "
                      f"t {table.loc[term, 't']:+6.2f}   "
                      f"p {table.loc[term, 'P>|t|']:.3g}")
                reg_rows.append({"layer": f"model:{term}", "n": len(d),
                                 "observed_peaks_per_element": np.nan,
                                 "null_mean": np.nan, "null_sd": np.nan,
                                 "z": float(table.loc[term, "t"]),
                                 "p_empirical": float(table.loc[term, "P>|t|"]),
                                 "draws": 0})
            common.write_result(pd.DataFrame(reg_rows), "regulatory_density.tsv",
                                index=False)

    # ---------------------------------------------------------------- 3
    print("\n=== Maturation axis, with the published markers on it (Figure 4B) ===")
    pseudobulk = cfg.REFERENCE / "growth_plate_pseudobulk.tsv"
    if not pseudobulk.exists():
        print(f"  Not available: {pseudobulk.name}")
        print("  This is the human growth-plate pseudobulk matrix (zones x genes).")
        print("  It is not a public download; see docs/DATA_SOURCES.md for how it")
        print("  is derived. Sections 1 and 2 above do not depend on it.")
        print("\nStage 04 complete. The figures are built from results/ by the")
        print("scripts in 05_figures/.")
        return 0

    M = pd.read_csv(pseudobulk, sep="\t", index_col=0)
    zones = [z for z in ZONE_ORDER if z in M.index]
    if len(zones) < 3:
        raise SystemExit(f"Pseudobulk matrix needs the zones {ZONE_ORDER}; "
                         f"found {list(M.index)}")
    cpm = M.loc[zones].div(M.loc[zones].sum(axis=1), axis=0) * 1e6
    expressed = cpm.max(axis=0) > 5
    lg = np.log2(cpm.loc[:, expressed] + 1)

    rank = np.arange(len(zones), dtype=float)
    centred = lg - lg.mean(axis=0)
    r = ((centred.T @ (rank - rank.mean()))
         / (np.sqrt((centred ** 2).sum(axis=0))
            * np.sqrt(((rank - rank.mean()) ** 2).sum())))
    fold = lg.loc[zones[-1]] - lg.loc[zones[0]]
    axis = pd.DataFrame({"axis_r": r, "log2_fold_resting_to_mature": fold}).dropna()
    axis["loading"] = axis.axis_r * axis.log2_fold_resting_to_mature.abs()
    axis["pole"] = np.where(axis.axis_r > 0, "mature", "resting")
    axis["in_core"] = axis.index.isin(core).astype(int)

    print(f"  genes on the axis: {len(axis):,}   zones: {', '.join(zones)}")
    print(f"  core genes on the axis: {int(axis.in_core.sum())}")

    print(f"\n  {'marker':10s}{'loading':>10s}{'pole':>10s}{'expected':>11s}"
          f"{'agrees':>9s}")
    marker_rows = []
    for gene, expected in EXTERNAL_MARKERS.items():
        if gene not in axis.index:
            print(f"  {gene:10s}      not expressed above threshold")
            continue
        row = axis.loc[gene]
        agrees = row.pole == expected
        print(f"  {gene:10s}{row.loading:+10.2f}{row.pole:>10s}{expected:>11s}"
              f"{'yes' if agrees else 'no':>9s}")
        marker_rows.append({"gene": gene, "loading": round(float(row.loading), 4),
                            "axis_r": round(float(row.axis_r), 4),
                            "pole": row.pole, "expected_pole": expected,
                            "agrees": bool(agrees)})
    if marker_rows:
        agree = sum(m["agrees"] for m in marker_rows)
        print(f"\n  {agree} of {len(marker_rows)} independently published markers")
        print("  land on the pole they were reported at. The axis was built with")
        print("  no knowledge of them.")
    common.write_result(axis, "maturation_axis.tsv")
    common.write_result(pd.DataFrame(marker_rows), "maturation_axis_markers.tsv",
                        index=False)

    print("\nStage 04 complete. The figures are built from results/ by the")
    print("scripts in 05_figures/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
