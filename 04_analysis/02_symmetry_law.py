#!/usr/bin/env python3
"""Acts widely, breaks widely: breadth predicts how far the phenotype extends.

What it does
    The second and central claim. Among genes that carry a Mendelian skeletal
    disease, it tests whether breadth on the common-variant side predicts the
    number of body regions the Mendelian phenotype touches - the same six
    anatomical regions on both sides of the comparison.

    Three things are computed:

    1. The law itself. Spearman correlation and the mean number of affected
       regions at each breadth, plus the extreme contrast (breadth <= 1 against
       breadth >= 5) by Mann-Whitney.
    2. The ascertainment control, which is what makes the result readable.
       A better-studied gene has more HPO terms. If the target were a count of
       terms, that alone would produce the law. It is not: the target is a
       proportion, regions out of six, and the correlation between breadth and
       total term count is essentially zero. Both the raw correlation and an
       OLS with log(term count) as a covariate are reported.
    3. The clinical reading. What share of genes at each breadth produce a
       "generalised" phenotype, four or more regions of six.

Reads
    results/gene_layers.tsv       (stage 02, script 01)
    results/truth_hpo.tsv         (stage 03, script 01)

Writes
    results/symmetry_law.tsv           per-gene table the law runs on
    results/symmetry_by_breadth.tsv    mean regions and generalised share, 0-6
    results/symmetry_models.tsv        correlations and the OLS coefficients

Numbers in the paper
    This script feeds Figure 2A. The panel tables it ends up in are
    fig2a_symmetry_law.csv, fig2a_symmetry_stats.csv,
    fig2a_symmetry_regression.csv.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Seconds.

Usage
    python 04_analysis/02_symmetry_law.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def standardise(X: pd.DataFrame) -> pd.DataFrame:
    return (X - X.mean()) / X.std()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--generalised-threshold", type=int, default=4,
                    help="regions of six that count as generalised (default 4)")
    args = ap.parse_args()

    layers = common.load_layers()
    hpo = common.load_truth("hpo")

    d = layers.join(hpo[["n_hpo_terms", "regions_skeletal",
                         "has_skeletal_disease"]], how="inner")
    S = d[(d.n_hpo_terms > 0) & (d.has_skeletal_disease == 1)].copy()

    common.banner("Claim 2: the symmetry law")
    print(f"  genes with a breadth value and an HPO disease: {len(d):,}")
    print(f"  of those, with a Mendelian skeletal disease:   {len(S):,}")
    print("  The law is asked only of genes that have a skeletal disease at all;")
    print("  it is about how far the phenotype reaches, not whether there is one.")

    # ---------------------------------------------------------------- 1
    print("\n=== The law ===")
    rho, p_rho = stats.spearmanr(S.breadth, S.regions_skeletal)
    print(f"  Spearman(breadth, affected body regions) = {rho:+.3f}  "
          f"p {p_rho:.2e}  n {len(S):,}")

    print(f"\n{'breadth':>8s}{'genes':>8s}{'mean regions':>15s}"
          f"{'generalised':>14s}")
    by_breadth = []
    for k in range(7):
        s = S[S.breadth == k]
        if len(s) < 15:
            continue
        gen = float((s.regions_skeletal >= args.generalised_threshold).mean())
        print(f"{k:8d}{len(s):8,d}{s.regions_skeletal.mean():15.2f}{gen:14.1%}")
        by_breadth.append({"breadth": k, "n_genes": len(s),
                           "mean_regions": round(float(s.regions_skeletal.mean()), 4),
                           "generalised_share": round(gen, 4),
                           "median_hpo_terms": float(s.n_hpo_terms.median())})
    common.write_result(pd.DataFrame(by_breadth), "symmetry_by_breadth.tsv",
                        index=False)

    low, high = S[S.breadth <= 1], S[S.breadth >= 5]
    u, p_mw = stats.mannwhitneyu(high.regions_skeletal, low.regions_skeletal)
    print(f"\n  breadth <= 1 : {low.regions_skeletal.mean():.2f} regions  "
          f"(n {len(low):,})")
    print(f"  breadth >= 5 : {high.regions_skeletal.mean():.2f} regions  "
          f"(n {len(high):,})")
    print(f"  Mann-Whitney p {p_mw:.2e}")
    print(f"  generalised share: {(low.regions_skeletal >= args.generalised_threshold).mean():.0%}"
          f"  ->  {(high.regions_skeletal >= args.generalised_threshold).mean():.0%}")

    # ---------------------------------------------------------------- 2
    print("\n=== Ascertainment control: is this just research intensity? ===")
    rho_asc, p_asc = stats.spearmanr(S.breadth, S.n_hpo_terms)
    print(f"  Spearman(breadth, total HPO terms) = {rho_asc:+.3f}  p {p_asc:.3f}")
    print("  Near zero, and that is structural rather than lucky: the target is")
    print("  a proportion - regions out of six - so research effort that adds")
    print("  terms within a region cannot move it.")

    X = S[["breadth", "n_hpo_terms"]].copy()
    X["n_hpo_terms"] = np.log1p(X.n_hpo_terms)
    model = sm.OLS(S.regions_skeletal, sm.add_constant(standardise(X))).fit()
    table = model.summary2().tables[1]
    print(f"\n  OLS: affected regions ~ breadth + log(HPO terms)   n {len(S):,}")
    for term in ("breadth", "n_hpo_terms"):
        print(f"    {term:14s} coef {table.loc[term, 'Coef.']:+.4f}   "
              f"t {table.loc[term, 't']:+6.2f}   p {table.loc[term, 'P>|t|']:.4g}")
    print(f"    R2 {model.rsquared:.4f}")

    # A second form of the same control: the proportion as the target.
    S["region_share"] = S.regions_skeletal / 6.0
    model_share = sm.OLS(S.region_share, sm.add_constant(standardise(X))).fit()
    t_share = model_share.summary2().tables[1]

    models = [
        {"model": "spearman", "term": "breadth", "estimate": rho, "statistic": np.nan,
         "p": p_rho, "n": len(S)},
        {"model": "spearman_ascertainment", "term": "n_hpo_terms",
         "estimate": rho_asc, "statistic": np.nan, "p": p_asc, "n": len(S)},
        {"model": "mannwhitney_extremes", "term": "breadth>=5 vs <=1",
         "estimate": float(high.regions_skeletal.mean() - low.regions_skeletal.mean()),
         "statistic": u, "p": p_mw, "n": len(low) + len(high)},
    ]
    for term in ("breadth", "n_hpo_terms"):
        models.append({"model": "ols_regions", "term": term,
                       "estimate": float(table.loc[term, "Coef."]),
                       "statistic": float(table.loc[term, "t"]),
                       "p": float(table.loc[term, "P>|t|"]), "n": len(S)})
        models.append({"model": "ols_region_share", "term": term,
                       "estimate": float(t_share.loc[term, "Coef."]),
                       "statistic": float(t_share.loc[term, "t"]),
                       "p": float(t_share.loc[term, "P>|t|"]), "n": len(S)})
    common.write_result(pd.DataFrame(models), "symmetry_models.tsv", index=False)

    # ---------------------------------------------------------------- 3
    print("\n=== What it would mean at a first presentation ===")
    base = float((S.regions_skeletal >= args.generalised_threshold).mean())
    print(f"  base rate of a generalised phenotype: {base:.1%}")
    print(f"\n{'rule':34s}{'genes':>8s}{'generalised':>14s}{'factor':>9s}")
    for label, mask in [("breadth 0 - no shared signal", S.breadth == 0),
                        ("breadth >= 3", S.breadth >= 3),
                        ("breadth >= 5", S.breadth >= 5),
                        ("breadth 6 - all six sites", S.breadth == 6)]:
        s = S[mask]
        if not len(s):
            continue
        share = float((s.regions_skeletal >= args.generalised_threshold).mean())
        print(f"{label:34s}{len(s):8,d}{share:14.1%}{share / base:9.2f}x")

    common.write_result(
        S[["breadth", "peak", "layer", "n_hpo_terms", "regions_skeletal",
           "region_share"]], "symmetry_law.tsv")
    print("\nNow run:  python 04_analysis/03_cross_organ_control.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
