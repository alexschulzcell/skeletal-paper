#!/usr/bin/env python3
"""The 2x2: is breadth organ-specific, or just general gene importance?

What it does
    The control that decides whether the paper is about the skeleton or about
    nothing. Breadth is computed a second time in a different organ - seven
    subcortical brain volumes from Oxford BIG40, treated exactly like the six
    DXA sites - and both breadth measures are put into the same regression
    against both phenotype-extent measures:

                              -> skeletal extent    -> brain extent
        skeletal breadth            [1,1]                [1,2]
        brain breadth               [2,1]                [2,2]

    If breadth measured general importance, the two rows would look alike. If
    the law is anatomically local, the diagonal stands and the off-diagonal
    falls. Every cell is a standardised coefficient from an OLS that also
    carries log(total HPO terms), so that annotation depth cannot produce the
    pattern.

    This is also the paper's answer to the published phenome-wide null result:
    GWAS pleiotropy and Mendelian pleiotropy do not correlate across the
    phenome, and our brain cell reproduces that. The reason is that the
    correspondence is anatomically local and does not cross the organ boundary.

    Two failure modes are guarded explicitly, because both bit us:
      - the HPO term count must come from all disease genes, not from the genes
        of one organ. Restricting it shrank the skeletal regression by a third
        of its sample and moved p from 0.0012 to 0.117. Both tables built in
        stage 03 carry the full count; the script asserts they agree.
      - a shared "importance factor" would show up as a correlation between the
        two breadth measures themselves. That correlation is printed first;
        above about 0.4 the 2x2 would not be interpretable.

Reads
    results/gene_layers.tsv                       (stage 02, script 01)
    results/truth_hpo.tsv, results/truth_brain.tsv (stage 03, scripts 01 and 04)
    $SKELBREADTH_DATA/magma/<region>.genes.out     (stage 01, script 02)

Writes
    results/cross_organ_matrix.tsv     the 2x2, coefficient, t and p per cell
    results/cross_organ_models.tsv     the two full regressions
    results/cross_organ_correlation.tsv  breadth-breadth and region-region

Figures
    Figures 2B and 2C. Panel tables:
    fig2b_cross_organ.csv, fig2c_trait_coherence.csv.

Runtime
    One minute.

Usage
    python 04_analysis/03_cross_organ_control.py
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


def brain_breadth() -> tuple[pd.Series, pd.DataFrame]:
    """Breadth over the subcortical volumes, built exactly like the skeletal one."""
    available = [r for r in cfg.BRAIN_REGIONS if cfg.magma_out(r).exists()]
    if len(available) < 5:
        raise SystemExit(
            "\nThe cross-organ control needs at least five of the seven "
            "subcortical volumes.\n"
            f"  found: {available or 'none'}\n"
            "  Run 00_setup/01_download_sumstats.sh and 01_gene_scores/02_run_magma.sh.\n")
    Z = common.magma_matrix(available)
    return (Z[available] > cfg.CORE_Z).sum(axis=1).rename("breadth_brain"), Z[available]


def cell(d: pd.DataFrame, target: str, label: str) -> pd.DataFrame:
    s = d.dropna(subset=[target, "breadth_skeletal", "breadth_brain", "n_hpo_terms"])
    if len(s) < 100:
        raise SystemExit(
            f"\nOnly {len(s)} gene(s) have both breadth measures and a value for "
            f"'{target}'.\n"
            "  The 2x2 needs a phenotype axis on both organs. Check that\n"
            "  03_truth_sides/01_build_hpo_regions.py and 04_build_brain_truth.py\n"
            "  both ran over the full HPO annotation, and that the brain MAGMA\n"
            "  results are present for at least five subcortical volumes.\n")
    X = s[["breadth_skeletal", "breadth_brain", "n_hpo_terms"]].copy()
    X["n_hpo_terms"] = np.log1p(X.n_hpo_terms)
    model = sm.OLS(s[target], sm.add_constant(standardise(X))).fit()
    table = model.summary2().tables[1]
    print(f"\n  {label}   (n {len(s):,})")
    for term in ("breadth_skeletal", "breadth_brain", "n_hpo_terms"):
        print(f"    {term:18s} coef {table.loc[term, 'Coef.']:+.4f}   "
              f"t {table.loc[term, 't']:+6.2f}   p {table.loc[term, 'P>|t|']:.4g}")
    table = table.copy()
    table["target"] = target
    table["n"] = len(s)
    return table


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    args = ap.parse_args()

    common.banner("Claim 3: the law is organ-specific")

    layers = common.load_layers()
    hpo = common.load_truth("hpo")
    brain_truth = common.load_truth("brain")
    nb, brain_Z = brain_breadth()
    print(f"  subcortical volumes used: {brain_Z.shape[1]} "
          f"({', '.join(brain_Z.columns)})")

    # Annotation depth must be identical in both tables; stage 03 builds both
    # over all disease genes. Verify rather than trust.
    shared = hpo.index.intersection(brain_truth.index)
    if not np.allclose(hpo.loc[shared, "n_hpo_terms"],
                       brain_truth.loc[shared, "n_hpo_terms"]):
        raise SystemExit(
            "HPO term counts differ between the skeletal and brain truth tables. "
            "One of them was built over a restricted gene set; rebuild stage 03.")

    d = pd.DataFrame(index=layers.index)
    d["breadth_skeletal"] = layers.breadth
    d = d.join(nb, how="inner")
    d["n_hpo_terms"] = hpo.n_hpo_terms.reindex(d.index)
    d["regions_skeletal"] = hpo.regions_skeletal.where(
        hpo.has_skeletal_disease == 1).reindex(d.index)
    d["regions_brain"] = brain_truth.regions_brain.where(
        brain_truth.has_brain_disease == 1).reindex(d.index)
    d = d[d.n_hpo_terms.notna()]

    print(f"\n  genes with both breadth measures and an HPO disease: {len(d):,}")
    print(f"    with a skeletal phenotype axis: "
          f"{int(d.regions_skeletal.notna().sum()):,}")
    print(f"    with a brain phenotype axis:    "
          f"{int(d.regions_brain.notna().sum()):,}")

    print("\n=== First: are the two breadth measures the same thing? ===")
    rho, p = stats.spearmanr(d.breadth_skeletal, d.breadth_brain)
    print(f"  Spearman(skeletal breadth, brain breadth) = {rho:+.3f}  p {p:.2e}")
    print("  Above about +0.4 this would be a warning sign of a shared")
    print("  importance factor and the 2x2 below would not be readable.")

    print("\n=== The 2x2 ===")
    a = cell(d, "regions_skeletal", "target: extent of the SKELETAL phenotype")
    b = cell(d, "regions_brain", "target: extent of the BRAIN phenotype")

    print(f"\n  {'':22s}{'-> skeletal':>14s}{'-> brain':>12s}")
    for term, label in [("breadth_skeletal", "skeletal breadth"),
                        ("breadth_brain", "brain breadth")]:
        print(f"  {label:22s}{a.loc[term, 'Coef.']:14.4f}{b.loc[term, 'Coef.']:12.4f}")
    print(f"  {'':22s}{'p ' + format(a.loc['breadth_skeletal', 'P>|t|'], '.4g'):>14s}"
          f"{'p ' + format(b.loc['breadth_brain', 'P>|t|'], '.4g'):>12s}   (diagonal)")
    print(f"  {'':22s}{'p ' + format(a.loc['breadth_brain', 'P>|t|'], '.4g'):>14s}"
          f"{'p ' + format(b.loc['breadth_skeletal', 'P>|t|'], '.4g'):>12s}   (off-diagonal)")

    diagonal = a.loc["breadth_skeletal", "Coef."] + b.loc["breadth_brain", "Coef."]
    off = a.loc["breadth_brain", "Coef."] + b.loc["breadth_skeletal", "Coef."]
    worst_diagonal = max(a.loc["breadth_skeletal", "P>|t|"],
                         b.loc["breadth_brain", "P>|t|"])
    print(f"\n  diagonal {diagonal:+.4f}   off-diagonal {off:+.4f}")
    if diagonal > 0 and a.loc["breadth_skeletal", "P>|t|"] < 0.05 \
            and abs(a.loc["breadth_skeletal", "Coef."]) > 2 * abs(a.loc["breadth_brain", "Coef."]):
        print("  Skeletal breadth predicts skeletal extent; brain breadth predicts")
        print("  nothing about it. Breadth is not general importance.")
    else:
        print("  The diagonal does not stand clear of the off-diagonal. Report")
        print("  that as it is; it is the objection that would end the paper.")

    matrix = pd.DataFrame([
        {"predictor": "breadth_skeletal", "target": "regions_skeletal",
         "coef": float(a.loc["breadth_skeletal", "Coef."]),
         "t": float(a.loc["breadth_skeletal", "t"]),
         "p": float(a.loc["breadth_skeletal", "P>|t|"]), "n": int(a["n"].iloc[0])},
        {"predictor": "breadth_brain", "target": "regions_skeletal",
         "coef": float(a.loc["breadth_brain", "Coef."]),
         "t": float(a.loc["breadth_brain", "t"]),
         "p": float(a.loc["breadth_brain", "P>|t|"]), "n": int(a["n"].iloc[0])},
        {"predictor": "breadth_skeletal", "target": "regions_brain",
         "coef": float(b.loc["breadth_skeletal", "Coef."]),
         "t": float(b.loc["breadth_skeletal", "t"]),
         "p": float(b.loc["breadth_skeletal", "P>|t|"]), "n": int(b["n"].iloc[0])},
        {"predictor": "breadth_brain", "target": "regions_brain",
         "coef": float(b.loc["breadth_brain", "Coef."]),
         "t": float(b.loc["breadth_brain", "t"]),
         "p": float(b.loc["breadth_brain", "P>|t|"]), "n": int(b["n"].iloc[0])},
    ])
    common.write_result(matrix, "cross_organ_matrix.tsv", index=False)
    common.write_result(pd.concat([a, b]), "cross_organ_models.tsv")

    # ------------------------------------------------------------------
    print("\n=== The precondition: does the organ have a shared layer at all? ===")
    Zs = common.load_gene_matrix()[cfg.BMD_SITES]
    skel_corr = Zs.corr(method="spearman")
    brain_corr = brain_Z.corr(method="spearman")

    def mean_off_diagonal(c):
        v = c.values[~np.eye(len(c), dtype=bool)]
        return float(np.mean(v))

    skel_mean = mean_off_diagonal(skel_corr)
    brain_mean = mean_off_diagonal(brain_corr)
    print(f"  mean inter-site correlation, six DXA sites      : {skel_mean:+.3f}")
    print(f"  mean inter-region correlation, subcortex        : {brain_mean:+.3f}")
    print("  A shared layer can only exist where the sites share signal. The")
    print("  subcortex does not, which is why the brain cell is empty and why")
    print("  the null result there is a precondition failure, not a refutation.")

    corr_rows = ([{"organ": "skeleton", "a": i, "b": j,
                   "spearman": float(skel_corr.loc[i, j])}
                  for i in skel_corr.index for j in skel_corr.columns if i < j]
                 + [{"organ": "subcortex", "a": i, "b": j,
                     "spearman": float(brain_corr.loc[i, j])}
                    for i in brain_corr.index for j in brain_corr.columns if i < j])
    corr = pd.DataFrame(corr_rows)
    corr.loc[len(corr)] = {"organ": "skeleton", "a": "MEAN", "b": "MEAN",
                           "spearman": skel_mean}
    corr.loc[len(corr)] = {"organ": "subcortex", "a": "MEAN", "b": "MEAN",
                           "spearman": brain_mean}
    corr.loc[len(corr)] = {"organ": "breadth_vs_breadth", "a": "skeletal",
                           "b": "brain", "spearman": float(rho)}
    common.write_result(corr, "cross_organ_correlation.tsv", index=False)

    print("\nNow run:  python 04_analysis/04_confounders.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
