#!/usr/bin/env python3
"""The four confounders that could produce breadth without any biology.

What it does
    Works through the alternative explanations in order of danger, each with a
    matched null rather than a covariate alone, because a covariate assumes the
    functional form and a matched null does not.

    1. Gene length and variant count. A long gene collects more variants and
       can reach Z > 2 at more sites for that reason alone. Tested three ways:
       the raw correlation between variant count and breadth, the separation
       within variant-count deciles, and a null in which the core is replaced
       2,000 times by genes matched on variant-count decile.
    2. Research intensity. Better-studied genes are better annotated and more
       often on panels. Tested with publication counts per gene from
       gene2pubmed, as a covariate in a logistic model and as a matched null.
    3. Constraint. Constrained genes are known disease genes, so if the core
       were simply constrained the two-layer result would be empty. Tested as
       a layer comparison of LOEUF and as a LOEUF-decile-matched null.
    4. General pleiotropy, which is handled in 03_cross_organ_control.py and
       only referenced here so the set of six is complete alongside linkage
       disequilibrium (02_layers/02) and the ClinVar benign control
       (04_analysis/01).

Reads
    results/gene_layers.tsv, results/gene_z_matrix.tsv       (stages 01-02)
    results/truth_panelapp.tsv, results/truth_hpo.tsv        (stage 03)
    $SKELBREADTH_DATA/magma/head.genes.out    (for NSNPS per gene)
    $SKELBREADTH_DATA/reference/gnomad_constraint.txt.bgz
    $SKELBREADTH_DATA/reference/gene2pubmed.gz    (optional; skipped if absent)

Writes
    results/confounder_summary.tsv    one row per confounder, with the verdict
    results/confounder_matched_nulls.tsv  observed, null mean, sd, z, empirical p
    results/confounder_models.tsv     the logistic models with covariates

Figures
    None; the confounder controls are reported in the text.

Runtime
    Two to four minutes; the matched nulls dominate.

Usage
    python 04_analysis/04_confounders.py
    python 04_analysis/04_confounders.py --draws 500     # faster, coarser
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def matched_null(d: pd.DataFrame, members: pd.Index, target: str,
                 bin_column: str, draws: int, seed: int) -> dict:
    """Replace `members` by genes matched on `bin_column`, `draws` times.

    Returns the observed rate, the null distribution's mean and sd, the z score
    and an empirical p. Matching is stratified: each draw takes, from every bin,
    as many non-member genes as the member set has in that bin.
    """
    member_set = set(members)
    sub = d.loc[list(members)]
    observed = float(sub[target].mean())
    rng = np.random.default_rng(seed)
    pool_by_bin = {b: g.index.values
                   for b, g in d[~d.index.isin(member_set)].groupby(bin_column)}
    null = []
    for _ in range(draws):
        picked = []
        for b, group in sub.groupby(bin_column):
            pool = pool_by_bin.get(b)
            if pool is None or len(pool) < len(group):
                continue
            picked.extend(rng.choice(pool, len(group), replace=False))
        if picked:
            null.append(float(d.loc[picked, target].mean()))
    null = np.asarray(null)
    sd = null.std() or np.nan
    return {"observed": observed, "null_mean": float(null.mean()),
            "null_sd": float(sd), "z": float((observed - null.mean()) / sd),
            "p_empirical": float((1 + (null >= observed).sum()) / (len(null) + 1)),
            "draws": len(null)}


def load_publication_counts(index) -> pd.Series | None:
    """log10(1 + publications) per gene symbol, from NCBI gene2pubmed."""
    path = cfg.REFERENCE / "gene2pubmed.gz"
    if not path.exists():
        return None
    counts: dict[int, int] = {}
    with gzip.open(path, "rt", errors="ignore") as fh:
        fh.readline()
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] != "9606":
                continue
            counts[int(parts[1])] = counts.get(int(parts[1]), 0) + 1
    loc = common.gene_location_table()
    by_symbol = pd.Series(counts).rename("n")
    by_symbol.index = by_symbol.index.map(loc["symbol"])
    by_symbol = by_symbol[by_symbol.index.notna()].groupby(level=0).max()
    return np.log10(by_symbol.reindex(index) + 1)


def load_loeuf(index) -> pd.Series | None:
    path = cfg.REFERENCE / "gnomad_constraint.txt.bgz"
    if not path.exists():
        return None
    con = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False,
                      usecols=["gene", "oe_lof_upper"])
    con = con.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
    return con["oe_lof_upper"].reindex(index)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--draws", type=int, default=cfg.N_PERMUTATIONS)
    args = ap.parse_args()

    d = common.load_layers()
    panel = common.load_truth("panelapp")
    hpo = common.load_truth("hpo")
    d["on_panel"] = d.index.isin(panel.index).astype(int)
    d["n_hpo_terms"] = hpo.n_hpo_terms.reindex(d.index)
    d["regions_skeletal"] = hpo.regions_skeletal.where(
        hpo.has_skeletal_disease == 1).reindex(d.index)

    # Variant count and gene length, as MAGMA saw them.
    head = pd.read_csv(cfg.require(cfg.magma_out("head"),
                                   "MAGMA output for the head site"),
                       sep=r"\s+", index_col=0)
    loc = common.gene_location_table()
    head["symbol"] = loc["symbol"].reindex(head.index).values
    head = head.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    d["n_snps"] = head["NSNPS"].reindex(d.index)
    d["gene_length"] = (d.gene_end - d.gene_start)

    common.banner("The confounders")
    print(f"  genes {len(d):,} | core {int((d.layer == 'core').sum())} "
          f"| panel genes {int(d.on_panel.sum())}")

    core = d.index[d.layer == "core"]
    summary, nulls, models = [], [], []

    # ---------------------------------------------------------------- 1
    print("\n=== 1. Gene length and variant count ===")
    s = d.dropna(subset=["n_snps"]).copy()
    for column, label in [("n_snps", "variant count"), ("gene_length", "gene length")]:
        rho, p = stats.spearmanr(s[column], s.breadth)
        rho_min, _ = stats.spearmanr(s[column], s.min_z)
        print(f"  Spearman({label:14s}, breadth) = {rho:+.3f}  p {p:.1e}"
              f"   (against the weakest-site Z: {rho_min:+.3f})")
        summary.append({"confounder": "length", "test": f"rho({column}, breadth)",
                        "estimate": rho, "p": p})
    s["snp_decile"] = pd.qcut(s.n_snps, 10, labels=False, duplicates="drop")
    result = matched_null(s, s.index.intersection(core), "on_panel",
                          "snp_decile", args.draws, cfg.SEED)
    print(f"\n  core panel share {result['observed']:.1%}   "
          f"variant-count-matched null {result['null_mean']:.1%} "
          f"+- {result['null_sd']:.1%}")
    print(f"  z {result['z']:+.2f}   empirical p {result['p_empirical']:.4f}  "
          f"({result['draws']} draws)")
    nulls.append({"confounder": "variant count", "matched_on": "variant-count decile",
                  **result})
    summary.append({"confounder": "length", "test": "variant-count-matched null",
                    "estimate": result["z"], "p": result["p_empirical"]})

    # ---------------------------------------------------------------- 2
    print("\n=== 2. Research intensity ===")
    rho, p = stats.spearmanr(d.breadth, d.n_hpo_terms, nan_policy="omit")
    print(f"  Spearman(breadth, total HPO terms) = {rho:+.3f}  p {p:.3f}")
    summary.append({"confounder": "research intensity",
                    "test": "rho(breadth, HPO terms)", "estimate": rho, "p": p})

    pubs = load_publication_counts(d.index)
    if pubs is None:
        print("  gene2pubmed not present - publication controls skipped.")
        print("  Run 00_setup/02_download_references.sh to include them.")
    else:
        d["log_pubs"] = pubs
        s = d.dropna(subset=["log_pubs"]).copy()
        print(f"  publication counts mapped for {len(s):,} genes")
        print(f"  median publications: panel genes "
              f"{10 ** s.loc[s.on_panel == 1, 'log_pubs'].median() - 1:.0f}  "
              f"against {10 ** s.loc[s.on_panel == 0, 'log_pubs'].median() - 1:.0f}")
        X = pd.DataFrame({
            "core": s.index.isin(core).astype(int),
            "log_pubs": (s.log_pubs - s.log_pubs.mean()) / s.log_pubs.std(),
            "log_length": np.log10(s.gene_length + 1),
        }, index=s.index)
        X["log_length"] = (X.log_length - X.log_length.mean()) / X.log_length.std()
        model = sm.Logit(s.on_panel, sm.add_constant(X)).fit(disp=0)
        table = model.summary2().tables[1]
        print(f"\n  logit(panel gene) ~ core + log(publications) + log(length)")
        for term in ("core", "log_pubs", "log_length"):
            print(f"    {term:12s} coef {table.loc[term, 'Coef.']:+.3f}   "
                  f"z {table.loc[term, 'z']:+6.2f}   p {table.loc[term, 'P>|z|']:.4g}")
            models.append({"model": "panel ~ core + pubs + length", "term": term,
                           "coef": float(table.loc[term, "Coef."]),
                           "z": float(table.loc[term, "z"]),
                           "p": float(table.loc[term, "P>|z|"]), "n": len(s)})
        summary.append({"confounder": "research intensity",
                        "test": "core coefficient with log(publications)",
                        "estimate": float(table.loc["core", "z"]),
                        "p": float(table.loc["core", "P>|z|"])})

        s["pub_decile"] = pd.qcut(s.log_pubs, 10, labels=False, duplicates="drop")
        result = matched_null(s, s.index.intersection(core), "on_panel",
                              "pub_decile", args.draws, cfg.SEED + 1)
        print(f"\n  core panel share {result['observed']:.1%}   "
              f"publication-matched null {result['null_mean']:.1%} "
              f"+- {result['null_sd']:.1%}   z {result['z']:+.2f}  "
              f"p {result['p_empirical']:.4f}")
        nulls.append({"confounder": "research intensity",
                      "matched_on": "publication decile", **result})

    # ---------------------------------------------------------------- 3
    print("\n=== 3. Constraint (LOEUF) ===")
    loeuf = load_loeuf(d.index)
    if loeuf is None:
        print("  gnomAD constraint file not present - skipped.")
        print("  Run 00_setup/02_download_references.sh to include it.")
    else:
        d["loeuf"] = loeuf
        s = d.dropna(subset=["loeuf"]).copy()
        print(f"  genes with LOEUF: {len(s):,}")
        print(f"\n{'layer':16s}{'n':>7s}{'median LOEUF':>15s}{'LOEUF < 0.35':>15s}")
        for layer in ("core", "specific", "rest"):
            g = s[s.layer == layer]
            if not len(g):
                continue
            print(f"{layer:16s}{len(g):7,d}{g.loeuf.median():15.3f}"
                  f"{(g.loeuf < 0.35).mean():15.1%}")
        rest = s[s.layer == "rest"]
        for layer in ("core", "specific"):
            g = s[s.layer == layer]
            if len(g) < 5:
                continue
            _, p = stats.mannwhitneyu(g.loeuf, rest.loeuf)
            print(f"  {layer} against rest: Mann-Whitney p {p:.4f}")
            summary.append({"confounder": "constraint",
                            "test": f"LOEUF, {layer} vs rest",
                            "estimate": float(g.loeuf.median() - rest.loeuf.median()),
                            "p": p})
        print("  The core is not, in itself, a constrained gene set. That is the")
        print("  point: if it were, the two-layer result would be a restatement")
        print("  of constraint.")

        s["loeuf_decile"] = pd.qcut(s.loeuf, 10, labels=False, duplicates="drop")
        result = matched_null(s, s.index.intersection(core), "on_panel",
                              "loeuf_decile", args.draws, cfg.SEED + 2)
        print(f"\n  core panel share {result['observed']:.1%}   "
              f"LOEUF-matched null {result['null_mean']:.1%} +- {result['null_sd']:.1%}")
        print(f"  z {result['z']:+.2f}   empirical p {result['p_empirical']:.4f}")
        nulls.append({"confounder": "constraint", "matched_on": "LOEUF decile",
                      **result})
        summary.append({"confounder": "constraint", "test": "LOEUF-matched null",
                        "estimate": result["z"], "p": result["p_empirical"]})

        sub = s.dropna(subset=["regions_skeletal", "n_hpo_terms"]).copy()
        X = pd.DataFrame({
            "breadth": sub.breadth,
            "neg_loeuf": -sub.loeuf,
            "log_hpo": np.log1p(sub.n_hpo_terms),
        }, index=sub.index)
        X = (X - X.mean()) / X.std()
        model = sm.OLS(sub.regions_skeletal, sm.add_constant(X)).fit()
        table = model.summary2().tables[1]
        print(f"\n  The symmetry law with constraint in the model   n {len(sub):,}")
        for term in ("breadth", "neg_loeuf", "log_hpo"):
            print(f"    {term:12s} coef {table.loc[term, 'Coef.']:+.4f}   "
                  f"t {table.loc[term, 't']:+6.2f}   p {table.loc[term, 'P>|t|']:.4g}")
            models.append({"model": "regions ~ breadth + (-LOEUF) + log(HPO)",
                           "term": term, "coef": float(table.loc[term, "Coef."]),
                           "z": float(table.loc[term, "t"]),
                           "p": float(table.loc[term, "P>|t|"]), "n": len(sub)})

    # ---------------------------------------------------------------- 4
    print("\n=== 4. The remaining three, handled elsewhere ===")
    print("  general pleiotropy     -> 04_analysis/03_cross_organ_control.py")
    print("  linkage disequilibrium -> 02_layers/02_locus_clumping.py")
    print("  curation artefact      -> ClinVar benign control in 04_analysis/01")

    common.write_result(pd.DataFrame(summary), "confounder_summary.tsv", index=False)
    common.write_result(pd.DataFrame(nulls), "confounder_matched_nulls.tsv",
                        index=False)
    if models:
        common.write_result(pd.DataFrame(models), "confounder_models.tsv",
                            index=False)
    print("\nNow run:  python 04_analysis/05_replication_gefos.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
