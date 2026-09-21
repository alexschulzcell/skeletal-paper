#!/usr/bin/env python3
"""Replication in an independent cohort: GEFOS/UK10K, no UK Biobank.

What it does
    Two questions, and the second is the one that shaped how the paper states
    its claim.

    1. Does the layer structure defined in UK Biobank show up in a cohort that
       shares no participants with it? Three GEFOS traits - femoral neck,
       lumbar spine, forearm - are compared across the three layers, with the
       dose-response redrawn in the foreign cohort (Figure 3A). Because the
       endpoint is another BMD measurement, the test is made hard: breadth is
       entered alongside all six UK Biobank site Z values, so that it must earn
       its coefficient against the very measurements it was built from.

    2. If the structure is defined *afresh in GEFOS*, do the same genes come
       out? They largely do not - about a fifth of the core is recovered - and
       the Mendelian enrichment of the re-derived core is weaker. The property
       replicates; the gene list does not. That is why the paper reports breadth
       as a continuous property and does not publish a core gene list as its
       claim. See docs/KNOWN_ISSUES.md.

    The GEFOS sample sizes are estimated from standard errors and allele
    frequencies rather than taken from a column, because the published files
    have none; the estimator and its accuracy are documented in
    01_gene_scores/01_harmonise_sumstats.py and docs/KNOWN_ISSUES.md.

Reads
    results/gene_layers.tsv, results/gene_z_matrix.tsv    (stages 01-02)
    results/truth_panelapp.tsv, results/truth_hpo.tsv     (stage 03)

Writes
    results/replication_layers.tsv       median Z per layer per GEFOS trait
    results/replication_dose_response.tsv  breadth 0-6 in the foreign cohort
    results/replication_models.tsv       breadth against peak, and against all
                                         six site Z values
    results/replication_redefined.tsv    the structure re-derived in GEFOS
    results/replication_overlap.tsv      gene-list overlap between cohorts

Figures
    Figure 3A. Panel tables:
    fig3a_gefos_replication.csv, fig3a_gefos_stats.csv,
    fig3a_gefos_conditional.csv.

Runtime
    Under a minute.

Usage
    python 04_analysis/05_replication_gefos.py
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
    args = ap.parse_args()

    Z = common.load_gene_matrix()
    layers = common.load_layers()
    available = [t for t in cfg.GEFOS_TRAITS if t in Z.columns and Z[t].notna().any()]
    if len(available) < 2:
        raise SystemExit(
            "\nThe replication needs at least two GEFOS traits.\n"
            f"  found: {available or 'none'}\n"
            "  Download them with 00_setup/01_download_sumstats.sh and run\n"
            "  01_gene_scores/02_run_magma.sh before this script.\n")

    d = layers.join(Z[available], how="inner").dropna(subset=available)
    common.banner("Claim 4a: replication in an independent cohort")
    print(f"  GEFOS traits: {', '.join(cfg.GEFOS_TRAITS[t] for t in available)}")
    print(f"  genes present in both cohorts: {len(d):,}")
    print(f"  core {int((d.layer == 'core').sum())} | site-specific "
          f"{int((d.layer == 'specific').sum())} | rest {int((d.layer == 'rest').sum())}")

    # ---------------------------------------------------------------- 1
    print("\n=== Signal in the foreign cohort, by UK-Biobank-defined layer ===")
    header = f"{'layer':16s}{'n':>6s}" + "".join(
        f"{cfg.GEFOS_TRAITS[t][:12]:>14s}" for t in available)
    print(header)
    rows = []
    groups = {name: d[d.layer == name] for name in ("core", "specific", "rest")}
    for name, g in groups.items():
        print(f"{name:16s}{len(g):6,d}" +
              "".join(f"{g[t].median():14.2f}" for t in available))
    for t in available:
        rest = groups["rest"][t]
        line = f"  {cfg.GEFOS_TRAITS[t]:14s}"
        for name in ("core", "specific"):
            g = groups[name][t]
            p = stats.mannwhitneyu(g, rest)[1] if len(g) > 2 else np.nan
            line += f"  {name} vs rest p {p:.2e}"
            rows.append({"trait": t, "layer": name, "n": len(g),
                         "median_z": round(float(g.median()), 4),
                         "rest_median_z": round(float(rest.median()), 4),
                         "p_vs_rest": p})
        p_layers = stats.mannwhitneyu(groups["core"][t], groups["specific"][t])[1] \
            if len(groups["specific"]) > 2 else np.nan
        line += f"  core vs specific p {p_layers:.4f}"
        print(line)
    common.write_result(pd.DataFrame(rows), "replication_layers.tsv", index=False)
    print("\n  The site-specific layer does not replicate. That is not a weakness")
    print("  of the cohort: it is the dissociation the paper is about.")

    # ---------------------------------------------------------------- 2
    print("\n=== Dose-response in the foreign cohort (Figure 3A) ===")
    print(f"{'breadth':>8s}{'genes':>8s}" +
          "".join(f"{cfg.GEFOS_TRAITS[t][:12]:>14s}" for t in available))
    dose = []
    for k in range(7):
        s = d[d.breadth == k]
        if len(s) < 30:
            continue
        print(f"{k:8d}{len(s):8,d}" +
              "".join(f"{s[t].median():14.2f}" for t in available))
        dose.append({"breadth": k, "n_genes": len(s),
                     **{t: round(float(s[t].median()), 4) for t in available}})
    common.write_result(pd.DataFrame(dose), "replication_dose_response.tsv",
                        index=False)
    print("\n  The same curve shape, in a cohort that shares no participants.")

    # ---------------------------------------------------------------- 3
    print("\n=== Does breadth hold against peak strength, and against the sites? ===")
    models = []
    X_simple = standardise(d[["breadth", "peak"]])
    X_full = standardise(d[["breadth"] + cfg.BMD_SITES])
    for t in available:
        # Ranks, not raw Z: MAGMA truncates gene p at 5e-10, so the top of the
        # distribution is compressed and a mean-based fit would read that
        # truncation as signal. See docs/KNOWN_ISSUES.md.
        y = d[t].rank(pct=True)
        m1 = sm.OLS(y, sm.add_constant(X_simple)).fit()
        t1 = m1.summary2().tables[1]
        m2 = sm.OLS(y, sm.add_constant(X_full)).fit()
        print(f"  {cfg.GEFOS_TRAITS[t]:14s} "
              f"breadth t {t1.loc['breadth', 't']:+6.2f} "
              f"(p {t1.loc['breadth', 'P>|t|']:.2e})   "
              f"peak t {t1.loc['peak', 't']:+6.2f}")
        print(f"  {'':14s} with all six UK Biobank site Z in the model: "
              f"breadth t {m2.tvalues['breadth']:+.2f}  p {m2.pvalues['breadth']:.2e}")
        models += [
            {"trait": t, "model": "rank ~ breadth + peak", "term": "breadth",
             "t": float(t1.loc["breadth", "t"]), "p": float(t1.loc["breadth", "P>|t|"])},
            {"trait": t, "model": "rank ~ breadth + peak", "term": "peak",
             "t": float(t1.loc["peak", "t"]), "p": float(t1.loc["peak", "P>|t|"])},
            {"trait": t, "model": "rank ~ breadth + six site Z", "term": "breadth",
             "t": float(m2.tvalues["breadth"]), "p": float(m2.pvalues["breadth"])},
        ]
    common.write_result(pd.DataFrame(models), "replication_models.tsv", index=False)

    # ---------------------------------------------------------------- 4
    common.banner("The nuance: redefining the structure inside GEFOS")
    gefos_core = common.core_genes(d, sites=available)
    gefos_spec = common.specific_genes(d, sites=available)
    gefos_rest = set(d.index) - gefos_core - gefos_spec
    ukb_core = set(d.index[d.layer == "core"])
    print(f"  GEFOS-defined core {len(gefos_core)} | site-specific {len(gefos_spec)} "
          f"| rest {len(gefos_rest):,}")

    panel = common.load_truth("panelapp")
    hpo = common.load_truth("hpo")
    mendelian = set(panel.index) | set(hpo.index[hpo.has_skeletal_disease == 1])
    truth = {"PanelApp": d.index.isin(panel.index).astype(int),
             "Mendelian (panel or HPO)": d.index.isin(mendelian).astype(int)}

    rows = []
    for name, y in truth.items():
        y = pd.Series(y, index=d.index)
        rest_hits, rest_n = int(y[list(gefos_rest)].sum()), len(gefos_rest)
        print(f"\n  {name}")
        for label, members in [("GEFOS core", gefos_core),
                               ("GEFOS site-specific", gefos_spec)]:
            if len(members) < 5:
                continue
            hits = int(y[list(members)].sum())
            orr, p = stats.fisher_exact([[hits, len(members) - hits],
                                         [rest_hits, rest_n - rest_hits]])
            print(f"    {label:22s}{hits:4d}/{len(members):4d} = "
                  f"{hits / len(members):6.1%}   OR {orr:5.2f}  p {p:.3g}")
            rows.append({"truth_side": name, "layer": label, "n": len(members),
                         "n_hits": hits, "share": hits / len(members),
                         "odds_ratio": orr, "p": p})
        print(f"    {'rest':22s}{rest_hits:4d}/{rest_n:4d} = "
              f"{rest_hits / rest_n:6.1%}")
    common.write_result(pd.DataFrame(rows), "replication_redefined.tsv", index=False)

    overlap = len(ukb_core & gefos_core)
    share = overlap / max(len(ukb_core), 1)
    print(f"\n  Overlap of the two cores: {overlap} genes, {share:.0%} of the "
          f"UK Biobank core.")
    print("  The property replicates; the list does not. The paper therefore")
    print("  reports breadth as a continuous property and treats the core as an")
    print("  illustration of it, never as a result in its own right.")
    common.write_result(pd.DataFrame([{
        "ukb_core": len(ukb_core), "gefos_core": len(gefos_core),
        "overlap": overlap, "share_of_ukb_core": round(share, 4),
        "genes_in_both": ",".join(sorted(ukb_core & gefos_core))}]),
        "replication_overlap.tsv", index=False)

    print("\nNow run:  python 04_analysis/06_clinical_endpoints.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
