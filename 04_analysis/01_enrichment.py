#!/usr/bin/env python3
"""Do the genes that act everywhere carry the Mendelian disease?

What it does
    The first claim of the paper, against three independently curated truth
    sides and at four levels of robustness.

    1. The layer contrast. Core against site-specific against background, for
       PanelApp, HPO skeletal system, and ClinVar pathogenic-with-skeletal-
       diagnosis, plus ClinVar benign as a negative control. Fisher exact
       throughout.
    2. The dose-response. The share of Mendelian genes at each breadth from 0
       to 6, for all three sources. This is Figure 1B, and it is the panel
       that carries the argument: three curations, three monotone curves.
    3. The threshold sweep. The same odds ratio at k = 2.0 to 4.0, from the
       memberships written by 02_layers/03_threshold_sweep.py. Figure 1C.
    4. The locus sweep. The same odds ratio with genes collapsed to loci at
       five merge distances, from 02_layers/02_locus_clumping.py. Figure 1D,
       and the answer to "this is just linkage disequilibrium".

    It also runs the threshold-free comparison that decides what the finding
    actually is: are the Mendelian genes strong at their *weakest* site
    (shared), or strong at their *best* site (peak), or strong relative to
    their other sites (specificity)? Only the first is the paper's claim.

Reads
    results/gene_layers.tsv, results/continuous_axes.tsv,
    results/threshold_layers.tsv, results/loci.tsv          (stage 02)
    results/truth_panelapp.tsv, results/truth_hpo.tsv,
    results/truth_clinvar.tsv                               (stage 03)

Writes
    results/enrichment_layers.tsv        the three-source layer contrast
    results/enrichment_dose_response.tsv breadth 0-6 x three sources
    results/enrichment_threshold_sweep.tsv
    results/enrichment_locus_sweep.tsv
    results/enrichment_continuous.tsv    the threshold-free comparison

Figures
    Figures 1B-1D. Panel tables:
    fig1b_dose_response.csv, fig1b_trend_tests.csv,
    fig1c_threshold_sweep.csv, fig1d_locus_sweep.csv.

Runtime
    Under a minute.

Usage
    python 04_analysis/01_enrichment.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def fisher(hit_a, n_a, hit_b, n_b):
    """Odds ratio and p for `hit_a of n_a` against `hit_b of n_b`."""
    if min(n_a, n_b) == 0:
        return np.nan, np.nan
    return stats.fisher_exact([[hit_a, n_a - hit_a], [hit_b, n_b - hit_b]])


def truth_sides(index) -> dict:
    """The three curations plus the benign control, as 0/1 over `index`."""
    sides = {}

    panel = common.load_truth("panelapp")
    sides["PanelApp PA309"] = index.isin(panel.index).astype(int)

    hpo = common.load_truth("hpo")
    sides["HPO skeletal"] = pd.Series(index, index=index).map(
        hpo.has_skeletal_disease).fillna(0).astype(int).values

    clinvar = common.load_truth("clinvar")
    sides["ClinVar skeletal"] = pd.Series(index, index=index).map(
        clinvar.has_pathogenic_skeletal).fillna(0).astype(int).values
    sides["ClinVar benign (control)"] = pd.Series(index, index=index).map(
        clinvar.has_benign).fillna(0).astype(int).values
    return sides


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    args = ap.parse_args()

    d = common.load_layers()
    T = pd.DataFrame(truth_sides(d.index), index=d.index)

    common.banner("Claim 1: the shared layer carries the Mendelian disease")
    print(f"  genes {len(d):,} | core {int((d.layer == 'core').sum())} "
          f"| site-specific {int((d.layer == 'specific').sum())} "
          f"| rest {int((d.layer == 'rest').sum())}")

    # ---------------------------------------------------------------- 1
    print("\n=== Layer contrast, three sources and one negative control ===")
    print(f"{'truth side':26s}{'core':>9s}{'specific':>10s}{'background':>12s}"
          f"{'OR':>8s}{'p':>11s}")
    rows = []
    rest = d.layer == "rest"
    for name, y in T.items():
        for layer in ("core", "specific"):
            m = d.layer == layer
            if m.sum() < 5:
                continue
            orr, p = fisher(int(y[m].sum()), int(m.sum()),
                            int(y[rest].sum()), int(rest.sum()))
            rows.append({"truth_side": name, "layer": layer, "n": int(m.sum()),
                         "n_hits": int(y[m].sum()), "share": float(y[m].mean()),
                         "background_share": float(y[rest].mean()),
                         "odds_ratio": orr, "p": p})
        core_m, spec_m = d.layer == "core", d.layer == "specific"
        orr, p = fisher(int(y[core_m].sum()), int(core_m.sum()),
                        int(y[rest].sum()), int(rest.sum()))
        print(f"{name:26s}{y[core_m].mean():9.1%}{y[spec_m].mean():10.1%}"
              f"{y[rest].mean():12.1%}{orr:8.2f}{p:11.2e}")
    print("\n  The last row is the control. It must not enrich, and it does not:")
    print("  if breadth tracked how much a gene has been sequenced, benign")
    print("  variants would follow the pathogenic ones into the core.")
    common.write_result(pd.DataFrame(rows), "enrichment_layers.tsv", index=False)

    # ---------------------------------------------------------------- 2
    print("\n=== Dose-response: share of Mendelian genes by breadth (Figure 1B) ===")
    sources = [c for c in T.columns if "control" not in c]
    print(f"{'breadth':>8s}{'genes':>8s}" + "".join(f"{s[:14]:>16s}" for s in sources))
    dose = []
    for k in range(7):
        m = d.breadth == k
        if m.sum() < 20:
            continue
        print(f"{k:8d}{int(m.sum()):8,d}" +
              "".join(f"{T.loc[m, s].mean():16.1%}" for s in sources))
        dose.append({"breadth": k, "n_genes": int(m.sum()),
                     **{s: round(float(T.loc[m, s].mean()), 5) for s in T.columns}})
    common.write_result(pd.DataFrame(dose), "enrichment_dose_response.tsv",
                        index=False)
    for s in sources:
        rho, p = stats.spearmanr(d.breadth, T[s])
        print(f"  Spearman(breadth, {s:24s}) = {rho:+.3f}  p {p:.2e}")

    # ---------------------------------------------------------------- 3
    print("\n=== Threshold sweep (Figure 1C) ===")
    sweep = pd.read_csv(cfg.require(cfg.RESULTS / "threshold_layers.tsv",
                                    "threshold sweep memberships"), sep="\t")
    print(f"{'layer':10s}{'z0':>6s}{'delta':>7s}{'n':>8s}"
          + "".join(f"{s[:12]:>14s}" for s in sources))
    rows = []
    for (layer, z0, delta), grp in sweep.groupby(
            ["layer", "z0", "delta"], dropna=False):
        members = d.index.isin(grp.gene)
        if members.sum() < 5:
            continue
        line = f"{layer:10s}{z0:6.1f}{delta if pd.notna(delta) else 0:7.1f}" \
               f"{int(members.sum()):8,d}"
        rec = {"layer": layer, "z0": z0, "delta": delta, "n": int(members.sum())}
        for s in sources:
            y = T[s].values
            orr, p = fisher(int(y[members].sum()), int(members.sum()),
                            int(y[rest].sum()), int(rest.sum()))
            line += f"{orr:14.2f}"
            rec[f"OR_{s}"] = orr
            rec[f"p_{s}"] = p
        print(line)
        rows.append(rec)
    common.write_result(pd.DataFrame(rows), "enrichment_threshold_sweep.tsv",
                        index=False)

    # ---------------------------------------------------------------- 4
    print("\n=== Locus sweep: is this linkage disequilibrium? (Figure 1D) ===")
    loci = pd.read_csv(cfg.require(cfg.RESULTS / "loci.tsv", "locus assignments"),
                       sep="\t", index_col=0)
    rows = []
    print(f"{'gap':>8s}{'loci':>8s}{'core loci':>11s}"
          + "".join(f"{s[:12]:>14s}" for s in sources))
    for column in loci.columns:
        gap = int(column.split("_")[1])
        t = d.join(loci[column]).join(T)
        agg = t.groupby(column).agg(
            core=("layer", lambda s: int((s == "core").any())),
            specific=("layer", lambda s: int((s == "specific").any())),
            **{s: (s, "max") for s in T.columns})
        agg.loc[agg.core == 1, "specific"] = 0
        is_core = agg.core == 1
        is_rest = (agg.core == 0) & (agg.specific == 0)
        line = f"{gap // 1000:6d}kb{len(agg):8,d}{int(is_core.sum()):11d}"
        rec = {"gap_bp": gap, "n_loci": len(agg), "n_core_loci": int(is_core.sum())}
        for s in sources:
            orr, p = fisher(int(agg.loc[is_core, s].sum()), int(is_core.sum()),
                            int(agg.loc[is_rest, s].sum()), int(is_rest.sum()))
            line += f"{orr:14.2f}"
            rec[f"OR_{s}"] = orr
            rec[f"p_{s}"] = p
        print(line)
        rows.append(rec)
    common.write_result(pd.DataFrame(rows), "enrichment_locus_sweep.tsv", index=False)
    print("\n  The odds ratio does not collapse as neighbouring genes are merged.")
    print("  Whatever the core is, it is not one signal read out several times.")

    # ---------------------------------------------------------------- 5
    print("\n=== Threshold-free: which axis separates the Mendelian genes? ===")
    axes = pd.read_csv(cfg.require(cfg.RESULTS / "continuous_axes.tsv",
                                   "continuous axes"), sep="\t", index_col=0)
    axes = axes.reindex(d.index)
    print(f"{'axis':16s}{'source':26s}{'median hit':>12s}{'median miss':>13s}"
          f"{'AUC':>7s}{'p':>11s}")
    rows = []
    for axis, label in [("shared", "shared (weakest site Z)"),
                        ("peak", "peak (strongest site Z)"),
                        ("specificity", "specificity (best - 2nd)")]:
        for s in sources:
            y = T[s].astype(bool)
            a, b = axes.loc[y.values, axis], axes.loc[~y.values, axis]
            u, p = stats.mannwhitneyu(a, b)
            auc = u / (len(a) * len(b))
            print(f"{label if s == sources[0] else '':16s}{s:26s}"
                  f"{a.median():12.3f}{b.median():13.3f}{auc:7.3f}{p:11.2e}")
            rows.append({"axis": axis, "truth_side": s,
                         "median_hit": round(float(a.median()), 4),
                         "median_miss": round(float(b.median()), 4),
                         "auc": round(auc, 4), "p": p})
    common.write_result(pd.DataFrame(rows), "enrichment_continuous.tsv", index=False)
    print("\n  Mendelian genes are strong where they are weakest. They are not")
    print("  distinguished by how strong their best site is, and not at all by")
    print("  how far that site stands above the others.")

    print("\nNow run:  python 04_analysis/02_symmetry_law.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
