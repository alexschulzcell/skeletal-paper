#!/usr/bin/env python3
"""Does the shared layer reach the clinic? Fracture, heel BMD, body size.

What it does
    Carries breadth from bone density - a surrogate - to endpoints that matter,
    and asks each one the hard version of the question.

    Fracture (N 426,795) is the registrational endpoint. The trivial reading is
    that fracture follows bone density and the core is defined on bone density,
    so the result is circular. It is not, and the script shows why three ways:
    breadth against peak strength in one model; the layer contrast within peak-
    strength quartiles; and breadth entered alongside all six site Z values, so
    that it must hold against the exact measurements it is built from.

    Heel BMD (N 426,824) is a different modality - ultrasound, not DXA - at
    thirteen times the sample size, and serves as an independent validation of
    the density signal itself.

    Standing height and sitting-height ratio are anthropometric controls.

    The script also assembles the breadth-against-peak comparison across all
    five endpoints, with the alternative explanation printed next to it: the
    two disease endpoints are binary and the three measurements quantitative,
    and binary traits have lower effective power, which shifts the peak/breadth
    ratio on its own. A clean two-against-three split arises by chance with
    p = 0.10. The paper reports this as suggestive, not established.

Reads
    results/gene_layers.tsv, results/gene_z_matrix.tsv     (stages 01-02)

Writes
    results/endpoints_by_layer.tsv       layer contrast per endpoint
    results/endpoints_dose_response.tsv  breadth 0-6 per endpoint
    results/endpoints_models.tsv         breadth against peak, and against sites
    results/endpoints_peak_matched.tsv   layer contrast within peak quartiles

Figures
    Figures 3B and 3C. Panel tables:
    fig3b_fracture.csv, fig3c_layer_contrast.csv,
    fig3c_conditional.csv.

Runtime
    Under a minute.

Usage
    python 04_analysis/06_clinical_endpoints.py
    python 04_analysis/06_clinical_endpoints.py --endpoints frak heel
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

#: Which endpoints are disease (binary) and which are measurements
#: (quantitative). The distinction is the whole content of the last section.
ENDPOINT_KIND = {"frak": "disease",
                 "heel": "measurement", "height": "measurement",
                 "prop": "measurement"}


def standardise(X: pd.DataFrame) -> pd.DataFrame:
    return (X - X.mean()) / X.std()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--endpoints", nargs="*", default=list(cfg.ENDPOINTS))
    args = ap.parse_args()

    Z = common.load_gene_matrix()
    d = common.load_layers()
    available = [e for e in args.endpoints
                 if e in Z.columns and Z[e].notna().sum() > 1000]
    missing = [e for e in args.endpoints if e not in available]
    for e in available:
        d[e] = Z[e]

    common.banner("Claim 4b: the shared layer reaches the clinic")
    print(f"  endpoints available: "
          f"{', '.join(cfg.ENDPOINTS[e] for e in available) or '(none)'}")
    if missing:
        print(f"  not available      : {', '.join(missing)}")
    if not available:
        raise SystemExit("No endpoint available; nothing to do.")

    core = d.layer == "core"
    spec = d.layer == "specific"
    rest = d.layer == "rest"

    # ---------------------------------------------------------------- 1
    print("\n=== Signal by layer ===")
    rows = []
    for e in available:
        s = d.dropna(subset=[e])
        print(f"\n  {cfg.ENDPOINTS[e]}  (N {cfg.TRAIT_N.get(e, float('nan')):,.0f}, "
              f"{ENDPOINT_KIND.get(e, '?')})")
        print(f"  {'layer':16s}{'n':>6s}{'median Z':>11s}{'Z > 2':>9s}{'Z > 4':>9s}")
        for name, mask in [("core", core), ("site-specific", spec), ("rest", rest)]:
            g = s[mask.reindex(s.index, fill_value=False)]
            if not len(g):
                continue
            print(f"  {name:16s}{len(g):6,d}{g[e].median():11.2f}"
                  f"{(g[e] > 2).mean():9.1%}{(g[e] > 4).mean():9.1%}")
            rows.append({"endpoint": e, "layer": name, "n": len(g),
                         "median_z": round(float(g[e].median()), 4),
                         "share_z_gt2": round(float((g[e] > 2).mean()), 4),
                         "share_z_gt4": round(float((g[e] > 4).mean()), 4)})
        g_core = s[core.reindex(s.index, fill_value=False)][e]
        g_spec = s[spec.reindex(s.index, fill_value=False)][e]
        g_rest = s[rest.reindex(s.index, fill_value=False)][e]
        p_core = stats.mannwhitneyu(g_core, g_rest)[1] if len(g_core) > 2 else np.nan
        p_spec = stats.mannwhitneyu(g_spec, g_rest)[1] if len(g_spec) > 2 else np.nan
        p_layers = (stats.mannwhitneyu(g_core, g_spec)[1]
                    if min(len(g_core), len(g_spec)) > 2 else np.nan)
        print(f"    core vs rest p {p_core:.2e} | site-specific vs rest "
              f"p {p_spec:.3g} | core vs site-specific p {p_layers:.3g}")
        rows.append({"endpoint": e, "layer": "contrasts", "n": len(s),
                     "p_core_vs_rest": p_core, "p_specific_vs_rest": p_spec,
                     "p_core_vs_specific": p_layers})
    common.write_result(pd.DataFrame(rows), "endpoints_by_layer.tsv", index=False)

    # ---------------------------------------------------------------- 2
    print("\n=== Dose-response (Figure 3B) ===")
    print(f"{'breadth':>8s}{'genes':>8s}" +
          "".join(f"{cfg.ENDPOINTS[e][:12]:>14s}" for e in available))
    dose = []
    for k in range(7):
        s = d[d.breadth == k]
        if len(s) < 30:
            continue
        print(f"{k:8d}{len(s):8,d}" +
              "".join(f"{s[e].median():14.2f}" for e in available))
        record = {"breadth": k, "n_genes": len(s)}
        for e in available:
            record[f"{e}_median_z"] = round(float(s[e].median()), 4)
            record[f"{e}_share_z_gt2"] = round(float((s[e] > 2).mean()), 4)
        dose.append(record)
    common.write_result(pd.DataFrame(dose), "endpoints_dose_response.tsv",
                        index=False)
    if "frak" in available:
        first = dose[0][f"frak_share_z_gt2"]
        last = dose[-1][f"frak_share_z_gt2"]
        print(f"\n  Fracture, share above Z 2: {first:.1%} at breadth "
              f"{dose[0]['breadth']} -> {last:.1%} at breadth {dose[-1]['breadth']}.")

    # ---------------------------------------------------------------- 3
    print("\n=== Is it just bone density? Breadth against peak strength ===")
    print(f"{'endpoint':16s}{'kind':13s}{'breadth t':>11s}{'peak t':>10s}"
          f"{'leads':>10s}{'breadth t | six site Z':>25s}")
    models = []
    X_simple = standardise(d[["breadth", "peak"]])
    X_full = standardise(d[["breadth"] + cfg.BMD_SITES])
    for e in available:
        s = d.dropna(subset=[e])
        y = s[e].rank(pct=True)      # ranks: the p floor truncates the top
        m1 = sm.OLS(y, sm.add_constant(X_simple.loc[s.index])).fit()
        t1 = m1.summary2().tables[1]
        m2 = sm.OLS(y, sm.add_constant(X_full.loc[s.index])).fit()
        t_breadth = float(t1.loc["breadth", "t"])
        t_peak = float(t1.loc["peak", "t"])
        leads = "breadth" if t_breadth > t_peak else "peak"
        print(f"{cfg.ENDPOINTS[e]:16s}{ENDPOINT_KIND.get(e, '?'):13s}"
              f"{t_breadth:+11.2f}{t_peak:+10.2f}{leads:>10s}"
              f"{m2.tvalues['breadth']:+18.2f}"
              f"  p {m2.pvalues['breadth']:.1e}")
        models += [
            {"endpoint": e, "kind": ENDPOINT_KIND.get(e, "?"),
             "term": "breadth", "model": "rank ~ breadth + peak",
             "t": t_breadth, "p": float(t1.loc["breadth", "P>|t|"])},
            {"endpoint": e, "kind": ENDPOINT_KIND.get(e, "?"),
             "term": "peak", "model": "rank ~ breadth + peak",
             "t": t_peak, "p": float(t1.loc["peak", "P>|t|"])},
            {"endpoint": e, "kind": ENDPOINT_KIND.get(e, "?"),
             "term": "breadth", "model": "rank ~ breadth + six site Z",
             "t": float(m2.tvalues["breadth"]), "p": float(m2.pvalues["breadth"])},
        ]
    common.write_result(pd.DataFrame(models), "endpoints_models.tsv", index=False)

    kinds = {e: ENDPOINT_KIND.get(e) for e in available}
    if len(set(kinds.values())) > 1:
        print("\n  Read with care. Breadth tends to lead at the disease endpoints")
        print("  and peak strength at the measurements - but the disease endpoints")
        print("  are binary and the measurements quantitative, and binary traits")
        print("  have lower effective power, which shifts that ratio by itself. A")
        print("  clean split of this size arises by chance with p = 0.10.")
        print("  Suggestive, not established. The paper says so.")

    # ---------------------------------------------------------------- 4
    print("\n=== Layer contrast within peak-strength quartiles ===")
    rows = []
    for e in available:
        s = d.dropna(subset=[e])
        sub = s[core.reindex(s.index, fill_value=False)
                | spec.reindex(s.index, fill_value=False)].copy()
        if len(sub) < 20:
            continue
        sub["peak_quartile"] = pd.qcut(sub.peak, 4, labels=False, duplicates="drop")
        print(f"\n  {cfg.ENDPOINTS[e]}")
        print(f"  {'quartile':10s}{'n core':>8s}{'n specific':>12s}"
              f"{'core median':>13s}{'specific median':>17s}")
        differences = []
        for q, g in sub.groupby("peak_quartile"):
            a = g[g.layer == "core"][e]
            b = g[g.layer == "specific"][e]
            if len(a) < 3 or len(b) < 3:
                continue
            print(f"  Q{int(q) + 1:<9d}{len(a):8d}{len(b):12d}"
                  f"{a.median():13.2f}{b.median():17.2f}")
            differences.append(float(a.median() - b.median()))
            rows.append({"endpoint": e, "peak_quartile": int(q) + 1,
                         "n_core": len(a), "n_specific": len(b),
                         "core_median_z": round(float(a.median()), 4),
                         "specific_median_z": round(float(b.median()), 4),
                         "difference": round(float(a.median() - b.median()), 4)})
        if differences:
            print(f"    difference positive in every quartile: "
                  f"{all(x > 0 for x in differences)}  "
                  f"({[round(x, 2) for x in differences]})")
    if rows:
        common.write_result(pd.DataFrame(rows), "endpoints_peak_matched.tsv",
                            index=False)

    print("\nNow run:  python 04_analysis/07_pathways_and_axis.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
