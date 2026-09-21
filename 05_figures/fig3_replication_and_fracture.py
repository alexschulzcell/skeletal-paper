"""Figure 3 - The same layer replicates in an independent cohort and carries fracture risk.

Reads   results/fig3a_gefos_replication.csv, fig3a_gefos_conditional.csv,
        fig3b_fracture.csv, fig3c_conditional.csv, fig3c_layer_contrast.csv
Writes  figures/figure3.pdf, figures/figure3.tif

Panels
  A  the dose-response curve repeated in GEFOS, a cohort with no UK Biobank samples
  B  breadth of action against fracture risk signal
  C  the two layers side by side, in replication and in the clinic

The P value quoted in panel A is the shared layer against all other genes, from
fig3c_layer_contrast.csv. That is the contrast the main text reports, and both
now read from the same table: fig3a_gefos_stats.csv holds the companion
contrast against every non-shared-layer gene, which is reported in the
supplement instead.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style as S

S.use_paper_style()

rep = pd.read_csv(S.RESULTS / "fig3a_gefos_replication.csv")
gcond = pd.read_csv(S.RESULTS / "fig3a_gefos_conditional.csv").set_index("endpoint")
frac = pd.read_csv(S.RESULTS / "fig3b_fracture.csv")
cond = pd.read_csv(S.RESULTS / "fig3c_conditional.csv").set_index("endpoint")
lay = pd.read_csv(S.RESULTS / "fig3c_layer_contrast.csv")

fig = plt.figure(figsize=(S.W2, 2.62))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.05], wspace=0.50,
                      left=0.072, right=0.988, top=0.800, bottom=0.300)

# ----------------------------------------------- A: independent replication
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "A", dx=-0.215, dy=1.23)
S.headline(ax, "The same curve in an independent cohort", y=1.050, x=-0.150)
for site, col in S.GEFOS.items():
    q = rep[rep.site == site]
    ax.plot(q.breadth, q.median_z, "o-", color=col, mfc="white", mew=1.2,
            ms=3.4, label=site)
ax.axhline(0, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.set_xticks(range(7))
ax.set_xlim(-0.35, 6.35)
ax.set_xlabel("Breadth of action, defined in UK Biobank")
ax.set_ylabel("Median GEFOS gene-level Z")
ax.set_ylim(-0.25, 3.05)
S.ygrid(ax)
ax.legend(loc="upper left", handlelength=1.2, handletextpad=0.4,
          labelspacing=0.22, borderpad=0.0,
          bbox_to_anchor=(-0.012, 1.015))

fn_p = float(lay[(lay.endpoint == "GEFOS femoral neck")
                 & (lay.layer == "Shared layer")].p_vs_rest.iloc[0])
ax.text(0.0, -0.300,
        f"femoral neck, shared layer vs all other genes  {S.pstar(fn_p)}\n"
        f"breadth holds beside all six UK Biobank site scores: "
        f"t = {gcond.loc['GEFOS Femoral neck', 't_breadth']:.1f}",
        transform=ax.transAxes, fontsize=S.SMALL, va="top", ha="left",
        color=S.MUTE, linespacing=1.55)

# -------------------------------------------------------- B: fracture risk
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "B", dx=-0.230, dy=1.23)
S.headline(ax, "And it carries fracture risk", y=1.050, x=-0.165)
cols = [S.SHARED_L] * 6 + [S.SHARED]
ax.bar(frac.breadth, frac.frac_z_gt2 * 100, width=0.68, color=cols, zorder=2)
for r in frac.itertuples():
    ax.text(r.breadth, r.frac_z_gt2 * 100 + 1.4, f"{r.frac_z_gt2 * 100:.0f}",
            ha="center", fontsize=S.SMALL,
            color=S.SHARED if r.breadth == 6 else S.MUTE,
            fontweight="bold" if r.breadth == 6 else "normal")
ax.set_xticks(range(7))
ax.set_xlim(-0.62, 6.62)
ax.set_xlabel("Breadth of action (number of skeletal sites)")
ax.set_ylabel("Genes with fracture-risk\nsignal, Z > 2 (%)")
ax.set_ylim(0, 70)
S.ygrid(ax)
ax.text(0.02, 0.975, "5% → 53%", transform=ax.transAxes,
        fontsize=S.BASE - 0.5, va="top", ha="left", color=S.SHARED,
        fontweight="bold")
ax.text(0.0, -0.300,
        f"breadth holds beside all six site scores: "
        f"t = {cond.loc['Fracture', 't_breadth']:.1f}, "
        f"{S.pstar(cond.loc['Fracture', 'p_breadth'])}",
        transform=ax.transAxes, fontsize=S.SMALL, va="top", ha="left",
        color=S.MUTE, linespacing=1.55)

# ------------------------------------------------ C: the two layers, side by side
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "C", dx=-0.200, dy=1.23)
S.headline(ax, "Only the shared layer travels", y=1.050, x=-0.135)
endpoints = ["GEFOS femoral neck", "GEFOS lumbar spine", "Fracture"]
layers = [("Shared layer", S.SHARED), ("Site-specific layer", S.SITE),
          ("All other genes", S.REST)]
w = 0.26
x = np.arange(len(endpoints))
for k, (name, col) in enumerate(layers):
    v = [lay[(lay.endpoint == e) & (lay.layer == name)].median_z.iloc[0]
         for e in endpoints]
    ax.bar(x + (k - 1) * w, v, width=w * 0.92, color=col, label=name, zorder=2)
ax.axhline(0, color=S.MUTE, lw=0.7)
ax.set_xticks(x)
ax.set_xticklabels(["GEFOS\nfemoral\nneck", "GEFOS\nlumbar\nspine", "Fracture\nrisk"])
ax.set_ylabel("Median gene-level Z")
ax.set_ylim(0, 3.25)
ax.set_xlim(-0.55, 2.55)
S.ygrid(ax)
# Headroom is added above the tallest bar so the legend can sit inside the
# panel without the bar reaching into its third entry.
ax.legend(loc="upper right", handlelength=0.9, handletextpad=0.35,
          labelspacing=0.22, borderpad=0.0, fontsize=S.SMALL,
          bbox_to_anchor=(1.01, 1.02))
ax.text(1.0, -0.375,
        "n = 103, 36 and 18,245 genes",
        transform=ax.transAxes, fontsize=S.SMALL, va="top", ha="right",
        color=S.MUTE)

S.save(fig, "figure3")
