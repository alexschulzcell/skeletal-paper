"""Figure 4 - What the shared layer is made of.

Reads   results/fig4a_pathways.csv, fig4b_regulatory_density.csv,
        fig4c_shared_layer_genes.csv
Writes  figures/figure4.pdf, figures/figure4.tif

Panels
  A  developmental signaling pathways in the shared layer
  B  regulatory density around each layer, against a length-matched expectation
  C  the shared layer gene by gene, across the two clinical ends

Panel C labels exactly the genes the main text names, and no others, so that
the sentence and the panel can be checked against each other.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import style as S

S.use_paper_style()

pw = pd.read_csv(S.RESULTS / "fig4a_pathways.csv")
rd = pd.read_csv(S.RESULTS / "fig4b_regulatory_density.csv").set_index("layer")
gn = pd.read_csv(S.RESULTS / "fig4c_shared_layer_genes.csv")

fig = plt.figure(figsize=(S.W2, 2.85))
gs = fig.add_gridspec(1, 3, width_ratios=[1.10, 0.72, 1.24], wspace=0.50,
                      left=0.168, right=0.988, top=0.845, bottom=0.205)

# ------------------------------------------------------- A: pathway forest
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "A", dx=-0.660, dy=1.155)
S.headline(ax, "A Wnt and ossification module", y=1.040, x=-0.600)
sh = pw[pw.layer == "Shared layer"].sort_values("odds_ratio")
y = np.arange(len(sh))
for yi, r in zip(y, sh.itertuples()):
    sig = r.p < 0.05
    ax.plot([r.ci_low, r.ci_high], [yi, yi], color=S.SHARED if sig else S.REST,
            lw=1.2, solid_capstyle="round", zorder=2)
    ax.plot(r.odds_ratio, yi, "o", color=S.SHARED if sig else S.REST, ms=4.2,
            mfc=S.SHARED if sig else "white", mew=1.2, zorder=3)
    ax.text(r.ci_high * 1.18, yi, f"{r.n_hits}/103", va="center",
            fontsize=S.SMALL, color=S.MUTE if sig else S.REST)
ax.axvline(1, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.set_yticks(y)
ax.set_yticklabels([p.replace("TGF-beta", "TGF-β") for p in sh.pathway],
                   fontsize=S.SMALL)
ax.set_xscale("log")
ax.set_xticks([0.3, 1, 3, 10, 30])
ax.set_xticklabels(["0.3", "1", "3", "10", "30"])
ax.minorticks_off()
ax.set_xlim(0.18, 145)
ax.set_ylim(-0.75, len(y) - 0.25)
ax.set_xlabel("Odds ratio, against all other genes")
ax.xaxis.grid(True, color=S.GRID, lw=0.5)
ax.set_axisbelow(True)

# ------------------------------------------------- B: regulatory density
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "B", dx=-0.600, dy=1.155)
S.headline(ax, "Denser regulation", y=1.040, x=-0.530)
# An observation against a null interval, drawn as an observation against a
# null interval. Bars were wrong here twice over: they implied a count from
# zero, and the axis does not start at zero, so their lengths meant nothing.
for i, (layer, col) in enumerate([("Shared layer", S.SHARED),
                                  ("Site-specific layer", S.SITE)]):
    r = rd.loc[layer]
    lo, hi = r.expected - 1.96 * r.sd, r.expected + 1.96 * r.sd
    ax.vlines(i, lo, hi, color=S.REST, lw=1.1, zorder=2)
    for y in (lo, hi):
        ax.hlines(y, i - 0.10, i + 0.10, color=S.REST, lw=1.1, zorder=2)
    ax.hlines(r.expected, i - 0.17, i + 0.17, color=S.MUTE, lw=1.4, zorder=3)
    ax.plot(i, r.observed_peaks_per_element, "o", color=col, ms=5.5,
            mec="white", mew=0.8, zorder=4)
    ax.annotate(f"Z = {r.z:+.1f}" if r.p < 0.05 else "n.s.",
                (i, r.observed_peaks_per_element), xytext=(0, 7),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=S.BASE - 1.0, fontweight="bold",
                color=col if r.p < 0.05 else S.MUTE)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Shared\nlayer", "Site-\nspecific"])
ax.set_xlim(-0.65, 1.65)
ax.set_ylim(9.0, 14.6)
ax.set_ylabel("Open-chromatin peaks per\nregulatory element")
S.ygrid(ax)
ax.text(0.5, -0.175, "grey: length-matched null,\nmean and 95% interval",
        transform=ax.transAxes, ha="center", va="top", fontsize=S.SMALL,
        color=S.MUTE, linespacing=1.5)

# ---------------------------------------- C: the shared layer, gene by gene
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "C", dx=-0.150, dy=1.155)
S.headline(ax, "One layer, two clinical ends", y=1.040, x=-0.095)
men = gn.mendelian == 1
# The background is drawn lighter than the named genes, so the labelled
# anchors stand out of a hundred points instead of competing with them.
ax.scatter(gn.loc[~men, "FRAK"], gn.loc[~men, "gefos_fn_z"], s=10,
           facecolor="white", edgecolor=S.SHARED, linewidth=0.6, alpha=0.75,
           zorder=3, label="no Mendelian disease")
ax.scatter(gn.loc[men, "FRAK"], gn.loc[men, "gefos_fn_z"], s=13,
           color=S.HIGHLIGHT, edgecolor="white", linewidth=0.35, alpha=0.85,
           zorder=4, label="Mendelian skeletal disease")
ax.axvline(2, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.axhline(2, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.set_xlabel("Fracture-risk signal, gene-level Z")
ax.set_ylabel("Independent replication,\nGEFOS femoral neck Z")
ax.set_xlim(-2.6, 9.2)
ax.set_ylim(-2.2, 10.4)
ax.legend(loc="upper left", handletextpad=0.2, labelspacing=0.22,
          borderpad=0.0, bbox_to_anchor=(-0.012, 1.015), fontsize=S.SMALL)
ax.text(1.0, -0.195,
        f"{int(men.sum())} of the {len(gn)} shared-layer genes\n"
        f"carry a Mendelian skeletal disease",
        transform=ax.transAxes, ha="right", va="top", fontsize=S.SMALL,
        color=S.MUTE, linespacing=1.5)

# The genes the text names: the osteoporosis-genetics textbook, including both
# halves of the OPG-RANK pair.
NAMED = ["ESR1", "LRP5", "WNT16", "WNT4", "WLS", "TNFRSF11B", "TNFRSF11A",
         "SOX6", "CPED1", "FGFRL1", "IDUA", "MDK", "NOTUM"]
lab = gn[gn.gene.isin(NAMED)]
missing = sorted(set(NAMED) - set(lab.gene))
if missing:
    raise SystemExit(f"figure4: genes named in the text are not in the shared "
                     f"layer table: {missing}")
# Ring the named genes so that the end of every leader is unmistakable, then
# set all the names in one ink colour: the marker already carries the Mendelian
# status, and two label colours on top of it only added noise.
ax.scatter(lab.FRAK, lab.gefos_fn_z, s=34, facecolor="none",
           edgecolor=S.INK, linewidth=0.7, zorder=5)
S.place_labels(ax, lab.FRAK.to_numpy(), lab.gefos_fn_z.to_numpy(),
               lab.gene.tolist(), [S.INK] * len(lab),
               fontsize=S.SMALL, marker_pt=3.4,
               obstacles=gn[["FRAK", "gefos_fn_z"]].to_numpy(),
               leader_color=S.MUTE)

S.save(fig, "figure4")
