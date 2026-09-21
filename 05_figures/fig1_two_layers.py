"""Figure 1 - Breadth of action across the skeleton marks the Mendelian layer.

Reads   results/fig1b_dose_response.csv, fig1b_trend_tests.csv,
        fig1c_threshold_sweep.csv, fig1d_locus_sweep.csv
Writes  figures/figure1.pdf, figures/figure1.tif

Panels
  A  how breadth is defined: six anatomically named DXA sites, one count per gene
  B  breadth against Mendelian skeletal disease, in three independent curations
  C  the same odds ratio across five definitions of the shared layer
  D  the same odds ratio after collapsing genes into independent loci
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch

import anatomy as A
import style as S

S.use_paper_style()

b = pd.read_csv(S.RESULTS / "fig1b_dose_response.csv")
btr = pd.read_csv(S.RESULTS / "fig1b_trend_tests.csv").set_index("source")
c = pd.read_csv(S.RESULTS / "fig1c_threshold_sweep.csv")
d = pd.read_csv(S.RESULTS / "fig1d_locus_sweep.csv")

fig = plt.figure(figsize=(S.W2, 4.25))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.30], height_ratios=[1.0, 0.92],
                      wspace=0.46, hspace=0.62,
                      left=0.088, right=0.986, top=0.912, bottom=0.115)

# ------------------------------------------------------------------ A: design
ax = fig.add_subplot(gs[0, 0])
ax.set_axis_off()
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
S.panel(ax, "A", dx=-0.070, dy=1.15)
S.headline(ax, "One number per gene: at how many sites does it act?",
           y=1.045, x=0.02)

# The body lives in its own axes so that it can hold a true 1:1 aspect ratio
# without forcing the rest of the panel into a square.
bax = ax.inset_axes([-0.04, 0.10, 0.52, 0.92])
bax.set_axis_off()
S.fit_aspect(bax, (-2, 102))

A.draw_body(bax)
A.draw_spine(bax)

# Labels sit in two columns clear of the body; a hairline connects each to its
# marker, so no label has to be read by proximity.
LABELS = {                        # site: (text x, text y, alignment)
    "Head":         (10.5, 95.0, "left"),
    "Lumbar spine": (10.5, 73.0, "left"),
    "Pelvis":       (10.5, 60.0, "left"),
    "Femoral neck": (10.5, 47.0, "left"),
    "Arms":         (-13.5, 66.0, "right"),
    "Legs":         (-13.5, 22.0, "right"),
}
for name in A.SITE_ORDER:
    mx, my = A.site_xy(name)
    tx, ty, ha = LABELS[name]
    step = 2.6 if ha == "left" else -2.6
    bax.plot([mx, mx + step, tx - step * 0.30], [my, ty, ty],
             color=S.REST, lw=0.5, zorder=3, solid_capstyle="round")
    bax.add_patch(Circle((mx, my), 1.7, fc=S.SHARED, ec="white", lw=0.7, zorder=5))
    bax.text(tx, ty, name, fontsize=S.SMALL, color=S.MUTE, ha=ha, va="center")

ax.text(0.22, 0.015, "UK Biobank DXA, n = 31,986 per site",
        transform=ax.transAxes, fontsize=S.SMALL, color=S.MUTE,
        ha="center", va="bottom")

# --- the two example genes --------------------------------------------------
# Markers are sized in points, so the columns stay circular whatever the
# panel's aspect ratio turns out to be.
for j, (label, filled, col, note) in enumerate(
        [("acts at all six", 6, S.SHARED, "breadth 6"),
         ("acts at one", 1, S.SITE, "breadth 1")]):
    x0 = 62.0 + j * 27.0
    ax.text(x0, 96.0, label, fontsize=S.BASE - 0.8, ha="center", va="center",
            color=col, fontweight="bold")
    ys = [84.0 - k * 9.6 for k in range(6)]
    ax.scatter([x0] * filled, ys[:filled], s=30, marker="o",
               facecolor=col, edgecolor="white", linewidth=0.7, zorder=3)
    ax.scatter([x0] * (6 - filled), ys[filled:], s=30, marker="o",
               facecolor="white", edgecolor=S.REST, linewidth=0.8, zorder=3)
    ax.add_patch(FancyBboxPatch((x0 - 11.0, 12.0), 22.0, 9.0,
                                boxstyle="round,pad=0,rounding_size=2.4",
                                fc=col, ec="none", alpha=0.17, zorder=2))
    ax.text(x0, 16.5, note, fontsize=S.BASE - 0.5, ha="center", va="center",
            color=col, fontweight="bold", zorder=3)

# ---------------------------------------------------- B: the dose-response
# The right third of the cell is reserved for the direct labels, so that they
# never overlap the curves and never run off the canvas.
sub = gs[0, 1].subgridspec(1, 2, width_ratios=[1.0, 0.48], wspace=0.0)
ax = fig.add_subplot(sub[0, 0])
S.panel(ax, "B", dx=-0.185, dy=1.15)
S.headline(ax, "Wider action, more Mendelian disease — all three curations agree",
           y=1.045, x=-0.105)

ax.set_yscale("log")
ax.set_ylim(0.82, 7.6)
ax.set_xlim(-0.28, 6.28)

ENDS = {}
for name, col in S.TRUTH.items():
    q = b[b.source == name].copy()
    base = q.fraction.iloc[0]
    ax.fill_between(q.breadth, q.ci_low / base, q.ci_high / base,
                    color=col, alpha=0.13, lw=0)
    ax.plot(q.breadth, q.fraction / base, "o-", color=col, mfc="white", mew=1.2,
            ms=3.4, clip_on=False, zorder=4)
    ENDS[name] = (q.fraction.iloc[-1] / base, q.fraction.iloc[-1] * 100)

# Two of the three curves end close together, so the labels are pushed apart
# in axes space and joined to their curve by a hairline. Nothing is nudged by
# hand: change a number and the labels re-space themselves.
names = list(ENDS)
to_axes = ax.transAxes.inverted().transform
ideal = [to_axes(ax.transData.transform((6, ENDS[n][0])))[1] for n in names]
ax_h_in = ax.get_position().height * fig.get_figheight()
gap = 2 * 1.38 * (S.BASE - 1.2) / 72 / ax_h_in
placed = S.spread_labels(ideal, gap, gap / 2, 1 - gap / 2)
for name, y0, y1 in zip(names, ideal, placed):
    ax.annotate(f"{name}\n{ENDS[name][1]:.0f}% at breadth 6",
                xy=(1.015, y0), xytext=(1.055, y1),
                xycoords="axes fraction", textcoords="axes fraction",
                color=S.TRUTH[name], fontweight="bold", fontsize=S.BASE - 1.2,
                va="center", ha="left", linespacing=1.35, annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color=S.TRUTH[name], lw=0.5,
                                shrinkA=0, shrinkB=1.5))

ax.axhline(1, color=S.REST, lw=0.8, ls=(0, (3, 2)), zorder=1)
ax.set_xticks(range(7))
ax.set_yticks([1, 1.5, 2, 3, 5, 7])
ax.set_yticklabels(["1", "1.5", "2", "3", "5", "7"])
ax.minorticks_off()
ax.set_xlabel("Breadth of action (number of skeletal sites)")
ax.set_ylabel("Mendelian disease, relative to\ngenes that act at no site")
S.ygrid(ax)

tx = "\n".join(
    f"{n}  OR {btr.loc[n, 'or_per_site']:.2f} per site   {S.pstar(btr.loc[n, 'p_trend'])}"
    for n in S.TRUTH)
ax.text(0.02, 0.985, tx, transform=ax.transAxes, fontsize=S.SMALL,
        va="top", ha="left", color=S.MUTE, linespacing=1.55, zorder=6)

# ------------------------------------------------ C: threshold sensitivity
ax = fig.add_subplot(gs[1, 0])
S.panel(ax, "C", dx=-0.205, dy=1.18)
S.headline(ax, "Not an artifact of where we cut", y=1.045)
ths = sorted(c.threshold.unique())
off = {"PanelApp": -0.11, "HPO": 0.0, "ClinVar": 0.11}
for name, col in S.TRUTH.items():
    q = c[c.source == name].sort_values("threshold")
    ax.errorbar(q.threshold + off[name], q.odds_ratio,
                yerr=[q.odds_ratio - q.ci_low, q.ci_high - q.odds_ratio],
                fmt="o", color=col, ecolor=col, elinewidth=1.0, capsize=1.8,
                mfc="white", mew=1.2, ms=3.4, label=name)
ax.axhline(1, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.set_yscale("log")
ax.set_yticks([1, 2, 5, 10, 20, 40])
ax.set_yticklabels(["1", "2", "5", "10", "20", "40"])
ax.minorticks_off()
ax.set_xticks(ths)
ax.set_xlim(1.82, 4.18)
ax.set_ylim(0.32, 120)
ax.set_xlabel("Z threshold defining action at a site")
ax.set_ylabel("Odds ratio for\nMendelian skeletal disease")
n_lab = c.drop_duplicates("threshold").set_index("threshold").n_core
for t in ths:
    ax.text(t, 0.38, f"n = {n_lab[t]}", ha="center", va="bottom",
            fontsize=S.SMALL, color=S.MUTE)
ax.legend(loc="upper left", ncol=3, columnspacing=0.8, handletextpad=0.2,
          bbox_to_anchor=(-0.015, 1.05), borderpad=0.0)
S.ygrid(ax)

# --------------------------------------------------------- D: locus level
ax = fig.add_subplot(gs[1, 1])
S.panel(ax, "D", dx=-0.128, dy=1.18)
S.headline(ax, "Not linkage: the signal survives collapsing genes into loci",
           y=1.045)
x = np.arange(len(d))
ax.errorbar(x, d.odds_ratio,
            yerr=[d.odds_ratio - d.ci_low, d.ci_high - d.odds_ratio],
            fmt="o", color=S.SHARED, ecolor=S.SHARED, elinewidth=1.1,
            capsize=2.2, mfc="white", mew=1.3, ms=4.0)
ax.axhline(1, color=S.REST, lw=0.8, ls=(0, (3, 2)))
ax.set_yscale("log")
ax.set_yticks([1, 2, 5, 10, 20, 50, 100])
ax.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"])
ax.minorticks_off()
ax.set_xticks(x)
ax.set_xticklabels([f"{k} kb" if k else "gene\nlevel" for k in d.merge_distance_kb])
ax.set_xlabel("Distance at which neighboring genes are merged into one locus")
ax.set_ylabel("Odds ratio for\nMendelian skeletal disease")
ax.set_ylim(0.42, 800)
ax.set_xlim(-0.5, len(d) - 0.5)
for xi, r in zip(x, d.itertuples()):
    ax.text(xi, 0.50, f"{r.n_core_loci} loci", ha="center", va="bottom",
            fontsize=S.SMALL, color=S.MUTE)
S.ygrid(ax)

S.save(fig, "figure1")
