"""Figure 2 - Breadth of action predicts breadth of phenotype, and stops at the organ boundary.

Reads   results/fig2a_symmetry_law.csv, fig2a_symmetry_stats.csv,
        fig2a_symmetry_regression.csv, fig2b_cross_organ.csv, fig2c_trait_coherence.csv
Writes  figures/figure2.pdf, figures/figure2.tif

Panels
  A  breadth of action against the number of body regions the Mendelian disease affects
  B  the 2x2 cross-organ control: skeleton and brain, breadth and phenotype breadth
  C  why the brain gives no answer: its measurement sites do not share a genetic layer
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import style as S

S.use_paper_style()

a = pd.read_csv(S.RESULTS / "fig2a_symmetry_law.csv")
ast = pd.read_csv(S.RESULTS / "fig2a_symmetry_stats.csv").set_index("statistic").value
areg = pd.read_csv(S.RESULTS / "fig2a_symmetry_regression.csv").iloc[0]
x2 = pd.read_csv(S.RESULTS / "fig2b_cross_organ.csv")
coh = pd.read_csv(S.RESULTS / "fig2c_trait_coherence.csv")

fig = plt.figure(figsize=(S.W2, 2.62))
gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 0.96, 0.96], wspace=0.56,
                      left=0.086, right=0.988, top=0.795, bottom=0.295)

# ------------------------------------------------- A: the symmetry law
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "A", dx=-0.190, dy=1.24)
S.headline(ax, "Genes that act widely cause disease widely", y=1.055, x=-0.130)

# The fitted trend is the regression the text reports: the adjusted slope of
# 0.078 regions per site, drawn through the sample means. It is a line, not a
# connector between two group averages, so that what the eye measures off the
# panel is what the model estimated.
w = a.n_genes
x_bar = float((w * a.breadth).sum() / w.sum())
y_bar = float((w * a.mean_regions).sum() / w.sum())
xs = np.array([0.0, 6.0])
ax.plot(xs, y_bar + areg.beta * (xs - x_bar), color=S.SHARED_L, lw=2.6,
        solid_capstyle="round", zorder=2)

ax.errorbar(a.breadth, a.mean_regions, yerr=a["sem"], fmt="o-", color=S.SHARED,
            ecolor=S.SHARED, elinewidth=1.0, capsize=2.0, mfc="white", mew=1.3,
            ms=3.6, zorder=3)
ax.set_xticks(range(7))
ax.set_xlim(-0.35, 6.35)
ax.set_xlabel("Breadth of action (number of skeletal sites)")
ax.set_ylabel("Body regions affected by the\nMendelian disease (of 6)")
ax.set_ylim(2.35, 4.32)
ax.set_yticks([2.5, 3.0, 3.5, 4.0])
S.ygrid(ax)

ax.text(0.03, 0.975,
        f"{ast['mean_regions_breadth_le1']:.2f} → "
        f"{ast['mean_regions_breadth_ge5']:.2f} regions\n"
        f"generalized disease {ast['frac_generalised_le1'] * 100:.0f}% → "
        f"{ast['frac_generalised_ge5'] * 100:.0f}%",
        transform=ax.transAxes, fontsize=S.BASE - 0.8, va="top", ha="left",
        fontweight="bold", color=S.SHARED, linespacing=1.5)

# The statistics belong under the panel, not inside it: inside, they collided
# with the error bar at breadth 2 and ran off the right edge.
rho = float(ast["spearman_breadth_vs_n_hpo_terms"])
ax.text(0.0, -0.290,
        f"n = {int(ast['n_genes']):,} genes    β = {areg.beta:.3f} per site    "
        f"{S.pstar(areg.p)}\n"
        f"research intensity ruled out: ρ(breadth, HPO terms) = "
        f"{S.signed(rho, 3)}",
        transform=ax.transAxes, fontsize=S.SMALL, va="top", ha="left",
        color=S.MUTE, linespacing=1.55)

# --------------------------------------------- B: cross-organ 2x2 control
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "B", dx=-0.300, dy=1.24)
S.headline(ax, "The rule is skeletal, not general", y=1.055, x=-0.220)
preds = ["Skeletal breadth", "Brain breadth"]
outs = ["Skeletal phenotype breadth", "Brain phenotype breadth"]
M = np.zeros((2, 2))
for _, r in x2.iterrows():
    M[preds.index(r.predictor), outs.index(r.outcome)] = r.t

# One hue for the whole paper: the shading is the same teal as every other
# panel, so a dark cell reads as "more of the same thing", not as a new colour.
ramp = LinearSegmentedColormap.from_list("shared", ["#FFFFFF", S.SHARED_L, S.SHARED])
im = ax.imshow(M, cmap=ramp, vmin=0, vmax=3.9, aspect="auto")
for i in range(2):
    for j in range(2):
        r = x2[(x2.predictor == preds[i]) & (x2.outcome == outs[j])].iloc[0]
        strong = r.p < 0.01
        ax.text(j, i - 0.14, f"β = {r.beta:+.3f}".replace("+", "+"),
                ha="center", va="center", fontsize=S.BASE - 0.5,
                fontweight="bold", color="white" if strong else S.INK)
        ax.text(j, i + 0.17, S.pstar(r.p), ha="center", va="center",
                fontsize=S.SMALL, color="white" if strong else S.MUTE)
ax.set_xticks([0, 1])
ax.set_xticklabels(["Skeletal\nphenotype", "Brain\nphenotype"])
ax.set_yticks([0, 1])
ax.set_yticklabels(["Skeletal\nbreadth", "Brain\nbreadth"])
ax.set_xlabel("Breadth of the Mendelian phenotype")
ax.set_ylabel("Breadth of GWAS action")
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)

# --------------------------------------- C: why the brain cannot answer
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "C", dx=-0.395, dy=1.24)
S.headline(ax, "A shared layer needs\ncoherent sites", y=1.055, x=-0.300)
cols = [S.SHARED, S.BRAIN]
for i, (_, r) in enumerate(coh.iterrows()):
    ax.bar(i, r.mean_r, width=0.52, color=cols[i], zorder=2)
    ax.vlines(i, r.min_r, r.max_r, color="white", lw=1.4, zorder=3)
    ax.text(i, r.mean_r + 0.025, f"r = {r.mean_r:.2f}", ha="center",
            fontsize=S.BASE - 0.5, fontweight="bold", color=cols[i])
ax.set_xticks([0, 1])
ax.set_xticklabels(["Skeleton\n6 DXA sites", "Brain\n7 volumes"])
ax.set_xlim(-0.62, 1.62)
ax.set_ylabel("Mean genetic correlation\nbetween measurement sites")
ax.set_ylim(0, 0.60)
S.ygrid(ax)
ax.text(0, -0.315, "103 genes act\nat every site", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=S.SMALL, color=S.SHARED, fontweight="bold")
ax.text(1, -0.315, "no gene acts\nat every volume", transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=S.SMALL, color=S.BRAIN, fontweight="bold")
S.save(fig, "figure2")
