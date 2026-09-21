"""Graphical abstract - drawn entirely from code, no photographic or generated imagery.

Writes  figures/graphical_abstract.tif, figures/graphical_abstract.pdf

Cell Press specification followed here
  * one single panel, 5.5 inches square at 300 dpi
  * Arial, 12-16 points throughout
  * reads top to bottom, with a clear start and a clear end
  * no data items of any kind; schematic only
  * general-purpose image generators are not permitted, so every element below
    is a matplotlib primitive

The figure uses one drawing for the human body, the Bezier outline in
`anatomy.py`, at four sizes. Repeating one shape rather than sketching a
different stick figure for each box is what makes the three outcomes read as
three states of the same thing.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Polygon

import anatomy as A
import style as S

S.use_paper_style()

TEAL, OCHRE, RED = S.SHARED, S.SITE, S.HIGHLIGHT
GREY_FILL, GREY_EDGE = "#EDF1F3", "#C3CCD1"
DIM_FILL, DIM_EDGE = "#F2F4F5", "#D8DEE1"
INK, MUTE = S.INK, S.MUTE

TITLE, LEAD, BODY, TAG = 16.0, 13.5, 12.0, 12.0   # all within the 12-16 pt rule

fig = plt.figure(figsize=(5.5, 5.5), dpi=300)
ax = fig.add_axes((0, 0, 1, 1))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_axis_off()
ax.set_facecolor("white")


def body(cx, feet, scale, lit, accent, fill=GREY_FILL, edge=GREY_EDGE,
         lw=1.0, dot=2.9, fracture=False):
    """The shared outline, with the named sites marked.

    Every decoration is sized in body units and then scaled with the figure,
    so a marker keeps the same proportion to the skull whether the body is
    40 mm tall or 12 mm tall.
    """
    A.draw_body(ax, x=cx, y=feet, scale=scale, face=fill, edge=edge, lw=lw)
    for name in lit:
        x, y = A.site_xy(name, x=cx, y=feet, scale=scale)
        ax.add_patch(Circle((x, y), dot * scale, fc=accent, ec="white",
                            lw=0.6 * max(scale / 0.29, 0.5), zorder=5))
    if fracture:
        # An impact burst on the femoral neck: where the fracture this
        # genetics predicts actually happens. At icon size an anatomical break
        # would be a smudge, and a ring round the hip reads as a prohibition
        # sign, so the mark is a burst instead.
        fx, fy = A.site_xy("Femoral neck", x=cx, y=feet, scale=scale)
        for k in range(8):
            a = np.pi / 8 + k * np.pi / 4
            r0, r1 = 3.4 * scale, (9.6 if k % 2 == 0 else 7.2) * scale
            ax.plot([fx + np.cos(a) * r0, fx + np.cos(a) * r1],
                    [fy + np.sin(a) * r0, fy + np.sin(a) * r1],
                    color=accent, lw=1.1, solid_capstyle="round", zorder=6)
        ax.add_patch(Circle((fx, fy), 2.6 * scale, fc=accent, ec="white",
                            lw=0.5, zorder=7))


def arrow(x0, y0, x1, y1, color, lw=1.8, head=2.6):
    """A plain shaft with a solid triangular head, pointing at (x1, y1)."""
    v = np.array([x1 - x0, y1 - y0], dtype=float)
    n = np.hypot(*v)
    u = v / n
    p = np.array([-u[1], u[0]])
    base = np.array([x1, y1]) - u * head
    ax.plot([x0, base[0]], [y0, base[1]], color=color, lw=lw,
            solid_capstyle="butt", zorder=3)
    ax.add_patch(Polygon([(x1, y1), tuple(base + p * head * 0.52),
                          tuple(base - p * head * 0.52)],
                         closed=True, fc=color, ec="none", zorder=3))


def card(x0, x1, y0, y1, color, alpha=0.09):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=2.4",
                                fc=color, ec="none", alpha=alpha, zorder=0))


# ----------------------------------------------------------------- headline
ax.text(50, 98.4, "Breadth of action across the skeleton", ha="center",
        va="top", fontsize=TITLE, fontweight="bold", color=INK)
ax.text(50, 92.6, "One gene property, two clinical ends", ha="center",
        va="top", fontsize=LEAD, fontweight="bold", color=TEAL)

# The rule separates the finding on the left from its control on the right.
ax.plot([66.8, 66.8], [4.0, 85.0], color="#E2E7E9", lw=1.1, zorder=0)

# ------------------------------------------------- the two kinds of gene
body(33.0, 59.0, 0.275, A.SITE_ORDER, TEAL)
ax.text(33.0, 56.0, "Acts at every\nskeletal site", ha="center", va="top",
        fontsize=LEAD, fontweight="bold", color=TEAL, linespacing=1.28)

body(84.0, 60.0, 0.245, ["Femoral neck"], OCHRE, fill=DIM_FILL, edge=DIM_EDGE,
     lw=0.9)
ax.text(84.0, 56.0, "Acts at\none site", ha="center", va="top",
        fontsize=LEAD, fontweight="bold", color=OCHRE, linespacing=1.28)

# ------------------------------------------------------- down to the ends
# One stem forks to the two clinical ends; the control drops straight down.
ax.plot([33.0, 33.0], [46.0, 44.2], color=TEAL, lw=1.8, solid_capstyle="round")
ax.plot([17.0, 49.0], [44.2, 44.2], color=TEAL, lw=1.8, solid_capstyle="round")
arrow(17.0, 44.2, 17.0, 41.0, TEAL)
arrow(49.0, 44.2, 49.0, 41.0, TEAL)
arrow(84.0, 46.0, 84.0, 41.0, "#B9C2C7", lw=1.5, head=2.3)

# ------------------------------------------------------------ the outcomes
CARD_TOP, CARD_BOT = 39.8, 3.2

# The whole body is tinted where the disease is generalized, and only the hip
# is marked where it is not: the two icons differ the way the diseases do.
card(3.2, 30.8, CARD_BOT, CARD_TOP, TEAL)
ax.text(17.0, 38.4, "severe mutation", ha="center", va="top", fontsize=TAG,
        color=MUTE)
body(17.0, 18.2, 0.155, [], RED, fill="#EBC8C2", edge="#C98B80", lw=0.9)
ax.text(17.0, 16.6, "Generalized\ndysplasia\nin childhood", ha="center",
        va="top", fontsize=BODY, color=INK, linespacing=1.34)

card(35.2, 62.8, CARD_BOT, CARD_TOP, TEAL)
ax.text(49.0, 38.4, "common variation", ha="center", va="top", fontsize=TAG,
        color=MUTE)
body(49.0, 18.2, 0.155, [], RED, lw=0.9, fracture=True)
ax.text(49.0, 16.6, "Low bone density\nand fracture\nin later life",
        ha="center", va="top", fontsize=BODY, color=INK, linespacing=1.34)

card(70.0, 97.0, CARD_BOT, CARD_TOP, "#9AA3A8", alpha=0.13)
ax.text(83.5, 33.2, "Neither", ha="center", va="center", fontsize=LEAD,
        fontweight="bold", color=MUTE)
ax.text(83.5, 26.6, "no Mendelian\ndisease,\nno fracture risk", ha="center",
        va="top", fontsize=BODY, color=MUTE, linespacing=1.34)

S.check_overflow(fig, "graphical_abstract")
S.check_text_collisions(fig, "graphical_abstract")

for ext in ("pdf", "tif"):
    out = S.FIGURES / f"graphical_abstract.{ext}"
    kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}
    fig.savefig(out, dpi=300, facecolor="white", **kw)
    if ext == "tif":
        S._flatten_tiff(out, 300)
    print(f"  wrote {out.relative_to(S.ROOT)}")
plt.close(fig)
