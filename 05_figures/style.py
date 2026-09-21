"""Single source of truth for figure style: colors, fonts, sizes, helpers.

Figure requirements encoded here:
  * column widths 85 mm (1 col), 114 mm (1.5 col), 174 mm (full width)
  * sans-serif face (Arial / Helvetica), minimum 6 pt at final size
  * panel labels: uppercase, bold, upper-left corner
  * TIFF, LZW, >= 600 dpi, no alpha channel
Every figure script imports from this module and nothing else style-related.

Why the canvas is never trimmed
-------------------------------
A tight bounding box crops to the ink, so the saved width drifts away from the
figure width and differs from figure to figure. Rescaling the result to the
column width afterwards silently rescales the type with it: a 6 pt label in a
canvas that came out 184 mm wide lands at 5.7 pt on a 174 mm column, below the
journal minimum. This module therefore fixes the canvas at exactly the column
width, reserves the margins explicitly, and saves untrimmed. What the script
asks for is what the journal gets.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# --------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

# ------------------------------------------------------------------- widths
MM = 1 / 25.4
W1, W15, W2 = 85 * MM, 114 * MM, 174 * MM  # inches

# Resolution of the submitted raster: 600 dpi, the usual floor for combination
# line-and-halftone artwork; these panels are line art throughout.
DPI_TIF = 600

# ------------------------------------------------------------------- colors
# Two layers carry the whole story, so they get the only two saturated colors.
SHARED = "#1F6F8B"        # shared layer  -- deep teal
SHARED_L = "#A8CEDB"
SITE = "#D97A34"          # site-specific layer -- warm ochre
SITE_L = "#F0C9A4"
REST = "#9AA3A8"          # background genes -- neutral grey
REST_L = "#D7DCDE"
INK = "#1A1E21"           # text
MUTE = "#6E787D"          # secondary text, axes
GRID = "#E6EAEC"
HIGHLIGHT = "#B03A2E"     # the one number that must be seen
BRAIN = "#7D6B9E"         # cross-organ control -- muted violet

# Three truth sides, one hue family so they read as one measurement.
TRUTH = {"PanelApp": "#0F4C5C", "HPO": "#2E8B9A", "ClinVar": "#6FB3BE"}
# Three replication cohort sites.
GEFOS = {"Femoral neck": "#1F6F8B", "Lumbar spine": "#4E9EB8", "Forearm": "#9AA3A8"}

# --------------------------------------------------------------------- rc
BASE = 7.0
SMALL = 6.0   # minimum in-figure font size


def use_paper_style() -> None:
    """Apply the manuscript rc settings. Call once at the top of a figure script."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        # Exponents are set with mathtext. Left at the default, they come out in
        # DejaVu Sans and the figure ships two typefaces; the journal asks for
        # one. Point every mathtext family at the body face instead.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.default": "regular",
        "font.size": BASE,
        "axes.titlesize": BASE + 0.5,
        "axes.labelsize": BASE,
        "xtick.labelsize": BASE - 0.5,
        "ytick.labelsize": BASE - 0.5,
        "legend.fontsize": BASE - 0.5,
        "axes.titleweight": "bold",
        "axes.labelcolor": INK,
        "text.color": INK,
        "axes.edgecolor": MUTE,
        "axes.linewidth": 0.6,
        "xtick.color": MUTE,
        "ytick.color": MUTE,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.unicode_minus": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        "figure.dpi": 150,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        # Never trim: see the module docstring.
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
        "savefig.dpi": DPI_TIF,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# ----------------------------------------------------------------- helpers
def panel(ax, letter: str, dx: float = -0.16, dy: float = 1.10) -> None:
    """Bold uppercase panel label in the upper-left corner."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=BASE + 2,
            fontweight="bold", va="top", ha="left", color=INK,
            clip_on=False, zorder=10)


def headline(ax, text: str, y: float = 1.02, x: float = 0.0,
             color: str | None = None) -> None:
    """One short, quantitative take-home line directly above a panel."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=BASE - 0.5,
            fontweight="bold", va="bottom", ha="left", color=color or INK,
            clip_on=False, zorder=10)


def ygrid(ax) -> None:
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.5)


def pstar(p: float) -> str:
    """p-value as it should appear in a panel: compact, scientific, unambiguous."""
    if p >= 0.01:
        return f"P = {p:.2f}"
    exp = int(f"{p:.0e}".split("e")[1])
    man = f"{p:.0e}".split("e")[0]
    return f"P = {man}×10$^{{{exp}}}$"


def signed(value: float, nd: int = 3) -> str:
    """A signed number with a real minus sign, not a hyphen."""
    return f"{value:+.{nd}f}".replace("-", "−")


def trim_spine(ax, axis: str = "x") -> None:
    """Shorten a spine to the span of the ticks.

    An axis whose limits are widened to make room for annotation would
    otherwise draw a rule out into empty space.
    """
    ticks = ax.get_xticks() if axis == "x" else ax.get_yticks()
    lo, hi = (ax.get_xlim() if axis == "x" else ax.get_ylim())
    inside = [t for t in ticks if lo <= t <= hi]
    if not inside:
        return
    spine = ax.spines["bottom" if axis == "x" else "left"]
    spine.set_bounds(min(inside), max(inside))


def fit_aspect(ax, ylim: tuple[float, float], x_center: float = 0.0) -> None:
    """Choose x limits so that one data unit is the same length on both axes.

    Used instead of ``set_aspect('equal')`` where the drawing must not be
    distorted but the panel box must keep the size the layout gave it.
    """
    fig = ax.figure
    box = ax.get_position()
    w_in = box.width * fig.get_figwidth()
    h_in = box.height * fig.get_figheight()
    span_x = (ylim[1] - ylim[0]) * w_in / h_in
    ax.set_ylim(*ylim)
    ax.set_xlim(x_center - span_x / 2, x_center + span_x / 2)


def spread_labels(ys: list[float], min_gap: float,
                  lo: float, hi: float) -> list[float]:
    """Push overlapping label anchors apart, keeping their order and centre.

    `ys` are the ideal positions in display or data units, `min_gap` the
    smallest separation two labels may have. The result stays inside
    [lo, hi]. This replaces per-label nudge factors, which have to be
    re-tuned by hand whenever a number changes.
    """
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = list(ys)
    # Forward pass: enforce the gap going up.
    for k in range(1, len(order)):
        i, j = order[k - 1], order[k]
        if out[j] - out[i] < min_gap:
            out[j] = out[i] + min_gap
    # Backward pass: pull the stack down if it overflowed the top.
    over = out[order[-1]] - hi
    if over > 0:
        for i in order:
            out[i] -= over
    under = lo - out[order[0]]
    if under > 0:
        for i in order:
            out[i] += under
    return out


def place_labels(ax, xs, ys, texts, colors, fontsize=None, pad_pt=1.8,
                 radii_pt=(4.5, 7.0, 10.0, 14.0, 19.0, 25.0), marker_pt=2.6,
                 obstacles=None, leader_color=None, verbose=False, halo=True,
                 sweeps=4):
    """Label a scatter without collisions, by search rather than by hand.

    Each label is scored over sixteen directions at six distances. The cost
    adds, in descending order of weight: running off the axes, overlapping
    another label, sitting on somebody else's labelled anchor, a leader that
    crosses another leader or passes behind another label, covering ordinary
    markers, and finally distance from its own point.

    A single greedy pass makes bad early choices — the first label placed does
    not know where the others will go — so the greedy result is followed by
    refinement sweeps in which every label is re-optimised against the final
    positions of all the others. A white halo lets a name rest against its own
    marker without the point cloud showing through, which keeps leaders short
    and the assignment obvious.

    `obstacles` is the full (x, y) array of the panel's markers; without it
    only the labelled points are avoided. This replaces a dictionary of
    hand-tuned offsets, which goes stale the moment a value changes.
    """
    import numpy as np
    from matplotlib.transforms import Bbox

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fs = fontsize if fontsize is not None else SMALL
    px_pt = fig.dpi / 72

    pts = ax.transData.transform(np.column_stack([xs, ys]))
    ob = pts if obstacles is None else ax.transData.transform(np.asarray(obstacles))
    axes_box = ax.get_window_extent(renderer)

    m = marker_pt * px_pt
    blocked = [Bbox.from_bounds(x - m, y - m, 2 * m, 2 * m) for x, y in ob]
    ma = (marker_pt + 3.0) * px_pt
    anchors = [Bbox.from_bounds(x - ma, y - ma, 2 * ma, 2 * ma) for x, y in pts]

    dirs = [(np.cos(a), np.sin(a))
            for a in np.linspace(0, 2 * np.pi, 16, endpoint=False)]

    # --- measure every label once -------------------------------------------
    arts, sizes = [], []
    for i, txt in enumerate(texts):
        t = ax.text(0, 0, txt, fontsize=fs, fontweight="bold", color=colors[i],
                    ha="center", va="center", transform=ax.transData,
                    clip_on=False, zorder=6)
        if halo:
            t.set_bbox(dict(boxstyle="round,pad=0.10", facecolor="white",
                            edgecolor="none", alpha=0.82))
        tb = t.get_window_extent(renderer)
        arts.append(t)
        sizes.append((tb.width + 2 * pad_pt * px_pt,
                      tb.height + 2 * pad_pt * px_pt))

    def overlap(a, b):
        dx = min(a.x1, b.x1) - max(a.x0, b.x0)
        dy = min(a.y1, b.y1) - max(a.y0, b.y0)
        return dx * dy if dx > 0 and dy > 0 else 0.0

    def crosses(p1, p2, p3, p4):
        def side(a, b, c):
            return ((b[0] - a[0]) * (c[1] - a[1])
                    - (b[1] - a[1]) * (c[0] - a[0]))
        d1, d2 = side(p3, p4, p1), side(p3, p4, p2)
        d3, d4 = side(p1, p2, p3), side(p1, p2, p4)
        return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))

    def score(i, cand, others):
        """Cost of putting label i at `cand`, given the others' positions."""
        (cx, cy), r_pt = cand
        w, h = sizes[i]
        unit = w * h
        px, py = pts[i]
        box = Bbox.from_bounds(cx - w / 2, cy - h / 2, w, h)
        cost = 0.0
        if (box.x0 < axes_box.x0 or box.x1 > axes_box.x1
                or box.y0 < axes_box.y0 or box.y1 > axes_box.y1):
            cost += 1e4
        for j, other in others.items():
            if j == i or other is None:
                continue
            (ox, oy), _ = other
            ow, oh = sizes[j]
            obox = Bbox.from_bounds(ox - ow / 2, oy - oh / 2, ow, oh)
            cost += 240.0 * overlap(box, obox) / unit
            if crosses((px, py), (cx, cy), pts[j], (ox, oy)):
                cost += 12.0
            for s in range(1, 10):
                f = s / 10.0
                if obox.contains(px + (cx - px) * f, py + (cy - py) * f):
                    cost += 14.0
                    break
        for k, b in enumerate(anchors):
            if k != i:
                cost += 14.0 * overlap(box, b) / unit
        for b in blocked:
            cost += 2.0 * overlap(box, b) / unit
        for s in range(1, 10):
            f = s / 10.0
            if any(b.contains(px + (cx - px) * f, py + (cy - py) * f)
                   for b in blocked):
                cost += 1.0
                break
        return cost + 0.95 * r_pt

    def candidates(i):
        w, h = sizes[i]
        px, py = pts[i]
        out = []
        for r_pt in radii_pt:
            r = r_pt * px_pt
            for dx, dy in dirs:
                out.append((((px + dx * (r + w / 2)),
                             (py + dy * (r + h / 2))), r_pt))
        return out

    centre = pts.mean(axis=0)
    order = sorted(range(len(texts)),
                   key=lambda i: -float(np.hypot(*(pts[i] - centre))))
    chosen: dict[int, tuple] = {i: None for i in range(len(texts))}
    for i in order:
        chosen[i] = min(candidates(i), key=lambda c: score(i, c, chosen))
    for _ in range(sweeps):
        moved = False
        for i in order:
            best = min(candidates(i), key=lambda c: score(i, c, chosen))
            if best != chosen[i]:
                chosen[i], moved = best, True
        if not moved:
            break

    # --- commit ------------------------------------------------------------
    inv = ax.transData.inverted()
    report = []
    for i, t in enumerate(arts):
        (cx, cy), r_pt = chosen[i]
        px, py = pts[i]
        w, h = sizes[i]
        t.set_position(tuple(inv.transform((cx, cy))))
        report.append((texts[i], round(score(i, chosen[i], chosen), 2), r_pt))
        if r_pt > radii_pt[0]:
            vx, vy = cx - px, cy - py
            n = max(float(np.hypot(vx, vy)), 1e-9)
            sx = px + vx / n * (marker_pt + 0.8) * px_pt
            sy = py + vy / n * (marker_pt + 0.8) * px_pt
            ex = cx - vx / n * (w / 2 + 1.0 * px_pt)
            ey = cy - vy / n * (h / 2 + 1.0 * px_pt)
            if float(np.hypot(ex - sx, ey - sy)) > 1.0:
                x0, y0 = inv.transform((sx, sy))
                x1, y1 = inv.transform((ex, ey))
                ax.plot([x0, x1], [y0, y1], color=leader_color or REST,
                        lw=0.55, zorder=2, solid_capstyle="round")

    worst = max(report, key=lambda r: r[1])
    if verbose or worst[1] > 5.0:
        print(f"    label placement: worst residual {worst[0]} cost {worst[1]}")
    return report


def check_overflow(fig, name: str, tol_mm: float = 0.25) -> None:
    """Fail if any ink falls outside the canvas.

    Because the canvas is never trimmed, a label that runs past the edge is
    silently cut off in the saved file instead of enlarging it. This compares
    the tight bounding box of everything drawn against the figure rectangle and
    names the offending artists, so the failure is loud and local.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    w_in, h_in = fig.get_size_inches()
    tol = tol_mm / 25.4
    bad = []
    # set_axis_off() stops the frame and the ticks from being drawn but leaves
    # their artists in the tree, still reporting an extent.
    hidden = set()
    for a in fig.axes:
        if not a.axison:
            for s in a.spines.values():
                hidden.add(id(s))
            for axis in (a.xaxis, a.yaxis):
                hidden.add(id(axis))
                hidden.add(id(axis.label))
                for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                    hidden.update({id(tick), id(tick.label1), id(tick.label2)})
    for art in fig.get_children():
        for a in [art, *getattr(art, "get_children", list)()]:
            if not a.get_visible() or id(a) in hidden:
                continue
            try:
                bb = a.get_window_extent(renderer)
            except Exception:
                continue
            if bb.width <= 0 or bb.height <= 0:
                continue
            # An artist that is clipped to its axes cannot print outside it, so
            # judge it on the part that actually reaches the page. Bars drawn
            # from a baseline below the visible y range are the common case.
            if a.get_clip_on():
                clip = a.get_clip_box()
                if clip is not None:
                    from matplotlib.transforms import Bbox
                    bb = Bbox.intersection(bb, clip) or bb
                    if bb.width <= 0 or bb.height <= 0:
                        continue
            bb = bb.transformed(fig.dpi_scale_trans.inverted())
            over = max(-bb.x0, -bb.y0, bb.x1 - w_in, bb.y1 - h_in)
            if over > tol:
                label = getattr(a, "get_text", lambda: "")() or type(a).__name__
                bad.append((over * 25.4, str(label)[:46]))
    if bad:
        bad.sort(reverse=True)
        lines = "\n".join(f"      {o:5.2f} mm  {t!r}" for o, t in bad[:8])
        raise ValueError(
            f"{name}: {len(bad)} artist(s) fall outside the {w_in * 25.4:.0f} x "
            f"{h_in * 25.4:.0f} mm canvas and would be cut off:\n{lines}"
        )


def check_text_collisions(fig, name: str, min_overlap_mm2: float = 0.35) -> None:
    """Fail if two pieces of text are drawn on top of each other.

    Overlapping labels are the defect a reader notices first and the one a
    layout is most likely to reintroduce when a number changes. Checking for
    it here means the figure cannot be written while it is present.
    """
    from matplotlib.text import Text

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    items = []

    # Matplotlib keeps a label artist for every tick it generated, including
    # ticks that fall outside the current view and are never drawn. Those would
    # register as phantom collisions, so they are excluded by position.
    ghosts = set()
    for ax in fig.axes:
        off = not ax.axison          # set_axis_off() leaves the label artists behind
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            if off:
                ghosts.add(id(axis.label))
            lo, hi = sorted(lim)
            for tick in axis.get_major_ticks():
                if off or not (lo - 1e-9 <= tick.get_loc() <= hi + 1e-9):
                    ghosts.add(id(tick.label1))
                    ghosts.add(id(tick.label2))

    def walk(art):
        for a in getattr(art, "get_children", list)():
            if isinstance(a, Text):
                if a.get_visible() and a.get_text().strip() and id(a) not in ghosts:
                    try:
                        bb = a.get_window_extent(renderer)
                    except Exception:
                        continue
                    if bb.width > 0 and bb.height > 0:
                        items.append((a, bb))
            else:
                walk(a)

    walk(fig)
    scale = (25.4 / fig.dpi) ** 2
    clashes = []
    for i in range(len(items)):
        ai, bi = items[i]
        for j in range(i + 1, len(items)):
            aj, bj = items[j]
            dx = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            dy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if dx <= 0 or dy <= 0:
                continue
            area = dx * dy * scale
            if area > min_overlap_mm2:
                clashes.append((area, ai.get_text()[:34], aj.get_text()[:34]))
    if clashes:
        clashes.sort(reverse=True)
        lines = "\n".join(f"      {a:6.2f} mm2  {x!r} / {y!r}"
                          for a, x, y in clashes[:8])
        raise ValueError(
            f"{name}: {len(clashes)} pair(s) of overlapping text:\n{lines}")


def save(fig, name: str, width_mm: float = 174.0) -> None:
    """Write the PDF (vector record) and the TIFF (submission raster).

    The canvas is saved untrimmed, so the file measures exactly what the figure
    was sized for and the type keeps the point size the script asked for. The
    TIFF is flattened onto white and written without an alpha channel, which
    production systems reject, at exactly `DPI_TIF`.
    """
    got_mm = fig.get_figwidth() * 25.4
    if abs(got_mm - width_mm) > 0.05:
        raise ValueError(
            f"{name}: figure is {got_mm:.2f} mm wide but {width_mm:.2f} mm was "
            f"requested. Set figsize to the column width; do not rescale afterwards."
        )
    check_overflow(fig, name)
    check_text_collisions(fig, name)

    pdf = FIGURES / f"{name}.pdf"
    fig.savefig(pdf)
    print(f"  wrote {pdf.relative_to(ROOT)}  ({got_mm:.1f} mm wide)")

    tif = FIGURES / f"{name}.tif"
    fig.savefig(tif, dpi=DPI_TIF, pil_kwargs={"compression": "tiff_lzw"})
    _flatten_tiff(tif, DPI_TIF)
    print(f"  wrote {tif.relative_to(ROOT)}  ({width_mm:.0f} mm at {DPI_TIF} dpi, RGB)")
    plt.close(fig)


def _flatten_tiff(path: Path, dpi: int) -> None:
    """Drop the alpha channel and stamp the resolution tag.

    Matplotlib writes RGBA. Most production pipelines reject TIFFs that carry
    transparency, so the image is composited onto white and rewritten as RGB.
    """
    from PIL import Image

    with Image.open(path) as im:
        im.load()
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            flat = Image.new("RGB", im.size, "white")
            flat.paste(im, mask=im.split()[-1])
        else:
            flat = im.convert("RGB")
        flat.save(path, compression="tiff_lzw", dpi=(dpi, dpi))
