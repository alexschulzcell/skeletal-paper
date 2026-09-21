"""Vector anatomy for the schematic panels: one body outline, drawn from code.

Two panels need a human figure — Figure 1A, where the six measurement sites are
named, and the graphical abstract, where the same figure appears twice. Both
draw it from here, so the two never drift apart.

The outline is a closed cubic Bezier curve through a fixed list of landmarks.
The landmarks describe the right half of an anterior view in a coordinate frame
where the feet are at y = 0, the crown at y = 100 and the midline at x = 0; the
left half is the mirror image. Nothing here is traced from an image and nothing
is generated: the numbers below are the drawing.
"""
from __future__ import annotations

import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# --------------------------------------------------------------- landmarks
# Right half of an anterior view, crown to crotch, in the order the pen travels.
# Eight-head canon: crown 100, chin 87.5, shoulder 82, crotch 46, knee 26, sole 0.
_HALF = [
    (0.0, 100.0), (2.9, 99.5),  # crown
    (5.0, 97.8), (6.5, 93.0),   # skull
    (6.0, 88.6), (4.0, 86.6),   # jaw, chin
    (3.3, 83.8),                # neck
    (7.8, 81.6), (12.1, 78.6),  # trapezius, deltoid
    (14.0, 74.8), (13.8, 67.3), # upper arm
    (13.1, 59.8),               # elbow
    (12.7, 52.3), (12.1, 45.3), # forearm, wrist
    (12.6, 41.0), (10.8, 39.6), # hand
    (9.8, 43.5), (9.9, 49.0),   # hand inner, forearm inner
    (9.7, 57.0), (9.4, 64.0),   # inner elbow
    (8.7, 72.2), (8.0, 76.8),   # axilla
    (6.6, 71.0), (6.2, 64.5),   # flank, waist
    (7.9, 58.5), (9.0, 53.0),   # iliac crest, greater trochanter
    (8.5, 47.0), (7.5, 37.0),   # thigh
    (6.5, 27.5), (5.6, 18.0),   # knee, calf
    (4.4, 6.5), (5.7, 1.6),     # ankle, heel
    (5.2, 0.4), (1.1, 0.4),     # sole
    (1.3, 7.5), (1.7, 19.0),    # medial ankle
    (2.1, 31.0), (1.3, 43.0),   # medial thigh
    (0.0, 46.5),                # crotch, on the midline
]

# --------------------------------------------------- the six skeletal sites
# Where each DXA site sits on the outline above. The x sign picks the side the
# marker is drawn on; 0 means on the midline.
SITES = {
    "Head":         (0.0, 93.5),
    "Arms":         (-11.4, 57.5),
    "Lumbar spine": (0.0, 62.5),
    "Pelvis":       (0.0, 51.5),
    "Femoral neck": (4.9, 48.5),
    "Legs":         (-5.6, 24.0),
}
SITE_ORDER = ["Head", "Arms", "Lumbar spine", "Pelvis", "Femoral neck", "Legs"]


def _catmull_rom(points: np.ndarray) -> Path:
    """Closed smooth curve through `points`, as cubic Beziers."""
    n = len(points)
    verts = [tuple(points[0])]
    codes = [Path.MOVETO]
    for i in range(n):
        p0 = points[(i - 1) % n]
        p1 = points[i]
        p2 = points[(i + 1) % n]
        p3 = points[(i + 2) % n]
        verts += [tuple(p1 + (p2 - p0) / 6.0),
                  tuple(p2 - (p3 - p1) / 6.0),
                  tuple(p2)]
        codes += [Path.CURVE4, Path.CURVE4, Path.CURVE4]
    return Path(verts, codes)


def body_path() -> Path:
    """The closed outline of the whole figure, both halves."""
    right = np.array(_HALF, dtype=float)
    # Mirror, dropping the two midline points so they are not repeated.
    left = right[1:-1][::-1] * np.array([-1.0, 1.0])
    return _catmull_rom(np.vstack([right, left]))


def draw_body(ax, x=0.0, y=0.0, scale=1.0, face="#EDF1F3", edge="#C3CCD1",
              lw=0.8, zorder=1):
    """Place the outline with its feet at (x, y) and its crown at y + 100*scale."""
    from matplotlib.transforms import Affine2D

    path = body_path().transformed(Affine2D().scale(scale).translate(x, y))
    patch = PathPatch(path, facecolor=face, edgecolor=edge, lw=lw,
                      joinstyle="round", zorder=zorder)
    ax.add_patch(patch)
    return patch


def site_xy(name: str, x=0.0, y=0.0, scale=1.0) -> tuple[float, float]:
    """Where a named site lands once the body has been placed."""
    sx, sy = SITES[name]
    return x + sx * scale, y + sy * scale


def draw_spine(ax, x=0.0, y=0.0, scale=1.0, color="#C3CCD1", lw=0.7, zorder=2):
    """A light axial line, so the midline markers read as skeletal, not random."""
    ax.plot([x, x], [y + 48.0 * scale, y + 84.0 * scale],
            color=color, lw=lw, zorder=zorder, solid_capstyle="round")
