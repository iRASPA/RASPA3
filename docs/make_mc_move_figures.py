"""Draws the illustrations of the Monte Carlo moves of docs/manual/mc_moves.md.

One schematic per move, named after the anchor of its section: docs/images/mc_move_<anchor>.png. The figures are
part of the repository; run this script to draw them again after a change (it needs matplotlib and numpy):

    python docs/make_mc_move_figures.py [anchor ...]

It also puts the figure, with its caption (CAPTIONS) and the paragraph that explains it (EXPLANATIONS), under the
'Keyword:' paragraph of the section of the move in the manual, in place of what a previous run put there.

The pictures are schematic and two-dimensional. What they share:

    grey molecules         the other molecules of the box, which the move leaves alone
    faded, dashed outline  the old configuration
    full colour            the new configuration
    thin dashed outline    trial positions or trial molecules that were generated and not chosen
    hatched                a fractional molecule, coupled to the rest with the parameter lambda
"""

import math
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

OUTPUT = pathlib.Path(__file__).resolve().parent / "images"

BLUE = "#2f6db5"
ORANGE = "#e08a1e"
GREEN = "#3c9a5f"
RED = "#c8443c"
GREY = "#c9ced6"
GREY_EDGE = "#9aa2ad"
INK = "#2b2f36"
SOFT = "#6b7380"

RADIUS = 0.034
BOND = 0.066

PANEL = 2.25  # the side of a box panel [inch]
GAP = 0.62  # the space between two panels, for an arrow
SMALL = 1.05  # the width of a reservoir panel
TITLE = 0.30
BAR = 0.62  # the height of a lambda bar
MARGIN = 0.08


# ----------------------------------------------------------------------------------------------------------------
# primitives
# ----------------------------------------------------------------------------------------------------------------


def chain(x, y, angle=0.0, n=3, bend=0.0):
    """The bead centres of a molecule of n beads around (x, y); 'bend' is the angle between successive bonds."""
    points = [np.zeros(2)]
    direction = math.radians(angle)
    for _ in range(n - 1):
        points.append(points[-1] + BOND * np.array([math.cos(direction), math.sin(direction)]))
        direction += math.radians(bend)
    points = np.array(points)
    return points - points.mean(axis=0) + np.array([x, y])


STYLES = {
    # face alpha, edge alpha, line style, line width, hatch
    "solid": (1.0, 1.0, "-", 1.0, None),
    "old": (0.16, 0.75, (0, (2.2, 1.6)), 0.9, None),
    "trial": (0.0, 0.8, (0, (1.6, 1.4)), 0.75, None),
    "fractional": (0.30, 1.0, "-", 1.0, "//////"),
}


def beads(ax, points, color=BLUE, style="solid", z=4, radius=RADIUS, bonds=True):
    face_alpha, edge_alpha, line_style, line_width, hatch = STYLES[style]
    rgb = matplotlib.colors.to_rgb(color)
    edge = tuple(0.62 * c for c in rgb) if style in ("solid", "fractional") else rgb
    if bonds and len(points) > 1:
        ax.plot(points[:, 0], points[:, 1], color=edge, alpha=edge_alpha, lw=1.6 if style == "solid" else 0.9,
                ls="-" if style in ("solid", "fractional") else line_style, zorder=z, solid_capstyle="round")
    for point in points:
        ax.add_patch(Circle(point, radius, facecolor=(*rgb, face_alpha), edgecolor=(*edge, edge_alpha), lw=line_width,
                            ls=line_style, hatch=hatch, zorder=z + 1))


def molecule(ax, x, y, angle=0.0, n=3, color=BLUE, style="solid", bend=0.0, z=4):
    points = chain(x, y, angle, n, bend)
    beads(ax, points, color, style, z)
    return points


def ion(ax, x, y, sign, color, style="solid", z=4):
    beads(ax, np.array([[x, y]]), color, style, z, radius=RADIUS * 1.25)
    if style in ("solid", "fractional"):
        ax.text(x, y - 0.002, sign, ha="center", va="center", fontsize=8.5, zorder=z + 3, fontweight="bold",
                color="white" if style == "solid" else INK)


def bystanders(ax, which="ring", seed_shift=0.0):
    """The molecules the move leaves alone, around the middle of the box."""
    places = {
        "ring": [(0.14, 0.14, 25), (0.86, 0.16, 100), (0.13, 0.86, -25), (0.84, 0.87, 55), (0.50, 0.92, 5),
                 (0.91, 0.52, 85), (0.09, 0.50, 80), (0.52, 0.08, 12)],
        "sparse": [(0.16, 0.16, 25), (0.85, 0.82, 60), (0.15, 0.84, -30), (0.86, 0.18, 100)],
        "left": [(0.13, 0.14, 25), (0.12, 0.86, -25), (0.09, 0.50, 80), (0.45, 0.92, 5), (0.45, 0.08, 12)],
        "dense": [(0.14, 0.14, 25), (0.86, 0.16, 100), (0.13, 0.86, -25), (0.84, 0.87, 55), (0.50, 0.92, 5),
                  (0.91, 0.52, 85), (0.09, 0.50, 80), (0.52, 0.08, 12), (0.30, 0.68, 40), (0.70, 0.66, -40),
                  (0.32, 0.30, -60), (0.69, 0.32, 20)],
    }[which]
    for x, y, angle in places:
        molecule(ax, x, y, angle + seed_shift, color=GREY, z=2)
        # (drawn as solid grey with a softer edge)
    for patch in ax.patches:
        if matplotlib.colors.to_hex(patch.get_facecolor()[:3]) == matplotlib.colors.to_hex(GREY):
            patch.set_edgecolor(GREY_EDGE)
    for line in ax.lines:
        if np.allclose(matplotlib.colors.to_rgb(line.get_color()), tuple(0.62 * c for c in matplotlib.colors.to_rgb(GREY))):
            line.set_color(GREY_EDGE)


def arrow(ax, start, end, color=INK, rad=0.0, lw=1.2, style="-|>", ls="-", z=8, scale=9, shrink=0.0):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=scale, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z, shrinkA=shrink, shrinkB=shrink))


def arc_arrow(ax, centre, radius, start, end, color=INK, lw=1.2, z=8):
    """An arrow along the circle around 'centre', from the angle 'start' to the angle 'end' [degrees]."""
    centre = np.asarray(centre, dtype=float)
    angles = np.radians(np.linspace(start, end, 40))
    points = centre + radius * np.column_stack([np.cos(angles), np.sin(angles)])
    ax.plot(points[:-2, 0], points[:-2, 1], color=color, lw=lw, zorder=z, solid_capstyle="round")
    arrow(ax, points[-4], points[-1], color, lw=lw, z=z, scale=8)


def label(ax, x, y, text, color=INK, size=7.5, ha="center", va="center", **keywords):
    ax.text(x, y, text, color=color, fontsize=size, ha=ha, va=va, zorder=12,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.6), **keywords)


def grow(ax, start, angles, color=BLUE, style="solid", alternatives=(-98, 92), z=4):
    """A molecule grown bead by bead from 'start': at every bead the trial directions that were not chosen are
    drawn as thin dashed beads next to the chosen one."""
    points = [np.array(start, dtype=float)]
    for angle in angles:
        for offset in alternatives:
            direction = math.radians(angle + offset)
            trial = points[-1] + BOND * np.array([math.cos(direction), math.sin(direction)])
            ax.plot(*zip(points[-1], trial), color=color, lw=0.7, ls=(0, (1.6, 1.4)), alpha=0.8, zorder=z - 1)
            beads(ax, np.array([trial]), color, "trial", z - 1, radius=0.8 * RADIUS)
        direction = math.radians(angle)
        points.append(points[-1] + BOND * np.array([math.cos(direction), math.sin(direction)]))
    points = np.array(points)
    beads(ax, points, color, style, z + 1)
    return points


def trial_positions(ax, places, color, radius=RADIUS * 1.25):
    for place in places:
        beads(ax, np.array([place]), color, "trial", 3, radius=radius)


def sphere(ax, centre, radius, text=r"$R_\mathrm{max}$"):
    ax.add_patch(Circle(centre, radius, facecolor=(0, 0, 0, 0.03), edgecolor=SOFT, lw=0.8, ls=(0, (3, 2)), zorder=1))
    edge = (centre[0] + radius * math.cos(math.radians(140)), centre[1] + radius * math.sin(math.radians(140)))
    ax.plot(*zip(centre, edge), color=SOFT, lw=0.7, zorder=1)
    label(ax, edge[0] - 0.02, edge[1] + 0.045, text, SOFT, 7)


# ----------------------------------------------------------------------------------------------------------------
# the layout of a figure: panels in a row, with titles above and lambda bars below
# ----------------------------------------------------------------------------------------------------------------


class Figure:
    def __init__(self, panels, connectors=None, bars=None, note=None):
        """panels: (kind, title) with kind 'box' or 'small'; connectors: the text between successive panels;
        bars: {panel index or (first, last): True}, a lambda bar under these panels; note: a line under it all."""
        widths = [PANEL if kind == "box" else SMALL for kind, _ in panels]
        self.bar_height = BAR if bars else 0.0
        self.note_height = 0.26 if note else 0.0
        width = 2 * MARGIN + sum(widths) + GAP * (len(panels) - 1)
        height = 2 * MARGIN + TITLE + PANEL + self.bar_height + self.note_height
        self.figure = plt.figure(figsize=(width, height), dpi=220)
        self.width, self.height = width, height
        self.axes, self.lefts, self.widths = [], [], widths
        bottom = MARGIN + self.bar_height + self.note_height
        left = MARGIN
        for (kind, title), panel_width in zip(panels, widths):
            ax = self.figure.add_axes([left / width, bottom / height, panel_width / width, PANEL / height])
            ratio = panel_width / PANEL
            ax.set_xlim(0.5 - 0.5 * ratio * 1.06, 0.5 + 0.5 * ratio * 1.06)
            ax.set_ylim(-0.03, 1.03)
            ax.set_aspect("equal")
            ax.axis("off")
            if title:
                self.figure.text((left + panel_width / 2) / width, (bottom + PANEL + 0.07) / height, title,
                                 ha="center", va="bottom", fontsize=8, color=INK)
            self.axes.append(ax)
            self.lefts.append(left)
            left += panel_width + GAP
        for index, text in enumerate(connectors or []):
            if text:
                x = (self.lefts[index] + widths[index] + GAP / 2) / width
                self.figure.text(x, (bottom + PANEL / 2) / height, text, ha="center", va="center", fontsize=15,
                                 color=INK)
        self.bars = {}
        for key in bars or {}:
            first, last = key if isinstance(key, tuple) else (key, key)
            start = self.lefts[first] + 0.12
            end = self.lefts[last] + widths[last] - 0.12
            ax = self.figure.add_axes([start / width, (MARGIN + self.note_height) / height, (end - start) / width,
                                       self.bar_height / height])
            ax.set_xlim(-0.08, 1.08)
            ax.set_ylim(-1.0, 1.0)
            ax.axis("off")
            self.bars[key] = ax
        if note:
            self.figure.text(0.5, (MARGIN + 0.02) / height, note, ha="center", va="bottom", fontsize=7.2, color=SOFT)

    def save(self, name):
        OUTPUT.mkdir(exist_ok=True)
        self.figure.savefig(OUTPUT / f"mc_move_{name}.png", facecolor="white")
        plt.close(self.figure)


def box(ax, x=0.0, y=0.0, w=1.0, h=1.0, style="-", color=INK, lw=1.1, fill="#fbfbfc"):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fill if fill else "none", edgecolor=color, lw=lw, ls=style,
                           zorder=0 if fill else 1))


def reservoir(ax, text="reservoir"):
    ratio = SMALL / PANEL
    ax.add_patch(Rectangle((0.5 - 0.5 * ratio, 0.22), ratio, 0.56, facecolor="#f3f5f8", edgecolor=SOFT, lw=0.9,
                           ls=(0, (3, 2)), zorder=0))
    label(ax, 0.5, 0.16, text, SOFT, 7)


def lambda_bar(ax, value, color=BLUE, name=r"$\lambda$", left_text=None, right_text=None, walls=False, step=0.13,
               track_label=None, text=None, directions=(-1, 1)):
    """A lambda axis with the present value, the range of a trial change, and what happens at the ends."""
    ax.plot([0, 1], [0, 0], color=INK, lw=1.2, zorder=2, solid_capstyle="butt")
    for x in (0, 1):
        ax.plot([x, x], [-0.22, 0.22], color=INK, lw=2.6 if walls else 1.2, zorder=2)
    for x in np.linspace(0.1, 0.9, 9):
        ax.plot([x, x], [-0.1, 0.1], color=SOFT, lw=0.6, zorder=2)
    ax.plot([value], [0], marker="o", ms=6.5, mfc=color, mec=INK, mew=0.8, zorder=6)
    for sign in directions:
        arrow(ax, (value + sign * 0.03, 0.36), (value + sign * (0.03 + step), 0.36), color, lw=1.1, scale=7)
    label(ax, value, 0.76, text or (name + r"$\,\pm\,\Delta\lambda$"), color, 7.5)
    label(ax, 0, -0.60, "0", INK, 7.5)
    label(ax, 1, -0.60, "1", INK, 7.5)
    if track_label:
        label(ax, 0.5, -0.60, track_label, SOFT, 6.8)
    if left_text:
        label(ax, 0.05, -0.60, left_text, SOFT, 6.6, ha="left")
    if right_text:
        label(ax, 0.95, -0.60, right_text, SOFT, 6.6, ha="right")


# ----------------------------------------------------------------------------------------------------------------
# single-molecule moves
# ----------------------------------------------------------------------------------------------------------------


def translation(name, random):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    bystanders(ax)
    if random:
        ax.plot([0.02, 0.98], [0.5, 0.5], color=SOFT, lw=0.7, ls=(0, (1, 2)), zorder=1)
        molecule(ax, 0.24, 0.5, 35, style="old")
        molecule(ax, 0.72, 0.5, 35)
        arrow(ax, (0.33, 0.56), (0.63, 0.56), BLUE, rad=-0.25)
        label(ax, 0.48, 0.70, "anywhere along one axis", BLUE, 7)
    else:
        molecule(ax, 0.40, 0.5, 35, style="old")
        molecule(ax, 0.56, 0.5, 35)
        arrow(ax, (0.40, 0.36), (0.56, 0.36), BLUE)
        ax.plot([0.26, 0.74], [0.28, 0.28], color=SOFT, lw=0.8, zorder=1)
        for x in (0.26, 0.74):
            ax.plot([x, x], [0.26, 0.30], color=SOFT, lw=0.8, zorder=1)
        label(ax, 0.50, 0.215, r"$\pm\,\Delta_\mathrm{max}$", SOFT, 7.5)
    figure.save(name)


def rotation(name, random):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    bystanders(ax)
    centre = np.array([0.5, 0.5])
    start, end = (20, 250) if random else (20, 62)
    molecule(ax, *centre, start, style="old")
    molecule(ax, *centre, end)
    radius = 0.15
    if random:
        ax.add_patch(Circle(centre, radius, facecolor="none", edgecolor=SOFT, lw=0.7, ls=(0, (1, 2)), zorder=1))
        label(ax, 0.5, 0.72, "any orientation", BLUE, 7)
    arc_arrow(ax, centre, radius, start, end, BLUE)
    if not random:
        label(ax, 0.74, 0.64, r"$\pm\,\Delta\theta_\mathrm{max}$", SOFT, 7.5)
    ax.plot(*centre, marker="+", color="white", ms=5, mew=0.9, zorder=9)
    figure.save(name)


def smart(name, translate, rotate, everything=False):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    if everything:
        places = [(0.20, 0.20, 25, 40), (0.78, 0.22, 100, 150), (0.20, 0.80, -25, -50), (0.80, 0.80, 55, 215),
                  (0.48, 0.46, 10, 80)]
    else:
        bystanders(ax)
        places = [(0.42, 0.46, 20, 30)]
    for x, y, angle, force_angle in places:
        force = np.array([math.cos(math.radians(force_angle)), math.sin(math.radians(force_angle))])
        normal = np.array([-force[1], force[0]])
        length = 0.11 if everything else 0.15
        shift = force * length if translate else np.zeros(2)
        # the random part of the proposal, next to the drift along the force
        extra = np.array([0.030, -0.035]) if (translate and not everything) else np.zeros(2)
        new = np.array([x, y]) + shift + extra
        if translate and not everything:
            for scale in (0.09, 0.15):
                ax.add_patch(Circle((x + shift[0], y + shift[1]), scale, edgecolor="none", zorder=1,
                                    facecolor=(*matplotlib.colors.to_rgb(BLUE), 0.07)))
        molecule(ax, x, y, angle, style="old")
        molecule(ax, *new, angle + (32 if rotate else 0))
        if translate:
            offset = normal * (0.085 if not everything else 0.07)
            arrow(ax, np.array([x, y]) + offset, np.array([x, y]) + offset + force * length * 1.25, RED, lw=1.4, scale=8)
            if not everything:
                label(ax, *(np.array([x, y]) + offset + force * length * 1.25 + np.array([0.035, 0.03])),
                      r"$\mathbf{F}$", RED, 9)
        if rotate:
            r = 0.135
            arc_arrow(ax, new, r, angle + 195, angle + 255, RED, lw=1.4)
            if not everything:
                middle = math.radians(angle + 225)
                label(ax, *(new + (r + 0.06) * np.array([math.cos(middle), math.sin(middle)])),
                      r"$\boldsymbol{\tau}$", RED, 9)
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# configurational-bias moves
# ----------------------------------------------------------------------------------------------------------------


def reinsertion(name, multi):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    bystanders(ax, "sparse")
    if multi:
        molecule(ax, 0.24, 0.50, 70, style="old")
        for x, y, angle in ((0.52, 0.80, 10), (0.50, 0.20, -15), (0.56, 0.50, 60)):
            molecule(ax, x, y, angle, style="trial")
        molecule(ax, 0.78, 0.50, 110)
        arrow(ax, (0.30, 0.62), (0.74, 0.64), BLUE, rad=-0.3)
        label(ax, 0.50, 0.93, r"$M$ trial molecules", BLUE, 7)
    else:
        old = grow(ax, (0.20, 0.36), (75, 95), style="old")
        new = grow(ax, (0.66, 0.42), (60, 100))
        arrow(ax, (0.30, 0.52), (0.60, 0.56), BLUE, rad=-0.3)
        label(ax, 0.22, 0.26, "retrace", SOFT, 7)
        label(ax, 0.70, 0.32, "grow", BLUE, 7)
    figure.save(name)


def partial_reinsertion(name):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    bystanders(ax, "sparse")
    fixed = chain(0.36, 0.50, 0, 3)
    beads(ax, fixed, INK)
    label(ax, 0.36, 0.40, "fixed", INK, 7)
    old = [fixed[-1]]
    for angle in (60, 25):
        old.append(old[-1] + BOND * np.array([math.cos(math.radians(angle)), math.sin(math.radians(angle))]))
    beads(ax, np.array(old), BLUE, "old", 3)
    ax.add_patch(Circle(fixed[-1], RADIUS, facecolor=INK, edgecolor=INK, zorder=7))
    grow(ax, fixed[-1], (-50, -15))
    ax.add_patch(Circle(fixed[-1], RADIUS, facecolor=INK, edgecolor=INK, zorder=9))
    label(ax, 0.66, 0.27, "regrown", BLUE, 7)
    figure.save(name)


def identity_change(name):
    figure = Figure([("box", "old"), ("box", "new")], connectors=["→"])
    for index, ax in enumerate(figure.axes):
        box(ax)
        bystanders(ax)
        molecule(ax, 0.30, 0.66, 40, n=2, color=ORANGE)
        molecule(ax, 0.70, 0.30, -30, color=BLUE)
        ax.add_patch(Circle((0.47, 0.47), 0.17, facecolor="none", edgecolor=SOFT, lw=0.8, ls=(0, (3, 2)), zorder=1))
    molecule(figure.axes[0], 0.47, 0.47, 20, color=BLUE)
    label(figure.axes[0], 0.47, 0.25, "component A", BLUE, 7)
    molecule(figure.axes[1], 0.47, 0.47, 20, color=BLUE, style="old")
    grow(figure.axes[1], (0.44, 0.44), (50,), color=ORANGE)
    label(figure.axes[1], 0.47, 0.25, "component B", ORANGE, 7)
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# insertion and deletion: one molecule, a pair, a group
# ----------------------------------------------------------------------------------------------------------------


def draw_members(ax, kind, centre, style, cbmc=False, multi=False):
    """The molecule, the pair or the group that a swap inserts, around 'centre'."""
    x, y = centre
    if kind == "single":
        if cbmc and not multi:
            grow(ax, (x - 0.05, y - 0.05), (55, 95), style=style)
        else:
            molecule(ax, x, y, 30, style=style)
        return
    satellites = [(0.13, 0.10)] if kind == "pair" else [(0.14, 0.09), (-0.10, -0.14)]
    sphere(ax, (x, y), 0.25)
    if cbmc:
        trial_positions(ax, [(x - 0.09, y + 0.07), (x + 0.02, y - 0.10)], BLUE)
        for dx, dy in satellites:
            trial_positions(ax, [(x + dx + 0.07, y + dy - 0.09), (x + dx - 0.09, y + dy + 0.08)], ORANGE)
    ion(ax, x, y, "+" if kind == "pair" else "2+", BLUE, style)
    for dx, dy in satellites:
        ion(ax, x + dx, y + dy, "−", ORANGE, style)


def swap(name, kind="single", cbmc=False, multi=False, cfcmc=False):
    bars = {1: True} if cfcmc else None
    figure = Figure([("small", None), ("box", None)], connectors=["⇄"], bars=bars)
    left, ax = figure.axes
    reservoir(left, r"reservoir ($\mu$, $T$)")
    if kind == "single":
        molecule(left, 0.5, 0.5, 60)
    elif kind == "pair":
        ion(left, 0.44, 0.47, "+", BLUE)
        ion(left, 0.57, 0.55, "−", ORANGE)
    else:
        ion(left, 0.50, 0.50, "2+", BLUE)
        ion(left, 0.61, 0.60, "−", ORANGE)
        ion(left, 0.39, 0.40, "−", ORANGE)
    box(ax)
    bystanders(ax, "ring" if kind == "single" else "sparse")
    style = "fractional" if cfcmc else "solid"
    if multi:
        for x, y, angle in ((0.30, 0.66, -20), (0.68, 0.70, 70), (0.66, 0.30, 15)):
            molecule(ax, x, y, angle, style="trial")
        molecule(ax, 0.36, 0.36, 40)
        label(ax, 0.50, 0.52, r"$M$ trial molecules", BLUE, 7)
    else:
        draw_members(ax, kind, (0.50, 0.50) if kind != "single" else (0.50, 0.48), style, cbmc)
    if cfcmc:
        lambda_bar(figure.bars[1], 0.62, left_text="deletion", right_text="insertion")
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# system moves
# ----------------------------------------------------------------------------------------------------------------


def volume(name, anisotropic):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    if anisotropic:
        old = (0.14, 0.20, 0.72, 0.60)
        new = (0.02, 0.26, 0.96, 0.48)
    else:
        old = (0.16, 0.16, 0.68, 0.68)
        new = (0.02, 0.02, 0.96, 0.96)
    box(ax, *new)
    box(ax, *old, style=(0, (3, 2)), color=SOFT, lw=0.9, fill=None)
    for u, v, angle in ((0.22, 0.25, 30), (0.75, 0.24, 100), (0.26, 0.76, -30), (0.74, 0.74, 60), (0.50, 0.50, 10)):
        molecule(ax, old[0] + u * old[2], old[1] + v * old[3], angle, style="old")
        molecule(ax, new[0] + u * new[2], new[1] + v * new[3], angle)
    if anisotropic:
        arrow(ax, (0.86, 0.50), (0.975, 0.50), INK)
        arrow(ax, (0.14, 0.50), (0.025, 0.50), INK)
        arrow(ax, (0.50, 0.80), (0.50, 0.745), INK)
        arrow(ax, (0.50, 0.20), (0.50, 0.255), INK)
        label(ax, 0.50, 0.90, "each box length on its own", INK, 7)
    else:
        for (x0, y0), (x1, y1) in (((0.84, 0.84), (0.97, 0.97)), ((0.16, 0.16), (0.03, 0.03)),
                                   ((0.84, 0.16), (0.97, 0.03)), ((0.16, 0.84), (0.03, 0.97))):
            arrow(ax, (x0, y0), (x1, y1), INK)
        label(ax, 0.50, 0.915, r"$V \rightarrow V + \Delta V$", INK, 7.5)
    figure.save(name)


def hybrid_mc(name):
    figure = Figure([("box", None)])
    ax = figure.axes[0]
    box(ax)
    places = [(0.18, 0.18, 25, 50, 0.4), (0.80, 0.20, 100, 150, -0.4), (0.20, 0.80, -25, -40, 0.4),
              (0.80, 0.82, 55, 215, -0.4), (0.46, 0.44, 10, 80, 0.5)]
    for x, y, angle, direction, rad in places:
        step = 0.15 * np.array([math.cos(math.radians(direction)), math.sin(math.radians(direction))])
        molecule(ax, x, y, angle, style="old")
        molecule(ax, x + step[0], y + step[1], angle + 35)
        arrow(ax, (x, y), (x + step[0], y + step[1]), RED, rad=rad, lw=1.2, scale=8, z=9)
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# Gibbs ensemble moves
# ----------------------------------------------------------------------------------------------------------------


def liquid(ax):
    for x, y, angle in ((0.14, 0.14, 25), (0.86, 0.16, 100), (0.13, 0.86, -25), (0.84, 0.87, 55), (0.50, 0.90, 5),
                        (0.90, 0.52, 85), (0.10, 0.50, 80), (0.50, 0.10, 12), (0.28, 0.72, 40), (0.74, 0.70, -40)):
        molecule(ax, x, y, angle, color=GREY, z=2)


def vapour(ax):
    for x, y, angle in ((0.16, 0.18, 25), (0.84, 0.84, 60), (0.18, 0.84, -30)):
        molecule(ax, x, y, angle, color=GREY, z=2)


def soften(ax):
    for patch in ax.patches:
        if matplotlib.colors.to_hex(patch.get_facecolor()[:3]) == matplotlib.colors.to_hex(GREY):
            patch.set_edgecolor(GREY_EDGE)
    for line in ax.lines:
        if np.allclose(matplotlib.colors.to_rgb(line.get_color()), tuple(0.62 * c for c in matplotlib.colors.to_rgb(GREY))):
            line.set_color(GREY_EDGE)


def gibbs_volume(name):
    figure = Figure([("box", "box A"), ("box", "box B")], connectors=["⇄"],
                    note=r"$V_A + V_B$ constant")
    a, b = figure.axes
    old_a, new_a = (0.14, 0.14, 0.72, 0.72), (0.02, 0.02, 0.96, 0.96)
    old_b, new_b = (0.04, 0.04, 0.92, 0.92), (0.18, 0.18, 0.64, 0.64)
    for ax, old, new, places in ((a, old_a, new_a, ((0.22, 0.24, 30), (0.78, 0.24, 100), (0.22, 0.78, -30),
                                                    (0.78, 0.78, 60), (0.50, 0.50, 10))),
                                 (b, old_b, new_b, ((0.25, 0.28, 30), (0.72, 0.70, 60), (0.30, 0.74, -30)))):
        box(ax, *new)
        box(ax, *old, style=(0, (3, 2)), color=SOFT, lw=0.9, fill=None)
        for u, v, angle in places:
            molecule(ax, old[0] + u * old[2], old[1] + v * old[3], angle, style="old")
            molecule(ax, new[0] + u * new[2], new[1] + v * new[3], angle)
    label(a, 0.50, 0.08, r"$+\Delta V$", INK, 7.5)
    label(b, 0.50, 0.10, r"$-\Delta V$", INK, 7.5)
    figure.save(name)


def gibbs_swap(name, multi=False, cfcmc=None, cbmc=True):
    """cfcmc: None, 'serial' (one fractional molecule, in one of the boxes) or 'parallel' (one in each box)."""
    bars = None
    if cfcmc == "serial":
        bars = {(0, 1): True}
    if cfcmc == "parallel":
        bars = {0: True, 1: True}
    figure = Figure([("box", "box A"), ("box", "box B")], connectors=["→" if not cfcmc else "⇄"], bars=bars)
    a, b = figure.axes
    for ax, fill in ((a, liquid), (b, vapour)):
        box(ax)
        fill(ax)
        soften(ax)
    if cfcmc is None:
        if multi:
            molecule(a, 0.50, 0.50, 20, style="old")
            for x, y, angle in ((0.50, 0.72, -10), (0.50, 0.28, 50)):
                molecule(a, x, y, angle, style="trial")
            for x, y, angle in ((0.30, 0.42, -20), (0.66, 0.70, 70), (0.70, 0.28, 15)):
                molecule(b, x, y, angle, style="trial")
            molecule(b, 0.42, 0.66, 40)
            label(b, 0.52, 0.08, r"$M$ trial molecules", BLUE, 7)
        else:
            grow(a, (0.45, 0.44), (60, 95), style="old")
            grow(b, (0.46, 0.42), (55, 100))
            label(a, 0.50, 0.30, "retrace", SOFT, 7)
            label(b, 0.52, 0.28, "grow", BLUE, 7)
    elif cfcmc == "serial":
        if cbmc:
            grow(a, (0.45, 0.44), (60, 95), style="fractional")
        else:
            molecule(a, 0.50, 0.50, 60, style="fractional")
        molecule(b, 0.50, 0.50, 20, style="trial")
        label(a, 0.50, 0.30, "the fractional molecule", BLUE, 6.8)
        label(b, 0.50, 0.36, "its place when it changes box", SOFT, 6.8)
        lambda_bar(figure.bars[(0, 1)], 0.40, walls=True,
                   track_label="the fractional molecule changes box at the same $\\lambda$")
    else:
        if cbmc:
            grow(a, (0.45, 0.44), (60, 95), style="fractional")
            grow(b, (0.46, 0.42), (55, 100), style="fractional")
        else:
            molecule(a, 0.50, 0.50, 60, style="fractional")
            molecule(b, 0.50, 0.50, 20, style="fractional")
        lambda_bar(figure.bars[0], 0.62, text=r"$\lambda_A + \nu$", directions=(1,), step=0.2)
        lambda_bar(figure.bars[1], 0.38, text=r"$\lambda_B - \nu$", directions=(-1,), step=0.2)
    figure.save(name)


def gibbs_identity_change(name):
    figure = Figure([("box", "box I"), ("box", "box II")], connectors=["⇄"])
    one, two = figure.axes
    for ax in (one, two):
        box(ax)
        bystanders(ax)
        ax.add_patch(Circle((0.47, 0.47), 0.17, facecolor="none", edgecolor=SOFT, lw=0.8, ls=(0, (3, 2)), zorder=1))
    molecule(one, 0.47, 0.47, 20, color=BLUE, style="old")
    grow(one, (0.44, 0.44), (50,), color=ORANGE)
    label(one, 0.47, 0.25, r"A $\rightarrow$ B", INK, 7.5)
    molecule(two, 0.47, 0.47, 50, n=2, color=ORANGE, style="old")
    grow(two, (0.41, 0.45), (15, 30), color=BLUE)
    label(two, 0.47, 0.25, r"B $\rightarrow$ A", INK, 7.5)
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# chemical-potential measurements, parallel tempering
# ----------------------------------------------------------------------------------------------------------------


def widom(name, cfcmc=False, cbmc=True):
    figure = Figure([("box", None)], bars={0: True} if cfcmc else None)
    ax = figure.axes[0]
    box(ax)
    bystanders(ax)
    if cfcmc:
        if cbmc:
            grow(ax, (0.45, 0.44), (55, 95), style="fractional")
        else:
            molecule(ax, 0.50, 0.50, 40, style="fractional")
        lambda_bar(figure.bars[0], 0.45, walls=True, left_text="no deletion", right_text="no insertion")
    else:
        grow(ax, (0.45, 0.44), (55, 95), style="old")
        label(ax, 0.50, 0.30, "ghost molecule", BLUE, 7)
        label(ax, 0.50, 0.225, r"weight $W$ recorded, removed", SOFT, 6.6)
    figure.save(name)


def parallel_tempering(name):
    figure = Figure([("box", r"replica $i$:  $T_i$"), ("box", r"replica $i+1$:  $T_{i+1}$")], connectors=["⇄"],
                    note="the two replicas exchange their configurations")
    cold = [(0.30, 0.30, 20), (0.44, 0.34, 20), (0.58, 0.38, 20), (0.32, 0.46, 20), (0.46, 0.50, 20),
            (0.60, 0.54, 20), (0.34, 0.62, 20), (0.48, 0.66, 20), (0.62, 0.70, 20)]
    hot = [(0.18, 0.20, 70), (0.80, 0.24, 140), (0.22, 0.80, -20), (0.76, 0.78, 50), (0.50, 0.52, 100),
           (0.52, 0.16, 5), (0.14, 0.50, 85), (0.86, 0.52, 30), (0.50, 0.86, -40)]
    for ax, places, color in ((figure.axes[0], cold, BLUE), (figure.axes[1], hot, RED)):
        box(ax)
        for x, y, angle in places:
            molecule(ax, x, y, angle, color=color)
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# reaction ensemble moves: A + B <-> C
# ----------------------------------------------------------------------------------------------------------------


def reactants(ax, style, cbmc=False):
    molecule(ax, 0.32, 0.56, 40, n=2, color=BLUE, style=style)
    molecule(ax, 0.40, 0.34, -20, n=1, color=ORANGE, style=style)


def products(ax, style, cbmc=False):
    if cbmc:
        grow(ax, (0.62, 0.40), (60, 95), color=GREEN, style=style)
    else:
        molecule(ax, 0.66, 0.50, 70, color=GREEN, style=style)


def reaction(name, mode, cbmc):
    """mode: 'cbmc' (whole molecules), 'parallel' (both sides fractional) or 'serial' (one side fractional)."""
    note = r"A + B $\rightleftharpoons$ C"
    if mode == "cbmc":
        figure = Figure([("box", "reactants"), ("box", "products")], connectors=["⇄"], note=note)
        for index, ax in enumerate(figure.axes):
            box(ax)
            bystanders(ax)
        reactants(figure.axes[0], "solid")
        reactants(figure.axes[1], "old")
        products(figure.axes[1], "solid", cbmc=True)
        label(figure.axes[0], 0.36, 0.72, "A", BLUE, 8)
        label(figure.axes[0], 0.50, 0.33, "B", ORANGE, 8)
        label(figure.axes[1], 0.80, 0.52, "C", GREEN, 8)
    elif mode == "parallel":
        figure = Figure([("box", None)], bars={0: True}, note=note)
        ax = figure.axes[0]
        box(ax)
        bystanders(ax)
        reactants(ax, "fractional")
        products(ax, "fractional", cbmc)
        label(ax, 0.34, 0.74, r"reactants: $1-\lambda$", SOFT, 6.8)
        label(ax, 0.68, 0.26, r"products: $\lambda$", SOFT, 6.8)
        lambda_bar(figure.bars[0], 0.55, color=GREEN, left_text="reactants", right_text="products")
    else:
        figure = Figure([("box", "reactant side fractional"), ("box", "product side fractional")],
                        connectors=["⇄"], bars={(0, 1): True}, note=note)
        for ax in figure.axes:
            box(ax)
            bystanders(ax)
        reactants(figure.axes[0], "fractional")
        products(figure.axes[1], "fractional", cbmc)
        lambda_bar(figure.bars[(0, 1)], 0.45, color=GREEN, walls=True,
                   track_label="one side is fractional at a time")
    figure.save(name)


# ----------------------------------------------------------------------------------------------------------------
# the moves, by the anchor of their section
# ----------------------------------------------------------------------------------------------------------------

FIGURES = {
    "translation": lambda n: translation(n, random=False),
    "random-translation": lambda n: translation(n, random=True),
    "rotation": lambda n: rotation(n, random=False),
    "random-rotation": lambda n: rotation(n, random=True),
    "translation-smart-mc": lambda n: smart(n, True, False),
    "rotation-smart-mc": lambda n: smart(n, False, True),
    "translation-rotation-smart-mc": lambda n: smart(n, True, True),
    "reinsertion-cbmc": lambda n: reinsertion(n, multi=False),
    "reinsertion-multi-trial-cbmc": lambda n: reinsertion(n, multi=True),
    "partial-reinsertion-cbmc": partial_reinsertion,
    "identity-change-cbmc": identity_change,
    "swap-conventional": lambda n: swap(n),
    "swap-cbmc": lambda n: swap(n, cbmc=True),
    "swap-multi-trial-cbmc": lambda n: swap(n, cbmc=True, multi=True),
    "swap-cfcmc": lambda n: swap(n, cfcmc=True),
    "swap-cbcfcmc": lambda n: swap(n, cbmc=True, cfcmc=True),
    "pair-swap-conventional": lambda n: swap(n, "pair"),
    "pair-swap-cbmc": lambda n: swap(n, "pair", cbmc=True),
    "pair-swap-cfcmc": lambda n: swap(n, "pair", cfcmc=True),
    "pair-swap-cbcfcmc": lambda n: swap(n, "pair", cbmc=True, cfcmc=True),
    "group-swap-conventional": lambda n: swap(n, "group"),
    "group-swap-cbmc": lambda n: swap(n, "group", cbmc=True),
    "group-swap-cfcmc": lambda n: swap(n, "group", cfcmc=True),
    "group-swap-cbcfcmc": lambda n: swap(n, "group", cbmc=True, cfcmc=True),
    "volume-change": lambda n: volume(n, anisotropic=False),
    "anisotropic-volume-change": lambda n: volume(n, anisotropic=True),
    "hybrid-mc": hybrid_mc,
    "translation-smart-mc-all": lambda n: smart(n, True, False, everything=True),
    "rotation-smart-mc-all": lambda n: smart(n, False, True, everything=True),
    "gibbs-volume-change": gibbs_volume,
    "gibbs-swap-cbmc": lambda n: gibbs_swap(n),
    "gibbs-swap-multi-trial-cbmc": lambda n: gibbs_swap(n, multi=True),
    "gibbs-swap-cfcmc": lambda n: gibbs_swap(n, cfcmc="serial", cbmc=False),
    "gibbs-swap-cbcfcmc": lambda n: gibbs_swap(n, cfcmc="serial", cbmc=True),
    "gibbs-conventional-cfcmc": lambda n: gibbs_swap(n, cfcmc="parallel", cbmc=False),
    "gibbs-conventional-cbcfcmc": lambda n: gibbs_swap(n, cfcmc="parallel", cbmc=True),
    "gibbs-identity-change-cbmc": gibbs_identity_change,
    "widom-insertion": lambda n: widom(n),
    "widom-cfcmc": lambda n: widom(n, cfcmc=True, cbmc=False),
    "widom-cbcfcmc": lambda n: widom(n, cfcmc=True, cbmc=True),
    "parallel-tempering": parallel_tempering,
    "reaction-cbmc": lambda n: reaction(n, "cbmc", True),
    "reaction-conventional-cfcmc": lambda n: reaction(n, "parallel", False),
    "reaction-conventional-cbcfcmc": lambda n: reaction(n, "parallel", True),
    "reaction-cfcmc": lambda n: reaction(n, "serial", False),
    "reaction-cbcfcmc": lambda n: reaction(n, "serial", True),
}


CAPTIONS = {
    "translation": "A molecule is displaced over a small random distance along one axis; its orientation is kept.",
    "random-translation": "A molecule is put anywhere along one axis of the box; its orientation is kept.",
    "rotation": "A molecule is rotated over a small random angle around an axis through its rotation centre.",
    "random-rotation": "A molecule is given a new orientation drawn uniformly from all orientations.",
    "translation-smart-mc": "The displacement is drawn around a drift along the force F on the molecule (shaded: where the proposals fall).",
    "rotation-smart-mc": "The rotation is drawn around a drift along the torque on the molecule.",
    "translation-rotation-smart-mc": "Displacement along the force and rotation along the torque in one trial.",
    "reinsertion-cbmc": "The molecule is retraced where it is and grown bead by bead at a random position; at every bead one of several trial positions (dashed) is chosen.",
    "reinsertion-multi-trial-cbmc": "M whole molecules are grown at random positions (dashed) and one of them is chosen as the new configuration.",
    "partial-reinsertion-cbmc": "A part of the molecule stays fixed (dark) and the rest is regrown with trial positions.",
    "identity-change-cbmc": "A molecule of component A is removed and a molecule of component B is grown in its place.",
    "swap-conventional": "A molecule is inserted from, or deleted to, a reservoir at a random position and orientation.",
    "swap-cbmc": "The inserted molecule is grown bead by bead from several trial positions (dashed); a deleted molecule is retraced the same way.",
    "swap-multi-trial-cbmc": "M whole trial molecules are grown (dashed) and one of them is inserted.",
    "swap-cfcmc": "A fractional molecule (hatched) is coupled to the system with lambda; a trial changes lambda, and an insertion or deletion happens when lambda leaves the interval.",
    "swap-cbcfcmc": "As the CFCMC swap, with the fractional molecule grown and retraced with trial positions (dashed).",
    "pair-swap-conventional": "A pair of molecules is inserted or deleted together; the second is placed within R_max of the first.",
    "pair-swap-cbmc": "The pair is inserted or deleted together, each molecule from several trial positions (dashed).",
    "pair-swap-cfcmc": "A fractional pair (hatched) shares one lambda, so that both molecules appear and disappear together.",
    "pair-swap-cbcfcmc": "As the CFCMC pair swap, with the fractional pair grown and retraced with trial positions (dashed).",
    "group-swap-conventional": "A group of molecules is inserted or deleted together; the satellites are placed within R_max of the central molecule.",
    "group-swap-cbmc": "The group is inserted or deleted together, each member from several trial positions (dashed).",
    "group-swap-cfcmc": "A fractional group (hatched) shares one lambda, so that the whole group appears and disappears together.",
    "group-swap-cbcfcmc": "As the CFCMC group swap, with the fractional group grown and retraced with trial positions (dashed).",
    "volume-change": "The box is scaled in all directions (dashed: the old box) and the molecules move with their centres of mass.",
    "anisotropic-volume-change": "The box lengths change independently (dashed: the old box) at fixed angles.",
    "hybrid-mc": "All molecules follow a short molecular-dynamics trajectory, which is accepted or rejected as one move.",
    "translation-smart-mc-all": "All molecules are displaced at once, each around a drift along the force on it.",
    "rotation-smart-mc-all": "All molecules are rotated at once, each around a drift along the torque on it.",
    "gibbs-volume-change": "Volume is exchanged between the two boxes at constant total volume (dashed: the old boxes).",
    "gibbs-swap-cbmc": "A molecule is retraced in one box and grown in the other, each bead from several trial positions (dashed).",
    "gibbs-swap-multi-trial-cbmc": "M trial molecules are grown in the receiving box and one is chosen; the donating box is weighed with the transferred molecule and M-1 reference molecules.",
    "gibbs-swap-cfcmc": "One fractional molecule (hatched) lives in one of the two boxes; a trial changes lambda or moves the fractional molecule to the other box at the same lambda.",
    "gibbs-swap-cbcfcmc": "As the serial CFCMC Gibbs swap, with the fractional molecule grown and retraced with trial positions (dashed).",
    "gibbs-conventional-cfcmc": "Each box holds a fractional molecule (hatched); one random change moves the two coupling parameters in opposite directions.",
    "gibbs-conventional-cbcfcmc": "As the parallel CFCMC Gibbs move, with the fractional molecules grown and retraced with trial positions (dashed).",
    "gibbs-identity-change-cbmc": "A molecule of A becomes B in one box while a molecule of B becomes A in the other.",
    "widom-insertion": "A ghost molecule is grown, its Rosenbluth weight is recorded, and it is removed again; the configuration does not change.",
    "widom-cfcmc": "The fractional molecule (hatched) only changes lambda; a trial beyond either end is rejected.",
    "widom-cbcfcmc": "As the CFCMC Widom move, with the fractional molecule grown with trial positions (dashed).",
    "parallel-tempering": "Two replicas at neighbouring state points exchange their configurations.",
    "reaction-cbmc": "The reactants (A and B) are retraced and removed and the products (C) are grown, or the reverse.",
    "reaction-conventional-cfcmc": "Fractional molecules (hatched) of both sides are present: the products are scaled with lambda and the reactants with 1 - lambda.",
    "reaction-conventional-cbcfcmc": "As the parallel reaction move, with the fractional molecules grown and retraced with trial positions (dashed).",
    "reaction-cfcmc": "Only one side of the reaction is fractional (hatched) at a time; trials change lambda or switch the fractional side.",
    "reaction-cbcfcmc": "As the serial reaction move, with the fractional molecules grown and retraced with trial positions (dashed).",
}

EXPLANATIONS = {
    "translation": r"""
The faded molecule with the dashed outline is the selected molecule at its old position, the fully coloured one the
trial position; the arrow is the displacement \f$\Delta r\f$ along the chosen axis. The bar under the molecule is
the range \f$[-\Delta_\text{max},+\Delta_\text{max}]\f$ from which the displacement is drawn. The orientation of
the molecule is the same before and after, and the grey molecules do not move.""",
    "random-translation": r"""
The dotted line is the chosen axis through the old position (faded). The new position (full colour) is drawn
uniformly along the whole box length on that line, so the arrow can span any distance; the orientation is kept.""",
    "rotation": r"""
The cross marks the rotation centre. The faded molecule is the old orientation and the fully coloured one the new
orientation; the curved arrow is the rotation angle, drawn from
\f$[-\Delta\theta_\text{max},+\Delta\theta_\text{max}]\f$. The position of the rotation centre does not change.""",
    "random-rotation": r"""
The dotted circle stands for all orientations of the molecule around its rotation centre (the cross). The new
orientation (full colour) is drawn uniformly from them, without reference to the old one (faded), so the arrow can
cover any angle.""",
    "translation-smart-mc": r"""
The red arrow is the force \f$\mathbf{F}\f$ on the molecule at its old position (faded). The shaded disc shows
where the proposals fall: a Gaussian around a point shifted along the force. The new position (full colour) is one
draw from that distribution, so it lies in the direction of the force but not exactly on the arrow.""",
    "rotation-smart-mc": r"""
The curved red arrow is the torque \f$\boldsymbol{\tau}\f$ on the molecule. The new orientation (full colour) is
drawn around a rotation along the torque, starting from the old orientation (faded); the centre of the molecule
stays where it is.""",
    "translation-rotation-smart-mc": r"""
Both biases act in one trial: the molecule is displaced into the shaded region around a drift along the force
\f$\mathbf{F}\f$ (straight red arrow) and rotated around a drift along the torque \f$\boldsymbol{\tau}\f$ (curved
red arrow). The faded molecule is the old configuration, the fully coloured one the new configuration.""",
    "reinsertion-cbmc": r"""
On the left the old configuration (faded) is retraced: next to every bead the thin dashed circles are the other
trial positions that are generated to weigh it, which gives the Rosenbluth weight \f$W(\mathbf{o})\f$. On the right
the molecule is grown at a random position: at every bead several trial positions are generated and one is chosen
(full colour), which gives \f$W(\mathbf{n})\f$. The arrow connects the old and the new configuration.""",
    "reinsertion-multi-trial-cbmc": r"""
Each thin dashed molecule is a complete trial molecule grown with CBMC at a random position; there are \f$M\f$ of
them in total. One is chosen with a probability proportional to its weight and becomes the new configuration (full
colour), which replaces the old one (faded). The acceptance uses the mean weight of the whole set.""",
    "partial-reinsertion-cbmc": r"""
The dark beads are the atoms that stay fixed. The faded beads are the old positions of the remaining atoms; the
fully coloured beads are regrown from the last fixed atom, each chosen from several trial positions (thin dashed
circles). Only the regrown part contributes to the Rosenbluth weights.""",
    "identity-change-cbmc": r"""
The dashed circle marks the selected molecule. In the left panel it is a molecule of component A (blue). In the
right panel that molecule is removed (faded) and a molecule of component B (orange) is grown from its position,
with trial positions shown as thin dashed circles. All other molecules, of both components, are unchanged.""",
    "swap-conventional": r"""
The dashed panel on the left is the reservoir at the imposed chemical potential and temperature; the double arrow
is the exchange of a molecule with the box. For an insertion the molecule is placed at a uniformly random position
with a random orientation, as drawn; for a deletion a randomly chosen molecule of the component is removed.""",
    "swap-cbmc": r"""
As in the conventional swap a molecule is exchanged with the reservoir (dashed panel), but an inserted molecule is
grown bead by bead: the thin dashed circles are the trial positions that were generated and not chosen at each
bead. For a deletion the selected molecule is retraced with the same number of trial positions.""",
    "swap-multi-trial-cbmc": r"""
For an insertion \f$M\f$ complete trial molecules are grown at random positions; the thin dashed molecules are the
ones not chosen and the fully coloured molecule is the one that is inserted. For a deletion the picture is read in
reverse: the selected molecule is retraced and \f$M-1\f$ reference molecules are grown to complete the set.""",
    "swap-cfcmc": r"""
The hatched molecule in the box is the fractional molecule, and the axis under the box is its coupling parameter
\f$\lambda\f$ with the present value as a dot. A trial changes \f$\lambda\f$ by a random amount (the two arrows).
As long as the new value stays between 0 and 1 only the coupling of the fractional molecule changes. Past the
right end the fractional molecule becomes a whole molecule and a new fractional molecule is inserted (an
insertion); past the left end the fractional molecule is removed and a whole molecule becomes fractional (a
deletion). The reservoir panel on the left indicates that the number of whole molecules changes this way.""",
    "swap-cbcfcmc": r"""
The picture is that of the CFCMC swap: a hatched fractional molecule and its \f$\lambda\f$ axis. The difference is
at the two ends of the axis, where a fractional molecule is created or removed: this is done with CBMC, and the
thin dashed circles are the trial positions of that growth or retrace.""",
    "pair-swap-conventional": r"""
The two ions stand for the selected component (blue) and its partner component (orange), which are inserted or
deleted together. The first molecule is placed at a random position; the dashed circle is the sphere of radius
\f$R_\text{max}\f$ around it, and the second molecule is placed uniformly inside this sphere. A deletion removes
a pair whose distance is within \f$R_\text{max}\f$.""",
    "pair-swap-cbmc": r"""
The pair is placed as in the conventional pair swap, inside the sphere of radius \f$R_\text{max}\f$ (dashed
circle), but both molecules are grown with CBMC: the thin dashed circles in the colour of each molecule are its
trial positions that were not chosen.""",
    "pair-swap-cfcmc": r"""
Both molecules of the pair are fractional (hatched) and share the single coupling parameter \f$\lambda\f$ of the
axis under the box, so they are coupled and decoupled together. Past the right end of the axis the fractional pair
becomes whole and a new fractional pair is inserted, the second molecule within \f$R_\text{max}\f$ of the first
(dashed circle); past the left end the fractional pair is removed and a whole pair becomes fractional.""",
    "pair-swap-cbcfcmc": r"""
The picture is that of the CFCMC pair swap, with a hatched fractional pair that shares one \f$\lambda\f$. Where a
fractional pair is created or removed, at the ends of the axis, the molecules are grown or retraced with CBMC; the
thin dashed circles are the trial positions.""",
    "group-swap-conventional": r"""
The central molecule (blue) and its satellites (orange) form the group that is inserted or deleted as a whole; a
divalent ion with two counter-ions is drawn as an example. The central molecule is placed at a random position and
every satellite uniformly inside the sphere of radius \f$R_\text{max}\f$ around it (dashed circle).""",
    "group-swap-cbmc": r"""
The group is placed as in the conventional group swap, with the satellites inside the sphere of radius
\f$R_\text{max}\f$ (dashed circle), and every member is grown with CBMC: the thin dashed circles are the trial
positions of each member that were not chosen.""",
    "group-swap-cfcmc": r"""
All members of the group are fractional (hatched) and share the single coupling parameter \f$\lambda\f$ of the
axis under the box. Past the right end of the axis the fractional group becomes whole and a new fractional group
is inserted, the satellites within \f$R_\text{max}\f$ of the central molecule (dashed circle); past the left end
the fractional group is removed and a whole group becomes fractional.""",
    "group-swap-cbcfcmc": r"""
The picture is that of the CFCMC group swap, with a hatched fractional group that shares one \f$\lambda\f$. Where
a fractional group is created or removed, at the ends of the axis, its members are grown or retraced with CBMC;
the thin dashed circles are the trial positions.""",
    "volume-change": r"""
The dashed square is the old box and the solid square the new one; the arrows show that it changes by the same
factor in all directions. Every molecule moves with its centre of mass from the faded to the fully coloured
position and keeps its orientation and internal geometry, so only the distances between molecules change.""",
    "anisotropic-volume-change": r"""
The dashed rectangle is the old box and the solid one the new box: here one box length grows while another
shrinks, each by its own factor, at fixed box angles. The molecules keep the fractional position of their centre
of mass (faded to full colour) and are not deformed.""",
    "hybrid-mc": r"""
Every molecule moves at the same time: the red arrows are the short molecular-dynamics trajectories that take the
old configuration (faded) to the new one (full colour), with translation and rotation together. The whole
trajectory is accepted or rejected as a single move on the change of the total energy.""",
    "translation-smart-mc-all": r"""
The red arrows are the forces on the molecules in the old configuration (faded). All molecules are displaced in
the same trial, each around a drift along its own force, to the fully coloured positions; one acceptance test is
applied to the whole set.""",
    "rotation-smart-mc-all": r"""
The curved red arrows are the torques on the molecules. All molecules are rotated in the same trial about their
own centres, each around a drift along its torque, from the faded to the fully coloured orientation; one
acceptance test is applied to the whole set.""",
    "gibbs-volume-change": r"""
The two panels are the two boxes of the Gibbs ensemble. Box A grows from its old size (dashed) to the new one
(solid) by \f$\Delta V\f$ while box B shrinks by the same amount, so the total volume is constant. In both boxes
the molecules move with their centres of mass, from the faded to the fully coloured positions.""",
    "gibbs-swap-cbmc": r"""
The left panel is the donating box (drawn dense) and the right panel the receiving box (drawn dilute). The
transferred molecule is retraced in the donating box (faded, with its trial positions) and grown bead by bead at a
random position in the receiving box (full colour); the thin dashed circles are the trial positions that were not
chosen. The direction of the transfer is chosen at random.""",
    "gibbs-swap-multi-trial-cbmc": r"""
In the receiving box (right) \f$M\f$ complete trial molecules are grown; the thin dashed ones are not chosen and
the fully coloured one is the transferred molecule. In the donating box (left) the transferred molecule (faded) is
retraced and \f$M-1\f$ reference molecules (thin dashed) are grown to complete its set.""",
    "gibbs-swap-cfcmc": r"""
There is one fractional molecule (hatched), here in box A; the dashed outline in box B is the place it takes when
it changes box. The axis under the two boxes is its coupling parameter \f$\lambda\f$: one sub-move changes
\f$\lambda\f$ inside the box, between the two ends of the axis. The other sub-moves take the fractional molecule to
the other box at the same \f$\lambda\f$, either by placing it there at a random position or by turning it whole
here while a whole molecule of the other box becomes fractional, which is what transfers a molecule.""",
    "gibbs-swap-cbcfcmc": r"""
The picture is that of the serial CFCMC Gibbs swap: one hatched fractional molecule, one \f$\lambda\f$ axis for
both boxes, and the dashed outline where it appears in the other box. Where a fractional configuration is created
or removed CBMC is used; the thin dashed circles are its trial positions.""",
    "gibbs-conventional-cfcmc": r"""
Each box has its own fractional molecule (hatched) and its own axis. One random number \f$\nu\f$ changes both
coupling parameters in opposite directions, as the two arrows show: \f$\lambda_A\f$ goes up by \f$\nu\f$ while
\f$\lambda_B\f$ goes down by \f$\nu\f$. When a coupling parameter passes the right end of its axis the fractional
molecule of that box becomes whole and a new one is inserted; when it passes the left end the fractional molecule
is removed and a whole molecule becomes fractional. Together this moves a molecule from one box to the other.""",
    "gibbs-conventional-cbcfcmc": r"""
The picture is that of the parallel CFCMC Gibbs move, with a hatched fractional molecule and an axis in each box
and one change \f$\nu\f$ that moves the two coupling parameters in opposite directions. The fractional molecules
that are created and removed when an end is passed are grown and retraced with CBMC; the thin dashed circles are
the trial positions.""",
    "gibbs-identity-change-cbmc": r"""
The dashed circles mark the selected molecule of each box. In box I a molecule of component A (blue, faded) is
replaced by a molecule of component B (orange) grown at its position; in box II a molecule of B (faded) is
replaced by a molecule of A. The two changes are made together, so the number of molecules of each component in
the two boxes together is the same before and after.""",
    "widom-insertion": r"""
The faded molecule is the ghost molecule: it is grown at a random position, with the trial positions drawn as thin
dashed circles, only to obtain its Rosenbluth weight \f$W\f$. It is then removed again, so the grey molecules, and
with them the configuration of the system, are never changed by this move.""",
    "widom-cfcmc": r"""
The hatched molecule is the fractional molecule and the axis under the box its coupling parameter \f$\lambda\f$.
The ends of the axis are drawn as walls: a trial that would take \f$\lambda\f$ past either end is rejected, so no
molecule is inserted or deleted and the number of whole molecules stays fixed. The chemical potential follows from
how often the two ends of the axis are visited.""",
    "widom-cbcfcmc": r"""
The picture is that of the CFCMC Widom move: a hatched fractional molecule whose \f$\lambda\f$ stays between the
two walls of the axis. The thin dashed circles are the trial positions of the configurational-bias growth of the
fractional molecule.""",
    "parallel-tempering": r"""
The two panels are two replicas of the same system at neighbouring state points, here a colder, ordered one
(blue) and a hotter, disordered one (red). The double arrow is the exchange: after an accepted move the
configuration of each replica continues at the state point of the other.""",
    "reaction-cbmc": r"""
The example is the reaction A + B ⇌ C. The left panel shows the reactant molecules A (blue) and B (orange) that
are selected. In the right panel they are removed (faded) and the product molecule C (green) is grown with CBMC;
the thin dashed circles are its trial positions. The backward reaction reads the picture from right to left.""",
    "reaction-conventional-cfcmc": r"""
The box holds fractional molecules (hatched) of both sides of the reaction A + B ⇌ C at the same time: the
reactants A and B, coupled with \f$1-\lambda\f$, and the product C, coupled with \f$\lambda\f$. Moving the dot
along the axis under the box couples one side while it decouples the other. Past the right end the products become
whole molecules and whole reactant molecules become the new fractional ones (a forward reaction); past the left
end the reverse happens.""",
    "reaction-conventional-cbcfcmc": r"""
The picture is that of the parallel reaction move, with fractional molecules of both sides and one \f$\lambda\f$
axis. The nearly decoupled fractional molecules that are created and removed when an end of the axis is passed are
grown and retraced with CBMC; the thin dashed circles are the trial positions.""",
    "reaction-cfcmc": r"""
Only one side of the reaction A + B ⇌ C is fractional at a time: the left panel is the state with fractional
reactants (hatched A and B), the right panel the state with a fractional product (hatched C). The axis under the
panels is the coupling parameter \f$\lambda\f$ of the fractional side; a trial that takes it past an end is
rejected. The double arrow stands for the two sub-moves that switch sides at the same \f$\lambda\f$: the
fractional molecules of one side are replaced by fractional molecules of the other side, or they become whole
while whole molecules of the other side become fractional.""",
    "reaction-cbcfcmc": r"""
The picture is that of the serial reaction move, with one side of the reaction fractional at a time and one
\f$\lambda\f$ axis. The fractional molecules that are inserted when the fractional side switches are grown with
CBMC, and the ones that are removed are retraced; the thin dashed circles are the trial positions.""",
}

MANUAL = pathlib.Path(__file__).resolve().parent / "manual" / "mc_moves.md"


BEGIN, END = "<!-- figure of the move: {} -->", "<!-- end of the figure of the move -->"


def figure_block(name):
    """The lines of the manual for a move: the figure with its caption, sized from the drawing (130 pixels and
    2.2 cm per inch), and the paragraph that explains it."""
    width_inch = plt.imread(OUTPUT / f"mc_move_{name}.png").shape[1] / 220.0
    caption = CAPTIONS[name]
    return [BEGIN.format(name),
            f'\\image html mc_move_{name}.png "{caption}" width={round(130 * width_inch)}px',
            f'\\image latex mc_move_{name}.png "{caption}" width={2.2 * width_inch:.1f}cm',
            "",
            *EXPLANATIONS[name].strip().split("\n"),
            END]


def update_manual():
    """Puts the figure of every move, with its explanation, under the 'Keyword:' paragraph of its section, in
    place of the one that is there."""
    lines, skipping = [], False
    for line in MANUAL.read_text(encoding="utf-8").split("\n"):
        if line.startswith("<!-- figure of the move:"):
            skipping = True
            # the blank line before the block went in with it
            if lines and lines[-1] == "":
                lines.pop()
        elif skipping:
            skipping = line != END
        elif not line.startswith(("\\image html mc_move_", "\\image latex mc_move_")):
            lines.append(line)
    out, pending, seen = [], None, set()
    index = 0
    while index < len(lines):
        line = lines[index]
        out.append(line)
        if line.startswith("#") and '<a name="' in line:
            anchor = line.split('<a name="')[1].split('"')[0]
            pending = anchor if anchor in FIGURES else None
        elif pending and line.startswith("Keyword:"):
            # to the end of the paragraph
            while index + 1 < len(lines) and lines[index + 1].strip():
                index += 1
                out.append(lines[index])
            out.append("")
            out.extend(figure_block(pending))
            seen.add(pending)
            pending = None
        index += 1
    missing = set(FIGURES) - seen
    if missing:
        raise SystemExit(f"no section with a 'Keyword:' paragraph for: {sorted(missing)}")
    # one blank line after a block
    tidy = []
    for line in out:
        if line == "" and len(tidy) >= 2 and tidy[-1] == "" and tidy[-2] == END:
            continue
        tidy.append(line)
    MANUAL.write_text("\n".join(tidy), encoding="utf-8")


if __name__ == "__main__":
    names = sys.argv[1:] or list(FIGURES)
    for name in names:
        FIGURES[name](name)
    print(f"{len(names)} figures in {OUTPUT}")
    update_manual()
    print(f"{MANUAL} shows them")
