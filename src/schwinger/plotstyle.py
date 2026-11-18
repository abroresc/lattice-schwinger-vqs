"""
Shared figure style. One place, so every figure in the report matches.

The categorical palette is colourblind-safe: it was checked with a Delta-E
separation validator under deuteranopia and tritanopia rather than eyeballed.
Worst adjacent pair separation is 11.0 (deutan), well above the floor of 8.
Series are also distinguished by marker shape, so identity never depends on
colour alone. That matters for a printed thesis and for a projected talk.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fixed order. Never cycled, never reassigned by rank.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#7B52AB", "#B8860B"]
MARKERS = ["o", "s", "^", "D", "v"]

INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#dcdcdc"


def use_report_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 4.5,
    })


def annotate(ax, text: str) -> None:
    """A short note under a panel, in muted ink. Used to state what a check proves."""
    ax.text(0.0, -0.26, text, transform=ax.transAxes, fontsize=7.5, color=MUTED,
            va="top", ha="left", wrap=True)


def save(fig, path) -> None:
    """
    Save a figure as PNG and PDF side by side.

    PNG for the README and for pasting into slides. PDF because it is vector,
    so it stays sharp at any zoom when included in the LaTeX thesis.
    """
    import pathlib

    p = pathlib.Path(path)
    fig.savefig(p.with_suffix(".png"))
    # No creation date in the PDF, so rerunning a script on unchanged data gives
    # a byte-identical file and git shows a figure as changed only when it is.
    fig.savefig(p.with_suffix(".pdf"), metadata={"CreationDate": None})
