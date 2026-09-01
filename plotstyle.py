"""House plotting style (colourblind-safe, Charter serif)."""
import matplotlib as mpl

BLUE, AMBER, PURPLE, TEAL = "#004C8C", "#B26B00", "#8E2F6E", "#006D77"
INK, GREY = "#222222", "#9A9A9A"


def apply_style():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Charter", "PT Serif", "Palatino", "Georgia", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 13, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "legend.fontsize": 9.5, "legend.frameon": False,
        "axes.linewidth": 1.0, "figure.dpi": 120, "savefig.dpi": 300,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": False, "ytick.right": True,
    })


def panel_label(ax, text, loc="upper left"):
    xy = {"upper left": (0, 1), "upper right": (1, 1),
          "lower left": (0, 0), "lower right": (1, 0)}[loc]
    off = {"upper left": (6, -6), "upper right": (-6, -6),
           "lower left": (6, 6), "lower right": (-6, 6)}[loc]
    ha = "left" if "left" in loc else "right"
    va = "top" if "upper" in loc else "bottom"
    ax.annotate(text, xy=xy, xycoords="axes fraction", textcoords="offset points",
                xytext=off, ha=ha, va=va, style="italic", fontsize=12,
                bbox=dict(fc="white", ec="none", alpha=0.72, pad=1.5))
