"""
10_figure.py - main figure
Every number is read from results/rev5_main_table.csv (Stata output of ldl_rev5_main.do).
A  % above own 2026 goal by Friedewald, Friedewald or paired non-HDL-C, Martin-Hopkins
   (Korn-Graubard 95% CI). Hollow marker = estimate does not meet NCHS presentation standards.
B  Martin-Hopkins minus Friedewald, and minus (Friedewald or non-HDL-C) (MOVER 95% CI);
   difference and unweighted up/down counts printed at right.
180 mm wide; text >= 7 pt; colours distinct for colour-blind readers plus distinct marker shapes.
Outputs: figures/Figure_LDL_equations.tiff (600 dpi), .pdf, .png
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results" / "rev5_main_table.csv")
d = d[d.block.str.startswith("A.")].set_index("label")

ROWS = [  # (csv label, shown label, kind)  kind: 'all' | 'head' | 'row'
    ("All adults with a defined goal", "All with a defined goal", "all"),
    (None, "Risk group", "head"),
    ("Primary prevention, borderline-high", "Primary prevention", "row"),
    ("ASCVD", "ASCVD", "row"),
    (None, "Treatment and TG", "head"),
    ("Treated, TG >=150", "Treated, TG \u2265150", "row"),
    ("Treated, TG <150", "Treated, TG <150", "row"),
    ("Not treated", "Not treated", "row"),
    (None, "Triglycerides, mg/dL", "head"),
    ("<100", "<100", "row"),
    ("100-149", "100\u2013149", "row"),
    ("150-199", "150\u2013199", "row"),
    ("200-399", "200\u2013399", "row"),
]
INK, INK2, GRID, BAND = "#111111", "#555555", "#e3e2de", "#f1f0ec"
C_F, C_E, C_M = "#eb6834", "#1baf7a", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
MINUS = "\u2212"
def sgn(v):
    return f"{v:+.1f}".replace("-", MINUS)
def bad(flag):
    return isinstance(flag, str) and flag.startswith("SUPPRESS")

W = 180 / 25.4
fig = plt.figure(figsize=(W, 6.1))
axA = fig.add_axes([0.315, 0.275, 0.245, 0.665])
axB = fig.add_axes([0.588, 0.275, 0.245, 0.665], sharey=axA)

n = len(ROWS)
ys = {i: n - 1 - i for i in range(n)}
off = 0.23
for i, (key, lab, kind) in enumerate(ROWS):
    y = ys[i]
    if kind == "head":
        for ax in (axA, axB):
            ax.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0, lw=0)
        fig.text(0.168, 0, "")  # placeholder to keep layout stable
        axA.text(-0.64, y, lab, transform=axA.get_yaxis_transform(), ha="left", va="center",
                 fontsize=7.5, fontweight="bold", color=INK)
        continue
    r = d.loc[key]
    axA.text(-0.64 if kind == "all" else -0.61, y, lab, transform=axA.get_yaxis_transform(), ha="left",
             va="center", fontsize=7.5, fontweight="bold" if kind == "all" else "normal", color=INK)
    for dy, p, lo, hi, c, m, flag in [(off, r.pF, r.lF, r.uF, C_F, "s", r.flagF),
                                      (0, r.pE, r.lE, r.uE, C_E, "^", "ok"),
                                      (-off, r.pM, r.lM, r.uM, C_M, "o", r.flagM)]:
        axA.plot([lo, hi], [y + dy] * 2, color=c, lw=1.1, solid_capstyle="butt", zorder=2)
        axA.plot(p, y + dy, marker=m, ms=4.6, color=c, mfc="white" if bad(flag) else c, mec=c, mew=1.0, zorder=3)
    axB.plot([r.gap_lb, r.gap_ub], [y + off / 1.6] * 2, color=INK, lw=1.1, zorder=2)
    axB.plot(r.gap, y + off / 1.6, marker="o", ms=4.6, color=INK, zorder=3)
    axB.plot([r.res_lb, r.res_ub], [y - off / 1.6] * 2, color=INK2, lw=1.0, zorder=2)
    axB.plot(r.res, y - off / 1.6, marker="D", ms=4.2, mfc="white", mec=INK2, mew=1.0, zorder=3)
    axB.text(1.03, y, f"{sgn(r.gap)}  ({int(r.n_up)}/{int(r.n_dn)})", transform=axB.get_yaxis_transform(),
             ha="left", va="center", fontsize=7.5, color=INK)

axB.text(1.03, n - 0.42, "MH \u2212 F (up/down)", transform=axB.get_yaxis_transform(), ha="left", va="center",
         fontsize=7, fontweight="bold", color=INK2)
axA.set_xlim(40, 95)
axA.set_xticks(range(40, 91, 10))
axB.set_xlim(-8, 22)
axB.set_xticks([-5, 0, 5, 10, 15, 20])
axB.axvline(0, color=INK2, lw=0.9, zorder=1)
axA.set_xlabel("Above own 2026 goal, % (95% CI)", fontsize=7.5)
axB.set_xlabel("Percentage points (95% CI)", fontsize=7.5)
for ax in (axA, axB):
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks([])
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=7.5, length=2.5)
axA.set_title("A  Above goal, by measure", loc="left", fontsize=8.5, fontweight="bold", color=INK, pad=6)
axB.set_title("B  Difference from Martin\u2013Hopkins", loc="left", fontsize=8.5, fontweight="bold", color=INK, pad=6)

from matplotlib.lines import Line2D
hA = [Line2D([0], [0], color=C_F, lw=1.1, marker="s", ms=4.6, label="Friedewald LDL-C"),
      Line2D([0], [0], color=C_E, lw=1.1, marker="^", ms=4.6, label="Friedewald LDL-C or non-HDL-C"),
      Line2D([0], [0], color=C_M, lw=1.1, marker="o", ms=4.6, label="Martin\u2013Hopkins LDL-C")]
hB = [Line2D([0], [0], color=INK, lw=1.1, marker="o", ms=4.6, label="Martin\u2013Hopkins \u2212 Friedewald"),
      Line2D([0], [0], color=INK2, lw=1.0, marker="D", ms=4.2, mfc="white", mec=INK2,
             label="Martin\u2013Hopkins \u2212 (Friedewald or non-HDL-C)")]
axA.legend(handles=hA, loc="upper left", bbox_to_anchor=(-0.10, -0.105), frameon=False, fontsize=7,
           handlelength=2.2, labelspacing=0.35)
axB.legend(handles=hB, loc="upper left", bbox_to_anchor=(-0.02, -0.105), frameon=False, fontsize=7,
           handlelength=2.2, labelspacing=0.35)
fig.text(0.03, 0.022,
         "Hollow marker: estimate does not meet NCHS presentation standards.\n"
         "Up/down: unweighted participants above goal by MH only / by F only. TG, triglycerides (mg/dL); MH, Martin\u2013Hopkins; F, Friedewald.\n"
         "Adults aged 30\u201379 y with a defined 2026 goal, NHANES August 2021\u2013August 2023.",
         fontsize=7, color=INK2, va="bottom", ha="left", linespacing=1.35)

(ROOT / "figures").mkdir(exist_ok=True)
out = str(ROOT / "figures" / "Figure_LDL_equations")
fig.savefig(out + ".png", dpi=300, facecolor="white")
fig.savefig(out + ".pdf", facecolor="white")
fig.savefig(out + ".tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("saved", out)
