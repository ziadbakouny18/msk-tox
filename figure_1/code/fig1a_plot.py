#!/usr/bin/env python3
"""Figure 1a - RAG/LLM pipeline for toxicity detection.
Lives in  <figure 1>/scripts/  and writes to  <figure 1>/results/fig1a.pdf

    python fig1a_plot.py
    ALLOW_FONT_FALLBACK=1 python fig1a_plot.py    
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle, FancyArrow

# ---- Hard-fail if Arial isn't actually resolved (no silent fallback) ----
FONT = "Arial"
_path = fm.findfont(FONT, fallback_to_default=False) if not os.environ.get(
    "ALLOW_FONT_FALLBACK") else None
if _path is not None and "rial" not in os.path.basename(_path).lower():
    raise RuntimeError(
        f"Arial not found -- matplotlib resolved to '{_path}' instead. "
        "Install Arial or set ALLOW_FONT_FALLBACK=1 for a proof render."
    )
if os.environ.get("ALLOW_FONT_FALLBACK"):
    FONT = "Liberation Sans"          # metric-compatible stand-in
else:
    print(f"Arial resolved to: {_path}")

matplotlib.rcParams["font.family"] = FONT
matplotlib.rcParams["pdf.fonttype"] = 42      # keep text editable in the PDF
matplotlib.rcParams["svg.fonttype"] = "none"  # keep text as text in the SVG

W, H = 310.0, 194.4
INK, MUTED, BLUE = "#1A1A1A", "#5C6672", "#1F4E79"
HDR, ZEBRA, RULE = "#2E6DA4", "#F5F7FA", "#9FB6D0"

OUT = Path(__file__).resolve().parent.parent / "results" / "main" / "fig1a.pdf"

# ---- layout -------------------------------------------------------------
C1, C1W = 5.0, 68.0          # clinical notes
G1 = (73.0, 97.0)            # 24 pt gutter -> RAG
C2, C2W = 97.0, 66.0         # retrieval
G2 = (163.0, 187.0)          # 24 pt gutter -> LLM
C3, C3W = 187.0, 45.0        # AE calls
G3 = (232.0, 248.0)          # 16 pt gutter
C4, C4W = 248.0, 57.0        # MSK-Tox dataset
TOP = 14.0                   # top edge of every table / card
ARROW_Y = 95.0
BODY, LEAD = 6.0, 10.3

fig = plt.figure(figsize=(W / 72, H / 72))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)          # invert y so coordinates read top-down
ax.axis("off")


def box(x, y, w, h, fc="white", ec=RULE, lw=0.6, z=1):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                           linewidth=lw, zorder=z))


def txt(x, y, s, size=BODY, weight="normal", color=INK, ha="left", z=3):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color,
            ha=ha, va="baseline", zorder=z)   # y is the baseline


def arrow(gutter, pad=4.0):
    x0, x1 = gutter[0] + pad, gutter[1] - pad
    ax.add_patch(FancyArrow(x0, ARROW_Y, x1 - x0, 0, width=0.9,
                            head_width=3.4, head_length=3.4,
                            length_includes_head=True, facecolor=INK,
                            edgecolor="none", zorder=3))


box(0.4, 0.4, W - 0.8, H - 0.8, fc="#F1F6FC", ec="#C4D6EA", z=0)

# ---------------- stage 1: clinical notes ----------------
txt(C1, 10, "Clinical notes", 7, "bold", BLUE)
box(C1, TOP, C1W, 168)
NOTE = [("Note 1:", "bold"), ("This is an 85 year old", ""),
        ("patient with non-small", ""), ("cell lung cancer here", ""),
        ("for follow-up.", ""), ("Family history:", "bold"),
        ("Mother with breast Ca", ""), ("Social history:", "bold"),
        ("No alcohol, plumber", ""), ("Physical exam:", "bold"),
        ("In no acute distress.", ""), ("Assessment and plan:", "bold"),
        ("Overall, appears to be", ""), ("responding but", ""),
        ("pneumonitis suspected", ""), ("due to pembrolizumab.", "")]
HL = 14                      # first NOTE index to highlight
box(C1 + 1, 22 + HL * LEAD - 7.5, C1W - 2, 2 * LEAD, fc="#FCEFB4", ec="none", z=1)
for i, (line, wt) in enumerate(NOTE):
    txt(C1 + 3, 22 + i * LEAD, line, BODY, wt or "normal")

# ---------------- stage 2: retrieval ----------------
txt(C2, 10, "Retrieval", 7, "bold", BLUE)
CHUNKS = [(TOP, "Chunk 4 \u00b7 cos 0.91", "pneumonitis suspected",
           "due to pembrolizumab", "#FCEFB4", "#D9C36A", INK),
          (54, "Chunk 11 \u00b7 cos 0.87", "bronchoscopy for",
           "suspected pneumonitis", "#FCEFB4", "#D9C36A", INK),
          (94, "Chunk 2 \u00b7 cos 0.24", "Mother with breast Ca",
           "No alcohol, plumber", "#F4F5F7", "#C9CFD8", MUTED)]
for y, head_, l1, l2, fc, ec, col in CHUNKS:
    box(C2, y, C2W, 34, fc=fc, ec=ec)
    txt(C2 + 3, y + 9, head_, 5.5, "bold", col)
    txt(C2 + 3, y + 19, l1, BODY, color=col)
    txt(C2 + 3, y + 29, l2, BODY, color=col)

for i, lab in enumerate(["0\u201390", "91\u2013180", "181\u2013270"]):
    bx = C2 + i * 23
    box(bx, 146, 20, 18, ec="#B8C4D4")
    txt(bx + 10, 157, lab, 5.0, ha="center")
txt(C2, 176, "90-day batches per patient", 5.5)

# ---------------- stage 3: AE calls ----------------
txt(C3, 10, "AE calls", 7, "bold", BLUE)
box(C3, TOP, C3W, 76)
box(C3, TOP, C3W, 16, fc=HDR, ec="none", z=2)
for yy in (30, 60):
    box(C3, yy, C3W, 15, fc=ZEBRA, ec="none", z=1)
XPT, XP, XCALL, XDIV = C3 + 3, 211.0, 226.0, 220.0
ax.plot([XDIV, XDIV], [TOP, 90], color="#DDE3EB", lw=0.5, zorder=2)
for x, s, ha in [(XPT, "Pt", "left"), (XP, "p(AE)", "center"),
                 (XCALL, "Call", "center")]:
    txt(x, 24.5, s, 5.0, "bold", "white", ha, z=3)
ROWS = [("Pt 1", "0.03", "0", ""), ("Pt 2", "0.91", "1", "bold"),
        ("Pt 3", "0.12", "0", ""), ("Pt 4", "0.07", "0", "")]
for i, (pt, p, call, wt) in enumerate(ROWS):
    y = 40.5 + i * 15
    txt(XPT, y, pt, BODY)
    txt(XP, y, p, BODY, wt or "normal", ha="center")
    txt(XCALL, y, call, BODY, wt or "normal", ha="center")

AX0, AX1 = C3 + 3, C3 + C3W - 3
ax.plot([AX0, AX1], [110, 110], color="#8A94A2", lw=0.6, zorder=2)
for x, lab in ((AX0, "0"), (AX1, "1")):
    ax.plot([x, x], [107.5, 112.5], color="#8A94A2", lw=0.6, zorder=2)
    txt(x, 120, lab, 5.5, ha="center")
TAU = AX0 + 0.53 * (AX1 - AX0)
ax.plot([TAU, TAU], [105, 115], color="#B03A2E", lw=0.6, ls=(0, (2, 1.6)), zorder=3)
txt(TAU, 101.5, "\u03c4 (max F1)", 5.5, "bold", "#B03A2E", ha="center")
for x in (AX0 + 2, AX0 + 6, AX0 + 9.5):
    ax.plot(x, 110, "o", ms=2.2, color="#8A94A2", zorder=3)
ax.plot(AX1 - 5, 110, "o", ms=2.5, color=INK, zorder=3)

CAP = ["Score = max across", "90-day intervals", "", "\u03c4 set per AE to",
       "maximise F1 vs", "the chart-reviewed", "gold standard"]
y = 136
for line in CAP:
    if line:
        txt(C3, y, line, 5.5)
    y += 8 if line else 2

# ---------------- stage 4: MSK-Tox dataset (table) ----------------
txt(C4, 10, "MSK-Tox dataset", 7, "bold", BLUE)
AES = ["Pneumonitis", "Colitis", "Hyperthyroidism", "Hypothyroidism",
       "Hepatotoxicity", "Adrenal"]
txt(C4, 28, "N = 55,406 patients", 5.5, "bold")
for i, name in enumerate(AES):
    y = 52 + i * 22
    ax.plot(C4 + 2, y - 2, "o", ms=1.7, color=INK, zorder=3)
    txt(C4 + 7, y, name, BODY)
txt(C4 + 7, 169.5, "insufficiency", BODY)

# ---------------- arrows ----------------
arrow(G1)
txt(sum(G1) / 2, 87, "RAG", 7, "bold", BLUE, ha="center")
arrow(G2)
txt(sum(G2) / 2, 87, "LLM", 7, "bold", BLUE, ha="center")
arrow(G3, pad=3.0)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, facecolor="none")
print(f"wrote {OUT}")
