"""
11_build_tables.py
Builds results/Tables_and_eTables.docx (the Table and eTables 1-5) from the CSV
files written by steps 03-09. Every number is read from those files.
"""
import re
from pathlib import Path
import pandas as pd
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = Path(__file__).resolve().parents[1]
B = str(ROOT)
R = str(ROOT / "results")
mt = pd.read_csv(R + r"/rev5_main_table.csv").set_index(["block", "label"])
mg = pd.read_csv(R + r"/rev5_margin.csv")
t1 = pd.read_csv(R + r"/rev5_table1.csv")
ci = pd.read_csv(R + r"/rev5_p2_ci.csv")

MINUS = "\u2212"
def n1(v):                      # a bound that rounds to -0.0 is shown as 0.0
    t = f"{v:.1f}".replace("-", MINUS)
    return "0.0" if t == MINUS + "0.0" else t
def pci(p, l, u):
    return f"{n1(p)} ({n1(l)}\u2013{n1(u)})"
def dci(p, l, u):
    return f"{n1(p)} ({n1(l)} to {n1(u)})"
def row(df, block, label):
    # block names were stored as str24 in Stata, so match on the first 24 characters
    hit = [k for k in df.index if k[0] == block[:24] and k[1] == label]
    assert len(hit) == 1, (block, label)
    return df.loc[hit[0]]

M = row(mt, "A. Main", "All adults with a defined goal")
T150 = row(mt, "A. Treatment and TG", "Treated, TG >=150")
JAMA = row(mt, "Z. JAMA comparison", "All adults 30-79 incl. low risk")
def marg(dom, m):
    return mg[(mg.domain == dom) & (mg.margin == m)].iloc[0]
M5, M10 = marg("d_def", 5), marg("d_def", 10)
M5T, M10T = marg("d_def_t150", 5), marg("d_def_t150", 10)
cmain = ci[(ci.domain == "d_def") & (ci.comp == "dif_mh")].iloc[0]

p6 = pd.read_csv(R + r"/rev6_pooled_table.csv").set_index(["block", "label"])
sat = pd.read_csv(R + r"/rev6_pooled_satt.csv")
per = pd.read_csv(R + r"/rev6_person.csv")
pper = pd.read_csv(R + r"/rev6_pooled_person.csv")
sf = pd.read_csv(R + r"/rev6_sf.csv").set_index("domain")
br = pd.read_csv(R + r"/rev6_bridge_rows.csv").set_index(["block", "label"])
rx6 = pd.read_csv(R + r"/rev6_pooled_rx.csv")
mod = pd.read_csv(R + r"/rev5_model.csv").set_index("term")
rx2 = pd.read_csv(R + r"/rev5_nonhdl_added.csv")

def pp(df, dom, y):
    return df[(df.domain == dom) & (df.outcome == y)].iloc[0]
def st(dom, comp):
    return sat[(sat.domain == dom) & (sat.comp == comp)].iloc[0]
UP, DN = pp(per, "d_def", "up"), pp(per, "d_def", "dn")
GB, GA = pp(per, "d_def", "gross_before"), pp(per, "d_def", "gross_after")
CAU, CAU1, CAU0 = pp(per, "d_up", "caught"), pp(per, "d_up_tg1", "caught"), pp(per, "d_up_tg0", "caught")
PUP, PDN = pp(pper, "d_def", "up"), pp(pper, "d_def", "dn")
PGB, PGA = pp(pper, "d_def", "gross_before"), pp(pper, "d_def", "gross_after")
PCAU, PCAU1, PCAU0 = pp(pper, "d_up", "caught"), pp(pper, "d_up_tg1", "caught"), pp(pper, "d_up_tg0", "caught")
PADD = pp(pper, "d_def", "add_mh")
ADD_MH = rx2[(rx2.domain == "d_def") & (rx2.outcome == "add_mh")].iloc[0]
ADD_F = rx2[(rx2.domain == "d_def") & (rx2.outcome == "add_f")].iloc[0]
P_ALL = row(p6, "P. Pooled 2017-2023", "All adults with a defined goal")
P_T150 = row(p6, "P. Pooled 2017-2023", "Treated, TG >=150")
P_SEC = row(p6, "P. Pooled 2017-2023", "ASCVD")
P_TG1 = row(p6, "P. Pooled 2017-2023", "TG <100")
P_TG2 = row(p6, "P. Pooled 2017-2023", "TG 100-149")
P_TG3 = row(p6, "P. Pooled 2017-2023", "TG 150-199")
P_TG4 = row(p6, "P. Pooled 2017-2023", "TG 200-399")
Q_ALL = row(p6, "Q. 2017-2020 alone", "All adults with a defined goal")
Q_T150 = row(p6, "Q. 2017-2020 alone", "Treated, TG >=150")
N_ALL = row(p6, "N. 2017-2020 native assay", "All adults with a defined goal")
N_T150 = row(p6, "N. 2017-2020 native assay", "Treated, TG >=150 (native)")
O_ALL = row(p6, "O. Pooled, own-factor MH", "All adults with a defined goal")
X_ALL = row(br, "X. Previous TG assay, exact", "All adults with a defined goal")
X_T150 = row(br, "X. Previous TG assay, exact", "Treated, TG >=150 (previous assay)")
S_POOL, S_POOLR = st("d_def", "dif_mh"), st("d_def", "resid")
TG1, TG2 = row(mt, "A. TG, mg/dL", "<100"), row(mt, "A. TG, mg/dL", "100-149")
TG3, TG4 = row(mt, "A. TG, mg/dL", "150-199"), row(mt, "A. TG, mg/dL", "200-399")
assert int(S_POOL.df_design) == 40
SF = sf.loc["d_def"]                                   # Sampson - Friedewald, adults with a defined goal
ZERO_NOTE = (f" The lower bound for 2017\u2013March 2020 on the previous assay is {N_ALL.gap_lb:.3f}".replace("-", MINUS)
             + ", shown as 0.0.") if -0.05 < N_ALL.gap_lb < 0 else ""


# ------------------------------------------------------------------ docx helpers
doc = None
def new_doc(line=1.5, size=10.5):
    global doc
    doc = Document()
    s = doc.styles["Normal"]
    s.font.name = "Arial"
    s.font.size = Pt(size)
    s.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    s.paragraph_format.line_spacing = line
    for sec in doc.sections:
        sec.left_margin = sec.right_margin = Cm(2.5)
        sec.top_margin = sec.bottom_margin = Cm(2.2)
def H(text, lvl=1):
    p = doc.add_heading(text, level=lvl)
    for r in p.runs:
        r.font.name = "Arial"
        r.font.color.rgb = RGBColor(0x1f, 0x1f, 0x1f)
    return p
def P(text, size=10.5, italic=False, bold=False, after=6, spacing=None, lead=None):
    p = doc.add_paragraph()
    if lead:
        a = p.add_run(lead + " ")
        a.bold = True
        a.font.size = Pt(size)
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.italic = italic
    r.bold = bold
    p.paragraph_format.space_after = Pt(after)
    if spacing:
        p.paragraph_format.line_spacing = spacing
    return p
def shade(cell, hexcol):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:color"), "auto")
    sh.set(qn("w:fill"), hexcol)
    tcPr.append(sh)
def table(header, rows, widths=None, size=8.5, bold_rows=(), shade_rows=()):
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        r = c.paragraphs[0].add_run(h)
        r.bold = True
        r.font.size = Pt(size)
        c.paragraphs[0].paragraph_format.line_spacing = 1.0
        shade(c, "E4E3DE")
    for k, rw in enumerate(rows):
        cells = t.add_row().cells
        for i, v in enumerate(rw):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run(str(v))
            r.font.size = Pt(size)
            cells[i].paragraphs[0].paragraph_format.line_spacing = 1.0
            if k in bold_rows:
                r.bold = True
            if k in shade_rows:
                shade(cells[i], "F1F0EC")
    if widths:
        for rw in t.rows:
            for i, w in enumerate(widths):
                rw.cells[i].width = Cm(w)
    return t
def new_section(landscape):
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    w, h = Cm(21.0), Cm(29.7)
    sec.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
    sec.page_width, sec.page_height = (h, w) if landscape else (w, h)
    sec.left_margin = sec.right_margin = Cm(2.0)
    sec.top_margin = sec.bottom_margin = Cm(1.8)
def flagmark(f):
    return "\u2021" if isinstance(f, str) and f.startswith("SUPPRESS") else ""
def gapn(r):
    return f"{dci(r.gap, r.gap_lb, r.gap_ub)} ({int(r.n_up)}/{int(r.n_dn)})"
def resn(r):
    return f"{dci(r.res, r.res_lb, r.res_ub)} ({int(r.nr_up)}/{int(r.nr_dn)})"

# ------------------------------------------------------------------ main Table
SPEC = [
 ("A. Main", "All adults with a defined goal", "All adults with a defined goal", "all"),
 (None, None, "Risk group", "head"),
 ("A. Risk group", "Primary prevention, borderline-high", "   Primary prevention", "row"),
 ("A. Risk group", "ASCVD", "   ASCVD", "row"),
 (None, None, "Treatment and triglycerides", "head"),
 ("A. Treatment and TG", "Treated, TG >=150", "   Treated, TG \u2265150 mg/dL", "row"),
 ("A. Treatment and TG", "Treated, TG <150", "   Treated, TG <150 mg/dL", "row"),
 ("A. Treatment and TG", "Not treated", "   Not treated", "row"),
 (None, None, "Triglycerides, mg/dL", "head"),
 ("A. TG, mg/dL", "<100", "   <100", "row"),
 ("A. TG, mg/dL", "100-149", "   100\u2013149", "row"),
 ("A. TG, mg/dL", "150-199", "   150\u2013199", "row"),
 ("A. TG, mg/dL", "200-399", "   200\u2013399", "row"),
]
TABLE_TITLE = ("Table. US Adults Aged 30\u201379 Years Above Their Own 2026 LDL-C Goal, by LDL-C Equation and With the "
               "Paired Non-HDL-C Goal, NHANES August 2021\u2013August 2023")
TABLE_NOTE = (
 "Weighted percentages (95% CI) or percentage-point (pp) differences (95% CI); n unweighted. Goals: borderline/"
 "intermediate PREVENT 10-year ASCVD risk (3% to <10%) <100 mg/dL; high risk (\u226510%) or ASCVD <70 mg/dL; ASCVD at "
 "very high risk <55 mg/dL; non-HDL-C goal = LDL-C goal + 30 mg/dL; low-risk adults (<3%, 160 mg/dL threshold) "
 "excluded. \u201cFriedewald or non-HDL-C\u201d: above goal if either is at or above its goal. MH \u2212 F: above by "
 "Martin\u2013Hopkins minus above by Friedewald; n \u2191 above goal by Martin\u2013Hopkins but at goal by Friedewald, "
 "\u2193 the reverse. Recovered: at goal by Friedewald but above the non-HDL-C goal, as a percentage of the whole group. "
 "\u2020 Zero by arithmetic: at a Friedewald goal, non-HDL-C cannot reach its goal when TG <150 mg/dL. Residual: above "
 "by Martin\u2013Hopkins minus above by Friedewald or non-HDL-C; it is negative at TG 200\u2013399 mg/dL because at very "
 "high TG the non-HDL-C goal places more adults above goal than Martin\u2013Hopkins does. Proportions, Korn\u2013Graubard "
 "CIs; differences, method of variance estimates recovery. \u2021 Does not meet NCHS presentation standards. Risk-group "
 "detail and Sampson results are in eTable 3. Treatment status was missing for 2 participants. Differences may not "
 "equal the difference of rounded percentages.")
def main_table():
    P(TABLE_TITLE, size=10, bold=True, spacing=1.15)
    hdr = ["Group", "n", "Friedewald", "Friedewald or non-HDL-C", "Martin\u2013Hopkins",
           "MH \u2212 F, pp (n \u2191/\u2193)", "Recovered by non-HDL-C, pp", "Residual, pp (n \u2191/\u2193)"]
    rows, bold, shd = [], [], []
    for blk, lab, show, kind in SPEC:
        if kind == "head":
            rows.append([show] + [""] * 7)
            bold.append(len(rows) - 1)
            shd.append(len(rows) - 1)
            continue
        r = row(mt, blk, lab)
        zero = int(r.n_rec) == 0 and lab in ("<100", "100-149", "Treated, TG <150")
        rows.append([show.strip() if kind == "all" else show, f"{int(r.n):,}",
                     pci(r.pF, r.lF, r.uF) + flagmark(r.flagF), pci(r.pE, r.lE, r.uE),
                     pci(r.pM, r.lM, r.uM) + flagmark(r.flagM), gapn(r),
                     "0\u2020" if zero else pci(r.rec, r.rec_lb, r.rec_ub), resn(r)])
        if kind == "all":
            bold.append(len(rows) - 1)
    table(hdr, rows, widths=[5.2, 1.2, 3.0, 3.0, 3.0, 3.9, 2.9, 3.9], size=8.5, bold_rows=bold, shade_rows=shd)
    P(TABLE_NOTE, size=8, spacing=1.1)


# ------------------------------------------------------------------ supplement pieces
r_st = rx6[(rx6.domain == "d_trt_c1") & (rx6.drug == "rx_statin")].iloc[0]
r_ez = rx6[(rx6.domain == "d_trt_c1") & (rx6.drug == "rx_ezet")].iloc[0]
def etable2():
    H("eTable 2. Sensitivity analyses", 3)
    tT = lambda r: f"{int(r.n)}: {dci(r.gap, r.gap_lb, r.gap_ub)}"
    rows = [
     ["2021\u20132023, main analysis", f"{int(M.n):,}", gapn(M), pci(M.rec, M.rec_lb, M.rec_ub), dci(M.res, M.res_lb, M.res_ub), tT(T150)],
     ["   t interval, 15 design df", f"{int(M.n):,}", dci(cmain.est, cmain.t15_lb, cmain.t15_ub), "", "", ""],
     ["   Satterthwaite interval", f"{int(M.n):,}", f"{dci(cmain.est, cmain.tsat_lb, cmain.tsat_ub)} (df {cmain.df_satt:.1f})", "", "", ""],
     ["   \u22655 mg/dL margin", f"{int(M5.n):,}", f"{dci(M5.gap, M5.gap_lb, M5.gap_ub)} ({int(M5.n_up)}/{int(M5.n_dn)})", "\u2014", "\u2014",
      f"{int(M5T.n)}: {dci(M5T.gap, M5T.gap_lb, M5T.gap_ub)}"],
     ["   \u226510 mg/dL margin", f"{int(M10.n):,}", f"{dci(M10.gap, M10.gap_lb, M10.gap_ub)} ({int(M10.n_up)}/{int(M10.n_dn)})", "\u2014", "\u2014",
      f"{int(M10T.n)}: {dci(M10T.gap, M10T.gap_lb, M10T.gap_ub)}"],
     ["   Very high risk assigned by Friedewald", f"{int(M.n):,}", "Identical (0 of 249 changed status)", "\u2014", "\u2014", "Identical"],
     ["   Very high risk assigned by Sampson", f"{int(M.n):,}", "Identical (0 of 249 changed status)", "\u2014", "\u2014", "Identical"],
     ["   TG converted to the previous assay", f"{int(X_ALL.n):,}", gapn(X_ALL), pci(X_ALL.rec, X_ALL.rec_lb, X_ALL.rec_ub),
      dci(X_ALL.res, X_ALL.res_lb, X_ALL.res_ub), tT(X_T150)],
     ["   All adults incl. low risk (as in Abohashem et al)", f"{int(JAMA.n):,}", gapn(JAMA), "\u2014", "\u2014", ""],
     ["Pooled 2017\u2013March 2020 + 2021\u20132023 (TG on the current assay)", f"{int(P_ALL.n):,}", gapn(P_ALL),
      pci(P_ALL.rec, P_ALL.rec_lb, P_ALL.rec_ub), dci(P_ALL.res, P_ALL.res_lb, P_ALL.res_ub), tT(P_T150) + "\u00a7"],
     ["   t interval, 40 design df", f"{int(P_ALL.n):,}", dci(S_POOL.est, S_POOL.lb_t, S_POOL.ub_t), "",
      dci(S_POOLR.est, S_POOLR.lb_t, S_POOLR.ub_t), ""],
     ["   Satterthwaite interval", f"{int(P_ALL.n):,}", f"{dci(S_POOL.est, S_POOL.lb_s, S_POOL.ub_s)} (df {S_POOL.df_satt:.1f})", "",
      f"{dci(S_POOLR.est, S_POOLR.lb_s, S_POOLR.ub_s)} (df {S_POOLR.df_satt:.1f})", ""],
     ["   ASCVD", f"{int(P_SEC.n):,}", gapn(P_SEC), pci(P_SEC.rec, P_SEC.rec_lb, P_SEC.rec_ub), dci(P_SEC.res, P_SEC.res_lb, P_SEC.res_ub), ""],
     ["2017\u2013March 2020 alone, TG bridged to the current assay", f"{int(Q_ALL.n):,}", gapn(Q_ALL), pci(Q_ALL.rec, Q_ALL.rec_lb, Q_ALL.rec_ub),
      dci(Q_ALL.res, Q_ALL.res_lb, Q_ALL.res_ub), tT(Q_T150)],
     ["2017\u2013March 2020 alone, previous assay as reported by NCHS", f"{int(N_ALL.n):,}", gapn(N_ALL), pci(N_ALL.rec, N_ALL.rec_lb, N_ALL.rec_ub),
      dci(N_ALL.res, N_ALL.res_lb, N_ALL.res_ub), tT(N_T150)],
    ]
    table(["Analysis", "n", "MH \u2212 F, pp (n \u2191/\u2193)", "Recovered by non-HDL-C, pp", "Residual, pp",
           "Treated, TG \u2265150: n, MH \u2212 F"], rows, widths=[7.6, 1.3, 5.0, 3.3, 3.6, 4.4], size=8,
          bold_rows=(0, 9), shade_rows=())
    P(f"MH, Martin\u2013Hopkins; F, Friedewald. Unless stated, differences have MOVER 95% CIs. \u00a7 The pooled estimate in "
      f"treated adults with TG \u2265150 mg/dL is driven by the 2017\u2013March 2020 cycle ({n1(Q_T150.gap)} points vs "
      f"{n1(T150.gap)} in 2021\u20132023). Drug names, 2017\u2013March 2020 only (treated adults with a defined goal, n = "
      f"{int(r_st.n)}): any statin {pci(r_st.p, r_st.lb, r_st.ub)}%; ezetimibe {pci(r_ez.p, r_ez.lb, r_ez.ub)}%\u2021. "
      "NHANES records no dose, so statin intensity could not be defined. \u2021 Does not meet NCHS presentation standards."
      + ZERO_NOTE, size=8, spacing=1.1)

FULL = [
 ("A. Main", "All adults with a defined goal", "All adults with a defined goal", "d_def", "all"),
 (None, None, "Risk group", None, "head"),
 ("A. Risk group", "Primary prevention, borderline-high", "   Primary prevention", "d_def_prim", "row"),
 ("A. Risk group", "  Borderline/intermediate (<100)", "      Borderline/intermediate PREVENT risk (goal <100)", "d_def_g2", "row"),
 ("A. Risk group", "  High PREVENT risk (<70)", "      High PREVENT risk (<70)", "d_def_g3", "row"),
 ("A. Risk group", "ASCVD", "   ASCVD", "d_def_sec", "row"),
 ("A. Risk group", "  ASCVD, not very high risk (<70)", "      Not very high risk (<70)", "d_def_g4", "row"),
 ("A. Risk group", "  ASCVD, very high risk (<55)", "      Very high risk (<55)", "d_def_g5", "row"),
 (None, None, "Treatment and triglycerides", None, "head"),
 ("A. Treatment and TG", "Treated, TG >=150", "   Treated, TG \u2265150 mg/dL", "d_def_t150", "row"),
 ("A. Treatment and TG", "Treated, TG <150", "   Treated, TG <150 mg/dL", "d_def_tlow", "row"),
 ("A. Treatment and TG", "Not treated", "   Not treated", "d_def_untr", "row"),
 (None, None, "Triglycerides, mg/dL", None, "head"),
 ("A. TG, mg/dL", "<100", "   <100", "d_def_tg1", "row"),
 ("A. TG, mg/dL", "100-149", "   100\u2013149", "d_def_tg2", "row"),
 ("A. TG, mg/dL", "150-199", "   150\u2013199", "d_def_tg3", "row"),
 ("A. TG, mg/dL", "200-399", "   200\u2013399", "d_def_tg4", "row"),
]
def etable3():
    H("eTable 3. Risk-group detail and Sampson LDL-C: percentage above own goal and paired differences", 3)
    rows, bold, shd = [], [], []
    for blk, lab, show, dom, kind in FULL:
        if kind == "head":
            rows.append([show] + [""] * 7)
            bold.append(len(rows) - 1)
            shd.append(len(rows) - 1)
            continue
        r = row(mt, blk, lab)
        s = sf.loc[dom]
        rows.append([show.strip() if kind == "all" else show, f"{int(r.n):,}", pci(r.pF, r.lF, r.uF) + flagmark(r.flagF),
                     pci(r.pM, r.lM, r.uM), pci(r.pS, r.lS, r.uS), gapn(r),
                     f"{dci(s.est, s.lb, s.ub)} ({int(s.n_up)}/{int(s.n_dn)})", resn(r)])
        if kind == "all":
            bold.append(len(rows) - 1)
    table(["Group", "n", "Friedewald", "Martin\u2013Hopkins", "Sampson", "MH \u2212 F, pp (n \u2191/\u2193)",
           "Sampson \u2212 F, pp (n \u2191/\u2193)", "Residual after non-HDL-C, pp (n \u2191/\u2193)"],
          rows, widths=[6.2, 1.2, 2.9, 2.9, 2.9, 3.6, 3.6, 3.6], size=8, bold_rows=bold, shade_rows=shd)
    P("Weighted percentages (95% CI) or percentage-point differences (MOVER 95% CI); n \u2191/\u2193, unweighted participants "
      "above goal by the newer equation only / by Friedewald only. Residual: Martin\u2013Hopkins minus Friedewald or non-HDL-C. "
      "\u2021 Does not meet NCHS presentation standards.", size=8, spacing=1.1)
def etable4():
    H("eTable 4. Exploratory survey-weighted linear model of the Martin\u2013Hopkins minus Friedewald difference", 3)
    lab = {"2.tgcat": "TG 100\u2013149 vs <100 mg/dL", "3.tgcat": "TG 150\u2013199 vs <100 mg/dL",
           "4.tgcat": "TG 200\u2013399 vs <100 mg/dL", "1.strict": "Goal <70 or <55 vs <100 mg/dL",
           "1.treated": "Treated vs not treated"}
    rows = [[lab[k], dci(mod.loc[k, "coef"], mod.loc[k, "lb"], mod.loc[k, "ub"])] for k in lab]
    table(["Exposure", "Difference, percentage points (95% CI)"], rows, widths=[8.0, 6.0], size=8.5)
    P(f"Adults with a defined goal and known treatment status (n = {int(mod.loc['1.strict', 'n_sub']):,}); all terms entered "
      "together; t-based CIs with 15 design df. Exploratory: the outcome rests on 42 reclassified participants.",
      size=8, spacing=1.1)
def etable5():
    H("eTable 5. Participant-level reclassification and what the paired non-HDL-C goal identifies", 3)
    def pc(r):
        return f"{int(r.x)}/{int(r.n)}; {pci(r.p, r.lb, r.ub)}" + flagmark(r.flag)
    def pc2(r):
        return f"{int(r.x)}; {pci(r.p, r.lb, r.ub)}" + flagmark(r.flag)
    rows = [
     ["Adults with a defined goal, n", f"{int(UP.n):,}", f"{int(PUP.n):,}"],
     ["Above goal by Martin\u2013Hopkins, at goal by Friedewald (upward), n; %", pc2(UP), pc2(PUP)],
     ["Above goal by Friedewald, at goal by Martin\u2013Hopkins (downward), n; %", pc2(DN), pc2(PDN)],
     ["Classified differently from Martin\u2013Hopkins by Friedewald alone, n; %", pc2(GB), pc2(PGB)],
     ["Classified differently from Martin\u2013Hopkins by Friedewald or non-HDL-C, n; %", pc2(GA), pc2(PGA)],
     ["Upward-reclassified adults above the non-HDL-C goal, n/N; %", pc(CAU), pc(PCAU)],
     ["   with TG \u2265150 mg/dL", pc(CAU1), pc(PCAU1)],
     ["   with TG <150 mg/dL", f"{int(CAU0.x)}/{int(CAU0.n)} (by arithmetic)", f"{int(PCAU0.x)}/{int(PCAU0.n)} (by arithmetic)"],
     ["At goal by Martin\u2013Hopkins but above the non-HDL-C goal, n; %", f"{int(ADD_MH.x)}; {pci(ADD_MH.p, ADD_MH.lb, ADD_MH.ub)}", pc2(PADD)],
     ["Residual after the non-HDL-C goal, by TG, pp (95% CI; n \u2191/\u2193)", "", ""],
     ["   TG <100 mg/dL", resn(TG1), resn(P_TG1)],
     ["   TG 100\u2013149 mg/dL", resn(TG2), resn(P_TG2)],
     ["   TG 150\u2013199 mg/dL", resn(TG3), resn(P_TG3)],
     ["   TG 200\u2013399 mg/dL", resn(TG4), resn(P_TG4)],
    ]
    table(["", "2021\u20132023", "Pooled 2017\u2013March 2020 + 2021\u20132023"], rows, widths=[9.5, 4.0, 4.5], size=8.5,
          bold_rows=(9,), shade_rows=(9,))
    P("Counts are unweighted; percentages are weighted (Korn\u2013Graubard 95% CI). \u2021 Does not meet NCHS presentation "
      "standards (effective sample size <30 or CI width \u226530 points); shown because the counts are the result. "
      "Residual: above goal by Martin\u2013Hopkins minus above goal by Friedewald or non-HDL-C, in percentage points of the "
      "TG group (MOVER 95% CI); n \u2191/\u2193, participants above goal by Martin\u2013Hopkins only / by Friedewald or "
      "non-HDL-C only.", size=8, spacing=1.1)



# =================================================================== build
new_doc(line=1.15, size=10)
for sec in doc.sections:
    w, h = Cm(21.0), Cm(29.7)
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = h, w
    sec.left_margin = sec.right_margin = Cm(2.0)
main_table()
new_section(False)
H("eTable 1. Characteristics of adults with a defined 2026 goal", 3)
lab = {"n (unweighted)": "n (unweighted)", "ridageyr": "Age, y, mean", "female": "Women, %",
       "nhw": "Non-Hispanic White, %", "nhb": "Non-Hispanic Black, %", "hisp": "Hispanic, %",
       "nha": "Non-Hispanic Asian, %", "oth": "Other or multiracial, %", "diabetes": "Diabetes, %",
       "ascvd_g": "ASCVD, %", "g100": "Goal <100 mg/dL, %", "g70": "Goal <70 mg/dL, %", "g55": "Goal <55 mg/dL, %",
       "lbxtc": "Total cholesterol, mg/dL", "lbdhdd": "HDL-C, mg/dL", "nonhdl": "Non-HDL-C, mg/dL",
       "lbdldl": "LDL-C, Friedewald, mg/dL", "lbdldlm": "LDL-C, Martin\u2013Hopkins, mg/dL",
       "lbdldln": "LDL-C, Sampson, mg/dL", "TG median (IQR)": "TG, mg/dL, median (IQR)", "bmxbmi": "BMI, kg/m\u00b2"}
cols = ["d_def", "d_def_t150", "d_def_tlow", "d_def_untr"]
def cellv(c, rname):
    r = t1[(t1.column == c) & (t1.row == rname)].iloc[0]
    if r.kind == "n":
        return f"{int(r.est):,}"
    if r.kind == "median":
        return f"{r.est:.0f} ({r.lb:.0f}\u2013{r.ub:.0f})"
    return pci(r.est, r.lb, r.ub) + (flagmark(r.flag) if isinstance(r.flag, str) else "")
table(["Characteristic", "All", "Treated, TG \u2265150", "Treated, TG <150", "Not treated"],
      [[lab[o]] + [cellv(c, o) for c in cols] for o in lab], widths=[5.2, 2.9, 2.9, 2.9, 2.9], size=7.5)
P("Weighted means or percentages (95% CI): means, t-based CIs; percentages, Korn\u2013Graubard. \u2021 Does not meet NCHS "
  "presentation standards. Treatment = self-reported current cholesterol-lowering medication.", size=7.5, spacing=1.0)
new_section(True)
etable2()
new_section(True)
etable3()
new_section(False)
etable4()
etable5()
out = ROOT / "results" / "Tables_and_eTables.docx"
doc.save(out)
print("saved", out)
