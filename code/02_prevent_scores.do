*=======================================================================*
* 02_prevent_scores.do
* AHA PREVENT 10-year ASCVD risk (base equations) for NHANES August
* 2021-August 2023 adults aged 30-79 years.
*
* EQUATIONS
*   Coefficients from the preventr R package (github.com/martingmayer/
*   preventr), which implements Khan SS et al. Circulation 2024;149:430-449.
*   The do-file checks them against the paper's worked example (Table S25)
*   and stops if the result differs.
*   Models: r10_ascvd (base, 10-year ASCVD; used in the paper), r10_cvd,
*   r30_ascvd, r10_ascvd_full (UACR + HbA1c; SDI not in public NHANES).
*
* INPUT DEFINITIONS
*   SBP          mean of available oscillometric readings (BPXO_L)
*   BP treatment bpq150 == 1 (asked only if told high BP; others = 0)
*   Statin       bpq101d == 1, any cholesterol-lowering medication
*                (RXQ_RX_L 2021-2023 has no drug names)
*   Diabetes     told diabetes, insulin, diabetes pills, or HbA1c >= 6.5
*   Smoking      smoked >= 100 cigarettes and now smokes every/some days
*   eGFR         CKD-EPI 2021 (race-free) from serum creatinine
*
* OUT-OF-RANGE VALUES
*   PREVENT valid ranges: TC 130-320, HDL 20-100, SBP 90-180, eGFR 15-140,
*   BMI 18.5-39.9, UACR 0.1-25000, HbA1c 4.5-15. Values are truncated to
*   the range; s_strict flags all required inputs inside the range.
*
* Risk categories (2026 ACC/AHA dyslipidemia guideline): low <3%,
*   borderline 3-<5%, intermediate 5-<10%, high >=10%.
* Inputs (.xpt, data/nhanes_2021_2023): DEMO_L TCHOL_L HDL_L BPXO_L BPQ_L
*   DIQ_L GHB_L BIOPRO_L ALB_CR_L SMQ_L MCQ_L BMX_L
* Output: nhanes_prevent.dta (one row per participant) in the data folder
*=======================================================================*

version 16
clear all
set more off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
cd "$DATA2123"

capture log close
log using "prevent_build_log.smcl", replace


*=======================================================================*
* PART 0. CHECK FILES
*=======================================================================*
local missing ""
foreach f in DEMO_L TCHOL_L HDL_L BPXO_L BPQ_L DIQ_L GHB_L BIOPRO_L ///
             ALB_CR_L SMQ_L MCQ_L BMX_L {
    capture confirm file "`f'.xpt"
    if _rc local missing "`missing' `f'"
}
if "`missing'" != "" {
    display as error "Files NOT FOUND:`missing'"
    log close
    exit 601
}


*=======================================================================*
* PART 1. PREVENT CALCULATOR (program)
*   Needs variables: female (0/1) p_age p_tc p_hdl p_sbp p_bptx p_statin
*                    p_dm p_smk p_egfr p_bmi p_uacr p_a1c
*   (p_uacr / p_a1c may be missing -> missing indicators in full model)
*   Cholesterol in mg/dL; converted to mmol/L with 0.02586 as in preventr.
*=======================================================================*
capture program drop prevent_calc
program define prevent_calc
    capture drop t_* lp_* r10_ascvd r10_cvd r30_ascvd r10_ascvd_full

    * --- centred / spline terms (Khan 2024, Table S12) ---
    gen double t_age      = (p_age - 55) / 10                if p_age < .
    gen double t_age2     = t_age^2
    gen double t_nhdl     = (p_tc - p_hdl) * 0.02586 - 3.5   if p_tc < . & p_hdl < .
    gen double t_hdl      = (p_hdl * 0.02586 - 1.3) / 0.3    if p_hdl < .
    gen double t_sbp_lo   = (min(p_sbp, 110) - 110) / 20     if p_sbp < .
    gen double t_sbp_hi   = (max(p_sbp, 110) - 130) / 20     if p_sbp < .
    gen double t_dm       = p_dm
    gen double t_smk      = p_smk
    gen double t_bmi_lo   = (min(p_bmi, 30) - 25) / 5        if p_bmi < .
    gen double t_bmi_hi   = (max(p_bmi, 30) - 30) / 5        if p_bmi < .
    gen double t_egfr_lo  = (min(p_egfr, 60) - 60) / -15     if p_egfr < .
    gen double t_egfr_hi  = (max(p_egfr, 60) - 90) / -15     if p_egfr < .
    gen double t_bptx     = p_bptx
    gen double t_statin   = p_statin
    gen double t_bptx_sbp    = t_bptx * t_sbp_hi
    gen double t_statin_nhdl = t_statin * t_nhdl
    gen double t_age_nhdl = t_age * t_nhdl
    gen double t_age_hdl  = t_age * t_hdl
    gen double t_age_sbp  = t_age * t_sbp_hi
    gen double t_age_dm   = t_age * t_dm
    gen double t_age_smk  = t_age * t_smk
    gen double t_age_bmi  = t_age * t_bmi_hi
    gen double t_age_egfr = t_age * t_egfr_lo
    * full-model extras (SDI needs ZIP code -> always "missing" in NHANES)
    gen double t_sdi46    = 0
    gen double t_sdi710   = 0
    gen double t_sdimiss  = 1
    gen double t_lnuacr   = cond(p_uacr < ., ln(p_uacr), 0)
    gen double t_uacrmiss = (p_uacr >= .)
    gen double t_a1c_dm   = cond(p_a1c < . & p_dm == 1, p_a1c - 5.3, 0)
    gen double t_a1c_nodm = cond(p_a1c < . & p_dm == 0, p_a1c - 5.3, 0)
    gen double t_a1cmiss  = (p_a1c >= .)

    * --- linear predictors (coefficients generated from preventr) ---
    * base_10yr - ascvd
    gen double lp_b10_ascvd = .
    replace lp_b10_ascvd = -3.819975 + (0.719883*t_age) + (0.1176967*t_nhdl) ///
        + (-0.151185*t_hdl) + (-0.0835358*t_sbp_lo) + (0.3592852*t_sbp_hi) ///
        + (0.8348585*t_dm) + (0.4831078*t_smk) + (0.4864619*t_egfr_lo) ///
        + (0.0397779*t_egfr_hi) + (0.2265309*t_bptx) + (-0.0592374*t_statin) ///
        + (-0.0395762*t_bptx_sbp) + (0.0844423*t_statin_nhdl) ///
        + (-0.0567839*t_age_nhdl) + (0.0325692*t_age_hdl) + (-0.1035985*t_age_sbp) ///
        + (-0.2417542*t_age_dm) + (-0.0791142*t_age_smk) + (-0.1671492*t_age_egfr) if female == 1
    replace lp_b10_ascvd = -3.500655 + (0.7099847*t_age) + (0.1658663*t_nhdl) ///
        + (-0.1144285*t_hdl) + (-0.2837212*t_sbp_lo) + (0.3239977*t_sbp_hi) ///
        + (0.7189597*t_dm) + (0.3956973*t_smk) + (0.3690075*t_egfr_lo) ///
        + (0.0203619*t_egfr_hi) + (0.2036522*t_bptx) + (-0.0865581*t_statin) ///
        + (-0.0322916*t_bptx_sbp) + (0.114563*t_statin_nhdl) ///
        + (-0.0300005*t_age_nhdl) + (0.0232747*t_age_hdl) + (-0.0927024*t_age_sbp) ///
        + (-0.2018525*t_age_dm) + (-0.0970527*t_age_smk) + (-0.1217081*t_age_egfr) if female == 0

    * base_10yr - total_cvd
    gen double lp_b10_cvd = .
    replace lp_b10_cvd = -3.307728 + (0.7939329*t_age) + (0.0305239*t_nhdl) ///
        + (-0.1606857*t_hdl) + (-0.2394003*t_sbp_lo) + (0.3600781*t_sbp_hi) ///
        + (0.8667604*t_dm) + (0.5360739*t_smk) + (0.6045917*t_egfr_lo) ///
        + (0.0433769*t_egfr_hi) + (0.3151672*t_bptx) + (-0.1477655*t_statin) ///
        + (-0.0663612*t_bptx_sbp) + (0.1197879*t_statin_nhdl) ///
        + (-0.0819715*t_age_nhdl) + (0.0306769*t_age_hdl) + (-0.0946348*t_age_sbp) ///
        + (-0.27057*t_age_dm) + (-0.078715*t_age_smk) + (-0.1637806*t_age_egfr) if female == 1
    replace lp_b10_cvd = -3.031168 + (0.7688528*t_age) + (0.0736174*t_nhdl) ///
        + (-0.0954431*t_hdl) + (-0.4347345*t_sbp_lo) + (0.3362658*t_sbp_hi) ///
        + (0.7692857*t_dm) + (0.4386871*t_smk) + (0.5378979*t_egfr_lo) ///
        + (0.0164827*t_egfr_hi) + (0.288879*t_bptx) + (-0.1337349*t_statin) ///
        + (-0.0475924*t_bptx_sbp) + (0.150273*t_statin_nhdl) ///
        + (-0.0517874*t_age_nhdl) + (0.0191169*t_age_hdl) + (-0.1049477*t_age_sbp) ///
        + (-0.2251948*t_age_dm) + (-0.0895067*t_age_smk) + (-0.1543702*t_age_egfr) if female == 0

    * base_30yr - ascvd
    gen double lp_b30_ascvd = .
    replace lp_b30_ascvd = -1.974074 + (0.4669202*t_age) + (-0.0893118*t_age2) ///
        + (0.1256901*t_nhdl) + (-0.1542255*t_hdl) + (-0.0018093*t_sbp_lo) ///
        + (0.322949*t_sbp_hi) + (0.6296707*t_dm) + (0.268292*t_smk) ///
        + (0.100106*t_egfr_lo) + (0.0499663*t_egfr_hi) + (0.1875292*t_bptx) ///
        + (0.0152476*t_statin) + (-0.0276123*t_bptx_sbp) ///
        + (0.0736147*t_statin_nhdl) + (-0.0521962*t_age_nhdl) ///
        + (0.0316918*t_age_hdl) + (-0.1046101*t_age_sbp) + (-0.2727793*t_age_dm) ///
        + (-0.1530907*t_age_smk) + (-0.1299149*t_age_egfr) if female == 1
    replace lp_b30_ascvd = -1.736444 + (0.3994099*t_age) + (-0.0937484*t_age2) ///
        + (0.1744643*t_nhdl) + (-0.120203*t_hdl) + (-0.0665117*t_sbp_lo) ///
        + (0.2753037*t_sbp_hi) + (0.4790257*t_dm) + (0.1782635*t_smk) ///
        + (-0.0218789*t_egfr_lo) + (0.0602553*t_egfr_hi) + (0.1421182*t_bptx) ///
        + (0.0135996*t_statin) + (-0.0218265*t_bptx_sbp) ///
        + (0.1013148*t_statin_nhdl) + (-0.0312619*t_age_nhdl) ///
        + (0.020673*t_age_hdl) + (-0.0920935*t_age_sbp) + (-0.2159947*t_age_dm) ///
        + (-0.1548811*t_age_smk) + (-0.0712547*t_age_egfr) if female == 0

    * full_10yr - ascvd
    gen double lp_f10_ascvd = .
    replace lp_f10_ascvd = -4.291503 + (0.7023067*t_age) + (0.0898765*t_nhdl) ///
        + (-0.1407316*t_hdl) + (-0.0256648*t_sbp_lo) + (0.314511*t_sbp_hi) ///
        + (0.4799217*t_dm) + (0.4062049*t_smk) + (0.3847744*t_egfr_lo) ///
        + (0.0495174*t_egfr_hi) + (0.2133861*t_bptx) + (-0.0678552*t_statin) ///
        + (-0.0451416*t_bptx_sbp) + (0.0788187*t_statin_nhdl) ///
        + (-0.0535985*t_age_nhdl) + (0.0291762*t_age_hdl) + (-0.0961839*t_age_sbp) ///
        + (-0.2001466*t_age_dm) + (-0.0586472*t_age_smk) + (-0.1537791*t_age_egfr) ///
        + (0.1413965*t_sdi46) + (0.228136*t_sdi710) + (0.1588908*t_sdimiss) ///
        + (0.1371824*t_lnuacr) + (0.0061613*t_uacrmiss) + (0.123192*t_a1c_dm) ///
        + (0.1410572*t_a1c_nodm) + (0.005866*t_a1cmiss) if female == 1
    replace lp_f10_ascvd = -3.969788 + (0.7128741*t_age) + (0.1465201*t_nhdl) ///
        + (-0.1125794*t_hdl) + (-0.3387216*t_sbp_lo) + (0.2980252*t_sbp_hi) ///
        + (0.399583*t_dm) + (0.3379111*t_smk) + (0.2582604*t_egfr_lo) ///
        + (0.0147769*t_egfr_hi) + (0.1686621*t_bptx) + (-0.1073619*t_statin) ///
        + (-0.0381038*t_bptx_sbp) + (0.1034169*t_statin_nhdl) ///
        + (-0.0228755*t_age_nhdl) + (0.0267453*t_age_hdl) + (-0.0897449*t_age_sbp) ///
        + (-0.1497464*t_age_dm) + (-0.077206*t_age_smk) + (-0.1198368*t_age_egfr) ///
        + (0.0651121*t_sdi46) + (0.2676683*t_sdi710) + (0.1388492*t_sdimiss) ///
        + (0.1375837*t_lnuacr) + (0.0652944*t_uacrmiss) + (0.101282*t_a1c_dm) ///
        + (0.1092726*t_a1c_nodm) + (-0.0112852*t_a1cmiss) if female == 0

    * --- risks in % ---
    gen double r10_ascvd      = 100 * invlogit(lp_b10_ascvd)
    gen double r10_cvd        = 100 * invlogit(lp_b10_cvd)
    gen double r30_ascvd      = 100 * invlogit(lp_b30_ascvd)
    gen double r10_ascvd_full = 100 * invlogit(lp_f10_ascvd)
    drop t_* lp_*
end


*=======================================================================*
* PART 2. SELF-TEST against Khan 2024 Table S25 / supplement
*   Woman 50 y, SBP 160 treated, TC 200, HDL 45, no statin, diabetes,
*   non-smoker, eGFR 90, BMI 35
*   Expected: 10-yr ASCVD 9.2%, 10-yr total CVD 14.7%, 30-yr ASCVD 35.4%
*=======================================================================*
preserve
clear
set obs 1
gen female = 1
gen p_age = 50
gen p_sbp = 160
gen p_bptx = 1
gen p_tc = 200
gen p_hdl = 45
gen p_statin = 0
gen p_dm = 1
gen p_smk = 0
gen p_egfr = 90
gen p_bmi = 35
gen p_uacr = .
gen p_a1c = .
prevent_calc
list r10_ascvd r10_cvd r30_ascvd, noobs
assert abs(r10_ascvd - 9.2)  < 0.06
assert abs(r10_cvd   - 14.7) < 0.06
assert abs(r30_ascvd - 35.4) < 0.06
display as result "SELF-TEST PASSED - equations reproduce the published example"
restore


*=======================================================================*
* PART 3. IMPORT AND MERGE
*=======================================================================*
tempfile t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11
import sasxport5 "TCHOL_L.xpt", clear
keep seqn lbxtc
save `t1'
import sasxport5 "HDL_L.xpt", clear
keep seqn lbdhdd
save `t2'
import sasxport5 "BPXO_L.xpt", clear
keep seqn bpxosy1 bpxosy2 bpxosy3
save `t3'
import sasxport5 "BPQ_L.xpt", clear
keep seqn bpq020 bpq150 bpq080 bpq101d
save `t4'
import sasxport5 "DIQ_L.xpt", clear
keep seqn diq010 diq050 diq070
save `t5'
import sasxport5 "GHB_L.xpt", clear
keep seqn lbxgh
save `t6'
import sasxport5 "BIOPRO_L.xpt", clear
keep seqn lbxscr
save `t7'
import sasxport5 "ALB_CR_L.xpt", clear
keep seqn urdact
save `t8'
import sasxport5 "SMQ_L.xpt", clear
keep seqn smq020 smq040
save `t9'
import sasxport5 "MCQ_L.xpt", clear
keep seqn mcq160b mcq160c mcq160d mcq160e mcq160f
save `t10'
import sasxport5 "BMX_L.xpt", clear
keep seqn bmxbmi bmxwaist
save `t11'

import sasxport5 "DEMO_L.xpt", clear
keep seqn riagendr ridageyr ridreth3 indfmpir ridexprg sdmvpsu sdmvstra wtmec2yr
assert _N == 11933
foreach t in t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 {
    merge 1:1 seqn using ``t'', keep(master match) nogen
}
assert _N == 11933


*=======================================================================*
* PART 4. INPUT VARIABLES (raw)
*=======================================================================*
label define yn 0 "No" 1 "Yes", replace

gen female = (riagendr == 2)
recode ridreth3 (1=1 "Mexican American") (2=2 "Other Hispanic") (3=3 "NH White") ///
    (4=4 "NH Black") (6=5 "NH Asian") (7=6 "Other/Multiracial"), gen(race6)
gen agecat = 1 if ridageyr >= 30 & ridageyr < 40
replace agecat = 2 if ridageyr >= 40 & ridageyr < 50
replace agecat = 3 if ridageyr >= 50 & ridageyr < 60
replace agecat = 4 if ridageyr >= 60 & ridageyr < 70
replace agecat = 5 if ridageyr >= 70 & ridageyr < 80
label define agel 1 "30-39" 2 "40-49" 3 "50-59" 4 "60-69" 5 "70-79", replace
label values agecat agel

* SBP: mean of available readings
egen sbp = rowmean(bpxosy1 bpxosy2 bpxosy3)

* BP treatment (bpq150 asked only of those told they have high BP)
gen bptx = .
replace bptx = 0 if bpq020 == 2
replace bptx = 1 if bpq020 == 1 & bpq150 == 1
replace bptx = 0 if bpq020 == 1 & bpq150 == 2

* "Statin" = any cholesterol-lowering medication (proxy)
gen statin = (bpq101d == 1) if inlist(bpq101d, 1, 2)

* Diabetes
gen dm = (diq010 == 1) if inlist(diq010, 1, 2, 3)
replace dm = 1 if diq050 == 1 | diq070 == 1 | (lbxgh >= 6.5 & lbxgh < .)

* Current cigarette smoking
gen smoker = 0 if inlist(smq020, 1, 2)
replace smoker = 1 if smq020 == 1 & inlist(smq040, 1, 2)

* eGFR, CKD-EPI 2021 race-free
gen double kappa = cond(female == 1, 0.7, 0.9)
gen double alpha = cond(female == 1, -0.241, -0.302)
gen double egfr = 142 * (min(lbxscr/kappa, 1))^alpha * (max(lbxscr/kappa, 1))^-1.200 ///
    * 0.9938^ridageyr * cond(female == 1, 1.012, 1) if lbxscr < .
drop kappa alpha

* Prior cardiovascular disease (CHF, CHD, angina, MI, stroke)
gen cvd_hx = .
replace cvd_hx = 0 if mcq160b == 2 & mcq160c == 2 & mcq160d == 2 & mcq160e == 2 & mcq160f == 2
replace cvd_hx = 1 if mcq160b == 1 | mcq160c == 1 | mcq160d == 1 | mcq160e == 1 | mcq160f == 1

label values bptx statin dm smoker cvd_hx yn


*=======================================================================*
* PART 5. SAMPLE FLOW
*=======================================================================*
gen s_age     = (ridageyr >= 30 & ridageyr <= 79 & ridexprg != 1)
gen s_nocvd   = (s_age & cvd_hx == 0)
gen s_complete = (s_nocvd & lbxtc < . & lbdhdd < . & sbp < . & bptx < . & ///
                  statin < . & dm < . & smoker < . & egfr < . & wtmec2yr > 0 & wtmec2yr < .)

display _n "=== SAMPLE FLOW ==="
count
count if s_age
display "  excluded, prior CVD:"
count if s_age & cvd_hx == 1
display "  excluded, CVD history missing:"
count if s_age & cvd_hx == .
count if s_nocvd
display "  missing inputs (each):"
foreach v in lbxtc lbdhdd sbp bptx statin dm smoker egfr {
    quietly count if s_nocvd & missing(`v')
    display "    `v': " r(N)
}
display "  FINAL analytic sample:"
count if s_complete


*=======================================================================*
* PART 6. TRUNCATE TO VALID RANGES AND CALCULATE
*=======================================================================*
gen oor = (lbxtc < 130 | lbxtc > 320 | lbdhdd < 20 | lbdhdd > 100 | ///
           sbp < 90 | sbp > 180 | egfr < 15 | egfr > 140) if s_complete
display "Outside PREVENT range (truncated in main analysis):"
tab oor if s_complete

gen p_age    = ridageyr
gen p_tc     = min(max(lbxtc, 130), 320)        if lbxtc < .
gen p_hdl    = min(max(lbdhdd, 20), 100)        if lbdhdd < .
gen p_sbp    = min(max(sbp, 90), 180)           if sbp < .
gen p_egfr   = min(max(egfr, 15), 140)          if egfr < .
gen p_bmi    = min(max(bmxbmi, 18.5), 39.9)     if bmxbmi < .
gen p_uacr   = min(max(urdact, 0.1), 25000)     if urdact < .
gen p_a1c    = min(max(lbxgh, 4.5), 15)         if lbxgh < .
gen p_bptx   = bptx
gen p_statin = statin
gen p_dm     = dm
gen p_smk    = smoker

prevent_calc
foreach v in r10_ascvd r10_cvd r30_ascvd r10_ascvd_full {
    replace `v' = . if !s_complete
}
* 30-year equations were developed for ages 30-59 only (preventr warns that
* estimates above 59 are "questionable") -> set to missing above 59
replace r30_ascvd = . if ridageyr > 59
drop p_*

label variable r10_ascvd      "PREVENT base 10-yr ASCVD risk (%)"
label variable r10_cvd        "PREVENT base 10-yr total CVD risk (%)"
label variable r30_ascvd      "PREVENT base 30-yr ASCVD risk (%), ages 30-59 only"
label variable r10_ascvd_full "PREVENT full 10-yr ASCVD risk (%), UACR+HbA1c, SDI missing"

* 2026 ACC/AHA categories
label define catl 1 "Low <3%" 2 "Borderline 3-<5%" 3 "Intermediate 5-<10%" 4 "High >=10%", replace
foreach m in r10_ascvd r10_ascvd_full {
    gen cat_`m' = 1 if `m' < 3
    replace cat_`m' = 2 if `m' >= 3 & `m' < 5
    replace cat_`m' = 3 if `m' >= 5 & `m' < 10
    replace cat_`m' = 4 if `m' >= 10 & `m' < .
    label values cat_`m' catl
}
rename cat_r10_ascvd      prevent_cat
rename cat_r10_ascvd_full prevent_cat_full

* Conditions that put people on a treatment pathway regardless of score
gen ckd3 = (egfr < 60) if egfr < .
label values ckd3 yn
label variable ckd3 "eGFR < 60 (CKD G3+) - check guideline: <60 or <=60"

gen s_strict = (s_complete & oor == 0)
label variable s_complete "PREVENT sample: 30-79, no prior CVD, complete inputs"
label variable s_strict   "PREVENT sample, all inputs inside valid range"

summarize r10_ascvd r10_cvd r30_ascvd r10_ascvd_full if s_complete, detail

save nhanes_prevent, replace


*=======================================================================*
* PART 7. DESCRIPTIVE RESULTS (survey-weighted, MEC weight)
*=======================================================================*
svyset sdmvpsu [pweight = wtmec2yr], strata(sdmvstra) vce(linearized)

gen one = 1
svy, subpop(s_complete): total one
display "Weighted population (millions): " %6.1f _b[one]/1e6

svy, subpop(s_complete): mean r10_ascvd r10_ascvd_full r30_ascvd
_pctile r10_ascvd [pweight = wtmec2yr] if s_complete, p(25 50 75)
display "10-yr ASCVD: P25 " %4.1f r(r1) "  median " %4.1f r(r2) "  P75 " %4.1f r(r3)

* Category distribution
svy, subpop(s_complete): tab prevent_cat, percent ci
svy, subpop(s_complete): tab prevent_cat_full, percent ci
svy, subpop(s_complete): tab prevent_cat female, col percent
svy, subpop(s_complete): tab prevent_cat agecat, col percent

* Base vs full model reclassification
tab prevent_cat prevent_cat_full if s_complete
svy, subpop(s_complete): tab prevent_cat prevent_cat_full, row percent

* Who is on cholesterol medication in each category (treatment teaser)
svy, subpop(s_complete): tab prevent_cat statin, row percent

* Automatic-pathway conditions inside each category
svy, subpop(s_complete): tab prevent_cat dm,   row percent
svy, subpop(s_complete): tab prevent_cat ckd3, row percent

* 30-year risk among those at low 10-year risk
gen s_low = s_complete & prevent_cat == 1
svy, subpop(s_low): mean r30_ascvd

* Sensitivity: strict in-range sample
svy, subpop(s_strict): tab prevent_cat, percent


*=======================================================================*
* PART 8. CATEGORY TABLE -> prevent_categories.csv
*=======================================================================*
tempname H
postfile `H' str16 model str24 category n pct lb ub millions using prevent_categories, replace
foreach m in prevent_cat prevent_cat_full {
    forvalues k = 1/4 {
        capture drop ind
        gen ind = (`m' == `k') if `m' < .
        quietly count if s_complete & `m' == `k'
        local n = r(N)
        quietly svy, subpop(s_complete): proportion ind, citype(logit)
        matrix T = r(table)
        local c = colsof(T)
        local p  = 100*T[1,`c']
        local lb = 100*T[5,`c']
        local ub = 100*T[6,`c']
        quietly svy, subpop(s_complete): total ind
        local mil = _b[ind]/1e6
        local lab : label catl `k'
        post `H' ("`m'") ("`lab'") (`n') (`p') (`lb') (`ub') (`mil')
    }
}
postclose `H'
drop ind one

preserve
use prevent_categories, clear
format pct lb ub millions %5.1f
list, noobs clean
export delimited using "prevent_categories.csv", replace
restore


log close
display "DONE - nhanes_prevent.dta, prevent_categories.csv, prevent_build_log.smcl"
*=======================================================================*
* END
*=======================================================================*
