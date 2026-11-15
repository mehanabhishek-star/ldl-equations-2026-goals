*=======================================================================*
* 05_main_table.do
* Adults aged 30-79 years with a defined 2026 goal (<100, <70 or <55 mg/dL):
*   A. Main table: % above goal by Friedewald, Friedewald or paired
*      non-HDL-C, Martin-Hopkins and Sampson; Martin-Hopkins minus
*      Friedewald with participants reclassified up / down; what the
*      non-HDL-C goal adds; residual (MOVER-KG intervals)
*   B. Sensitivity: very-high-risk criterion assigned by Friedewald
*   C. Sensitivity: 5 and 10 mg/dL margin rule
*   D. Sensitivity: TG on the previous assay (NCHS backward equation)
*   E. Table of characteristics with 95% CIs
*   F. Sample reconciliation with Abohashem et al
* Uses _kgprop.do and _moverkg.do.
* Outputs: results/rev5_main_table.csv, rev5_margin.csv, rev5_table1.csv,
*          rev5_flow.csv
*=======================================================================*
version 17
clear all
set more off
set varabbrev off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
local out "$RESULTS"
local code "$CODE"
capture log close
log using "`out'/rev5_main.log", text replace
do "`code'/_kgprop.do"
do "`code'/_moverkg.do"

use "`out'/rev3_analysis_data.dta", clear
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)

*--- domains
gen byte d_def      = (d_bm_all == 1 & formal_goal == 1)
gen byte d_def_prim = (d_def == 1 & inlist(goalgrp, 2, 3))
gen byte d_def_sec  = (d_def == 1 & inlist(goalgrp, 4, 5))
gen byte d_def_t150 = (d_def == 1 & treated == 1 & tg150 == 1)
gen byte d_def_tlow = (d_def == 1 & treated == 1 & tg150 == 0)
gen byte d_def_untr = (d_def == 1 & treated == 0)
forvalues k = 2/5 {
    gen byte d_def_g`k' = (d_def == 1 & goalgrp == `k')
}
forvalues c = 1/4 {
    gen byte d_def_tg`c' = (d_def == 1 & tgcat == `c')
}
gen byte d_def_dm1 = (d_def == 1 & diabetes == 1)
gen byte d_def_dm0 = (d_def == 1 & diabetes == 0)
gen byte d_def_a65 = (d_def == 1 & ridageyr >= 65)
gen byte d_def_a64 = (d_def == 1 & ridageyr < 65)

*--- "above goal" indicators (assigned goal, same for every equation)
gen byte fA = f_above
gen byte mA = above_mh
gen byte sA = above_s
gen byte nA = nh_miss

*=======================================================================*
* A. MAIN TABLE
*=======================================================================*
tempname T
postfile `T' str24 block str40 label n pF lF uF pE lE uE pM lM uM pS lS uS ///
    gap gap_lb gap_ub n_up n_dn rec rec_lb rec_ub n_rec res res_lb res_ub nr_up nr_dn ///
    str32 flagF str32 flagM using "`out'/rev5_main_table.dta", replace

rowstats `T' "A. Main" "All adults with a defined goal" d_def fA mA sA nA
rowstats `T' "A. Risk group" "Primary prevention, borderline-high" d_def_prim fA mA sA nA
rowstats `T' "A. Risk group" "  Borderline/intermediate (<100)" d_def_g2 fA mA sA nA
rowstats `T' "A. Risk group" "  High PREVENT risk (<70)" d_def_g3 fA mA sA nA
rowstats `T' "A. Risk group" "ASCVD" d_def_sec fA mA sA nA
rowstats `T' "A. Risk group" "  ASCVD, not very high risk (<70)" d_def_g4 fA mA sA nA
rowstats `T' "A. Risk group" "  ASCVD, very high risk (<55)" d_def_g5 fA mA sA nA
rowstats `T' "A. Treatment and TG" "Treated, TG >=150" d_def_t150 fA mA sA nA
rowstats `T' "A. Treatment and TG" "Treated, TG <150" d_def_tlow fA mA sA nA
rowstats `T' "A. Treatment and TG" "Not treated" d_def_untr fA mA sA nA
rowstats `T' "A. TG, mg/dL" "<100" d_def_tg1 fA mA sA nA
rowstats `T' "A. TG, mg/dL" "100-149" d_def_tg2 fA mA sA nA
rowstats `T' "A. TG, mg/dL" "150-199" d_def_tg3 fA mA sA nA
rowstats `T' "A. TG, mg/dL" "200-399" d_def_tg4 fA mA sA nA
rowstats `T' "A. Other" "Diabetes" d_def_dm1 fA mA sA nA
rowstats `T' "A. Other" "No diabetes" d_def_dm0 fA mA sA nA
rowstats `T' "A. Other" "Age 65-79" d_def_a65 fA mA sA nA
rowstats `T' "A. Other" "Age 30-64" d_def_a64 fA mA sA nA
* comparison with Abohashem (all adults incl. low-risk 160 threshold)
rowstats `T' "Z. JAMA comparison" "All adults 30-79 incl. low risk" d_bm_all fA mA sA nA

*=======================================================================*
* B. SENSITIVITY: very-high-risk criterion by Friedewald
*=======================================================================*
gen n_hrc_f = (ridageyr >= 65) + (smoker == 1) + (dm == 1) + (hf == 1) ///
    + (bpq020 == 1) + (egfr < 60) + (treated == 1 & lbdldl >= 100 & lbdldl < .)
gen vhr_f = (n_major >= 2 | (n_major == 1 & n_hrc_f >= 2)) if ascvd == 1
display _n "=== VHR status: Martin-Hopkins criterion (rows) vs Friedewald criterion (cols), ASCVD 30-79 ==="
tab vhr vhr_f if d_bm_all == 1 & ascvd == 1, missing
gen goalgrp_f = goalgrp
replace goalgrp_f = 4 if goalgrp < . & ascvd == 1 & vhr_f == 0
replace goalgrp_f = 5 if goalgrp < . & ascvd == 1 & vhr_f == 1
gen goal_f = goal
replace goal_f = 70 if goalgrp_f == 4
replace goal_f = 55 if goalgrp_f == 5
gen byte fA_v = (lbdldl  >= goal_f)      if d_bm_all == 1 & goal_f < .
gen byte mA_v = (lbdldlm >= goal_f)      if d_bm_all == 1 & goal_f < .
gen byte sA_v = (lbdldln >= goal_f)      if d_bm_all == 1 & goal_f < .
gen byte nA_v = (nonhdl  >= goal_f + 30) if d_bm_all == 1 & goal_f < .
gen byte d_def_v     = (d_bm_all == 1 & inrange(goalgrp_f, 2, 5))
gen byte d_def_v_sec = (d_def_v == 1 & inlist(goalgrp_f, 4, 5))
gen byte d_def_v_g5  = (d_def_v == 1 & goalgrp_f == 5)
gen byte d_def_v_t150 = (d_def_v == 1 & treated == 1 & tg150 == 1)
rowstats `T' "B. VHR by Friedewald" "All adults with a defined goal" d_def_v fA_v mA_v sA_v nA_v
rowstats `T' "B. VHR by Friedewald" "ASCVD" d_def_v_sec fA_v mA_v sA_v nA_v
rowstats `T' "B. VHR by Friedewald" "  ASCVD, very high risk (<55)" d_def_v_g5 fA_v mA_v sA_v nA_v
rowstats `T' "B. VHR by Friedewald" "Treated, TG >=150" d_def_v_t150 fA_v mA_v sA_v nA_v

*=======================================================================*
* D. SENSITIVITY: TG on the old assay (backward bridging, rev3 block E)
*=======================================================================*
gen byte fA_b = (f_b  >= goal) if d_bm_all == 1 & f_b  < .
gen byte mA_b = (mh_b >= goal) if d_bm_all == 1 & mh_b < .
gen byte sA_b = (s_b  >= goal) if d_bm_all == 1 & s_b  < .
gen byte d_def_bt150 = (d_def == 1 & treated == 1 & tg150_b == 1)
rowstats `T' "D. TG bridged to old assay" "All adults with a defined goal" d_def fA_b mA_b sA_b nA
rowstats `T' "D. TG bridged to old assay" "Treated, TG >=150 (bridged)" d_def_bt150 fA_b mA_b sA_b nA
postclose `T'

*=======================================================================*
* C. SENSITIVITY: margin rule
*   A reclassification counts only if the equation that puts the person
*   above goal does so by at least m mg/dL (m = 5, 10).
*=======================================================================*
tempname M
postfile `M' str14 domain margin n gap gap_lb gap_ub n_up n_dn ///
    using "`out'/rev5_margin.dta", replace
foreach m in 0 5 10 {
    gen byte up`m' = (lbdldl < goal & lbdldlm >= goal + `m') if d_bm_all == 1 & goal < .
    gen byte dn`m' = (lbdldlm < goal & lbdldl >= goal + `m') if d_bm_all == 1 & goal < .
    foreach d in d_def d_def_prim d_def_sec d_def_t150 {
        quietly count if `d' == 1
        local n = r(N)
        moverkg up`m' dn`m' `d' 0
        post `M' ("`d'") (`m') (`n') (r(est)) (r(lb)) (r(ub)) (r(x1)) (r(x2))
    }
}
postclose `M'

*=======================================================================*
* E. TABLE 1 with 95% CIs
*=======================================================================*
gen byte hisp  = inlist(race6, 1, 2) if race6 < .
gen byte nhw   = (race6 == 3) if race6 < .
gen byte nhb   = (race6 == 4) if race6 < .
gen byte nha   = (race6 == 5) if race6 < .
gen byte oth   = (race6 == 6) if race6 < .
gen byte ascvd_g = inlist(goalgrp, 4, 5) if goalgrp < .
gen byte g100 = (goal == 100) if goal < .
gen byte g70  = (goal == 70)  if goal < .
gen byte g55  = (goal == 55)  if goal < .
gen byte trt  = treated
gen byte tg150v = tg150

tempname E
postfile `E' str12 column str24 row str6 kind n est lb ub str32 flag ///
    using "`out'/rev5_table1.dta", replace
foreach c in d_def d_def_t150 d_def_tlow d_def_untr {
    quietly count if `c' == 1
    post `E' ("`c'") ("n (unweighted)") ("n") (r(N)) (r(N)) (.) (.) ("")
    foreach v in ridageyr lbxtc lbdhdd nonhdl lbdldl lbdldlm lbdldln bmxbmi {
        quietly svy, subpop(`c'): mean `v'
        local t = invttail(e(df_r), 0.025)
        post `E' ("`c'") ("`v'") ("mean") (e(N_sub)) (_b[`v']) (_b[`v'] - `t'*_se[`v']) ///
            (_b[`v'] + `t'*_se[`v']) ("")
    }
    _pctile lbxtlg [pweight = wtsaf2yr] if `c' == 1, percentiles(25 50 75)
    post `E' ("`c'") ("TG median (IQR)") ("median") (.) (r(r2)) (r(r1)) (r(r3)) ("")
    foreach v in female nhw nhb hisp nha oth diabetes ascvd_g trt tg150v g100 g70 g55 {
        kgprop `v', sub(`c')
        post `E' ("`c'") ("`v'") ("pct") (r(n)) (r(p)) (r(lb)) (r(ub)) ("`r(flag)'")
    }
}
postclose `E'

*=======================================================================*
* F. SAMPLE RECONCILIATION with Abohashem et al
*=======================================================================*
merge 1:1 seqn using "$DATA2123/nhanes_prevent.dta", keep(master match) nogen ///
    keepusing(sbp bptx statin wtmec2yr cvd_hx)
gen byte base = (ridageyr >= 30 & ridageyr <= 79 & wtsaf2yr > 0 & wtsaf2yr < . & lbdldlm < .)
tempname F
postfile `F' str60 step n using "`out'/rev5_flow.dta", replace
quietly count if base
post `F' ("JAMA-style base: 30-79, fasting weight >0, Martin-Hopkins LDL-C") (r(N))
quietly count if base & ridexprg == 1
post `F' ("  of which pregnant (we exclude)") (r(N))
gen byte b1 = base & ridexprg != 1
quietly count if b1 & (lbdldl >= . | lbdldln >= . | nonhdl >= .)
post `F' ("  of which missing Friedewald/Sampson/non-HDL-C") (r(N))
gen byte b2 = b1 & lbdldl < . & lbdldln < . & nonhdl < .
quietly count if b2
post `F' ("After those exclusions") (r(N))
quietly count if b2 & ascvd == 1
post `F' ("  ASCVD (CHD, angina, MI, stroke) - ours") (r(N))
quietly count if base & ascvd == 1
post `F' ("  ASCVD in JAMA-style base (before our exclusions)") (r(N))
quietly count if base & (mcq160c == 1 | mcq160d == 1 | mcq160e == 1 | mcq160f == 1)
post `F' ("  any ASCVD 'yes' in base (ignoring missing items)") (r(N))
quietly count if b2 & ascvd >= .
post `F' ("  ASCVD history missing/don't know (we exclude)") (r(N))
quietly count if b2 & ascvd == 0 & hf == 1
post `F' ("  heart failure without ASCVD (both exclude)") (r(N))
quietly count if b2 & ascvd == 0 & hf >= .
post `F' ("  no ASCVD, heart-failure item missing (we exclude)") (r(N))
quietly count if b2 & ascvd == 0 & hf == 0 & s_complete != 1
post `F' ("  no ASCVD/HF, missing a PREVENT input (both exclude)") (r(N))
foreach v in lbxtc lbdhdd sbp bptx statin dm smoker egfr {
    quietly count if b2 & ascvd == 0 & hf == 0 & s_complete != 1 & missing(`v')
    post `F' ("     missing `v'") (r(N))
}
quietly count if b2 & ascvd == 0 & hf == 0 & s_complete != 1 & !(wtmec2yr > 0 & wtmec2yr < .)
post `F' ("     MEC weight 0/missing") (r(N))
quietly count if d_bm_all == 1 & inrange(goalgrp, 1, 3)
post `F' ("Our primary prevention n (JAMA 2,031)") (r(N))
quietly count if d_bm_all == 1 & inrange(goalgrp, 4, 5)
post `F' ("Our ASCVD n (JAMA 282)") (r(N))
quietly count if d_bm_all == 1
post `F' ("Our total n (JAMA 2,313)") (r(N))
postclose `F'

*=======================================================================*
* EXPORT
*=======================================================================*
foreach f in rev5_main_table rev5_margin rev5_table1 rev5_flow {
    use "`out'/`f'.dta", clear
    list, noobs clean abbreviate(10)
    export delimited using "`out'/`f'.csv", replace
}
log close
display "DONE"
