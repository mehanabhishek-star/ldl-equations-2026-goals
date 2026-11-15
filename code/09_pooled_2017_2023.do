*=======================================================================*
* 09_pooled_2017_2023.do
* Pooled NHANES 2017-March 2020 + August 2021-August 2023.
*   (1) 2017-March 2020 TG placed on the current assay with the NCHS forward
*       equation, TG(new) = 12.45 + 1.022 x LBXTR(old), rounded to an
*       integer as NCHS reports TG; Friedewald and Sampson recomputed;
*   (2) Martin-Hopkins at the bridged TG looked up exactly from the NCHS
*       values (results/mh_vldl_lookup.csv, step 07); the participant's own
*       original factor is kept as ldlm_own for comparison;
*   (3) Satterthwaite effective-df intervals for the pooled estimates;
*   (4) participant-level recovery in the pooled sample;
*   (5) 2017-March 2020 on its own (previous) assay, as reported by NCHS.
* Weights: WTSAFPRP x 3.2/5.2 and WTSAF2YR x 2.0/5.2.
* Uses _kgprop.do, _moverkg.do, _prevent_programs.do, _satt.do.
* Inputs (.xpt, data/nhanes_2017_2020): P_DEMO P_TRIGLY P_TCHOL P_HDL P_BPXO
*   P_BPQ P_DIQ P_GHB P_BIOPRO P_SMQ P_MCQ P_BMX P_RXQ_RX
* Outputs: results/rev6_pooled_table.csv, rev6_pooled_rx.csv,
*          rev6_pooled_satt.csv, rev6_pooled_person.csv
*=======================================================================*
version 17
clear all
set more off
set varabbrev off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
local out  "$RESULTS"
local code "$CODE"
local P    "$DATA1720"
capture log close
log using "`out'/rev6_pooled.log", text replace
do "`code'/_kgprop.do"
do "`code'/_moverkg.do"
do "`code'/_prevent_programs.do"
do "`code'/_satt.do"

import delimited using "`out'/mh_vldl_lookup.csv", clear
rename tg tgi
tempfile lk
save `lk'

*=======================================================================*
* 1. Build 2017-March 2020
*=======================================================================*
tempfile tr tc hd bp bq dq gh bi sm mc bm rx
import sasxport5 "`P'/P_TRIGLY.xpt", clear
keep seqn lbxtr lbdldl lbdldlm lbdldln wtsafprp
save `tr'
import sasxport5 "`P'/P_TCHOL.xpt", clear
keep seqn lbxtc
save `tc'
import sasxport5 "`P'/P_HDL.xpt", clear
keep seqn lbdhdd
save `hd'
import sasxport5 "`P'/P_BPXO.xpt", clear
keep seqn bpxosy1 bpxosy2 bpxosy3
save `bp'
import sasxport5 "`P'/P_BPQ.xpt", clear
keep seqn bpq020 bpq040a bpq050a bpq080 bpq090d bpq100d
save `bq'
import sasxport5 "`P'/P_DIQ.xpt", clear
keep seqn diq010 diq050 diq070
save `dq'
import sasxport5 "`P'/P_GHB.xpt", clear
keep seqn lbxgh
save `gh'
import sasxport5 "`P'/P_BIOPRO.xpt", clear
keep seqn lbxscr
save `bi'
import sasxport5 "`P'/P_SMQ.xpt", clear
keep seqn smq020 smq040
save `sm'
import sasxport5 "`P'/P_MCQ.xpt", clear
keep seqn mcq160b mcq160c mcq160d mcq160e mcq160f
save `mc'
import sasxport5 "`P'/P_BMX.xpt", clear
keep seqn bmxbmi
save `bm'
import sasxport5 "`P'/P_RXQ_RX.xpt", clear
gen s = upper(rxddrug)
gen byte rx_statin = regexm(s, "ATORVASTATIN|ROSUVASTATIN|SIMVASTATIN|PRAVASTATIN|LOVASTATIN|FLUVASTATIN|PITAVASTATIN")
gen byte rx_ezet   = regexm(s, "EZETIMIBE")
gen byte rx_pcsk9  = regexm(s, "ALIROCUMAB|EVOLOCUMAB")
gen byte rx_hiint  = regexm(s, "ATORVASTATIN|ROSUVASTATIN")
collapse (max) rx_statin rx_ezet rx_pcsk9 rx_hiint, by(seqn)
save `rx'

import sasxport5 "`P'/P_DEMO.xpt", clear
keep seqn riagendr ridageyr ridreth3 ridexprg sdmvpsu sdmvstra wtmecprp
foreach t in tr tc hd bp bq dq gh bi sm mc bm rx {
    merge 1:1 seqn using ``t'', keep(master match) nogen
}

gen female = (riagendr == 2) if riagendr < .
gen nonhdl = lbxtc - lbdhdd if lbxtc < . & lbdhdd < .
gen treated = .
replace treated = 1 if bpq100d == 1
replace treated = 0 if bpq100d == 2
replace treated = 0 if bpq090d == 2
replace treated = 0 if bpq090d >= . & bpq080 == 2
gen diabetes = (diq010 == 1) if inlist(diq010, 1, 2, 3)
replace diabetes = 1 if lbxgh >= 6.5 & lbxgh < .
egen sbp = rowmean(bpxosy1 bpxosy2 bpxosy3)
gen bptx = .
replace bptx = 0 if bpq020 == 2
replace bptx = 1 if bpq020 == 1 & bpq050a == 1
replace bptx = 0 if bpq020 == 1 & (bpq040a == 2 | bpq050a == 2)
gen statin = treated
gen dm = (diq010 == 1) if inlist(diq010, 1, 2, 3)
replace dm = 1 if diq050 == 1 | diq070 == 1 | (lbxgh >= 6.5 & lbxgh < .)
gen smoker = 0 if inlist(smq020, 1, 2)
replace smoker = 1 if smq020 == 1 & inlist(smq040, 1, 2)
gen double kappa = cond(female == 1, 0.7, 0.9)
gen double alpha = cond(female == 1, -0.241, -0.302)
gen double egfr = 142 * (min(lbxscr/kappa, 1))^alpha * (max(lbxscr/kappa, 1))^-1.200 ///
    * 0.9938^ridageyr * cond(female == 1, 1.012, 1) if lbxscr < .
drop kappa alpha
gen cvd_hx = .
replace cvd_hx = 0 if mcq160b == 2 & mcq160c == 2 & mcq160d == 2 & mcq160e == 2 & mcq160f == 2
replace cvd_hx = 1 if mcq160b == 1 | mcq160c == 1 | mcq160d == 1 | mcq160e == 1 | mcq160f == 1
gen ascvd = .
replace ascvd = 0 if mcq160c == 2 & mcq160d == 2 & mcq160e == 2 & mcq160f == 2
replace ascvd = 1 if mcq160c == 1 | mcq160d == 1 | mcq160e == 1 | mcq160f == 1
gen hf = (mcq160b == 1) if inlist(mcq160b, 1, 2)

*--- forward bridging to the 2021-2023 assay: integer TG, exact Martin-Hopkins lookup
gen double tg = round(12.45 + 1.022 * lbxtr) if lbxtr < .
gen double mhfac = lbxtr / (nonhdl - lbdldlm) if lbdldlm < . & (nonhdl - lbdldlm) > 0
gen double ldlm_own = round(nonhdl - (12.45 + 1.022 * lbxtr)/mhfac) if mhfac < .
gen double ldlf = round(lbxtc - lbdhdd - tg/5)
gen double ldls = round(lbxtc/0.948 - lbdhdd/0.971 - (tg/8.56 + tg*nonhdl/2140 - tg^2/16100) - 9.44)
gen byte band = 0 if nonhdl <= 99
replace band = 1 if nonhdl >= 100 & nonhdl <= 129
replace band = 2 if nonhdl >= 130 & nonhdl <= 159
replace band = 3 if nonhdl >= 160 & nonhdl <= 189
replace band = 4 if nonhdl >= 190 & nonhdl <= 219
replace band = 5 if nonhdl >= 220 & nonhdl < .
gen tgi = tg
merge m:1 band tgi using `lk', keep(master match) nogen
gen double ldlm = nonhdl - v if v < . & lbdldlm < .
display "Exact lookup minus own-factor approximation, mg/dL (2017-2020, forward bridge):"
gen dmh_own = ldlm - ldlm_own
tab dmh_own
* native values (previous assay, as reported by NCHS, no bridging)
gen double tg_nat   = lbxtr
gen double ldlf_nat = lbdldl
gen double ldlm_nat = lbdldlm
gen double ldls_nat = lbdldln

gen byte analytic = (ridageyr >= 20 & ridexprg != 1 & lbdldl < . & lbdldlm < . & lbdldln < . ///
    & nonhdl < . & wtsafprp > 0 & wtsafprp < . & tg < 400 & ldlm < .)
display "2017-2020 analytic (20+): " _c
count if analytic
display "  excluded because bridged TG >= 400 (original TG < 400):"
count if ridageyr >= 20 & ridexprg != 1 & lbdldl < . & tg >= 400 & tg < .

*--- PREVENT (same rules as prevent_build.do)
gen s_complete = (ridageyr >= 30 & ridageyr <= 79 & ridexprg != 1 & cvd_hx == 0 & lbxtc < . ///
    & lbdhdd < . & sbp < . & bptx < . & statin < . & dm < . & smoker < . & egfr < . ///
    & wtmecprp > 0 & wtmecprp < .)
gen p_age    = ridageyr
gen p_tc     = min(max(lbxtc, 130), 320)        if lbxtc < .
gen p_hdl    = min(max(lbdhdd, 20), 100)        if lbdhdd < .
gen p_sbp    = min(max(sbp, 90), 180)           if sbp < .
gen p_egfr   = min(max(egfr, 15), 140)          if egfr < .
gen p_bmi    = min(max(bmxbmi, 18.5), 39.9)     if bmxbmi < .
gen p_uacr   = .
gen p_a1c    = min(max(lbxgh, 4.5), 15)         if lbxgh < .
gen p_bptx   = bptx
gen p_statin = statin
gen p_dm     = dm
gen p_smk    = smoker
prevent_calc
replace r10_ascvd = . if !s_complete
capture drop p_*
capture drop t_*
capture drop lp_*
gen prevent_cat = 1 if r10_ascvd < 3
replace prevent_cat = 2 if r10_ascvd >= 3 & r10_ascvd < 5
replace prevent_cat = 3 if r10_ascvd >= 5 & r10_ascvd < 10
replace prevent_cat = 4 if r10_ascvd >= 10 & r10_ascvd < .

gen byte cycle = 1
gen double wt0 = wtsafprp
keep seqn cycle ridageyr ridexprg female ridreth3 sdmvpsu sdmvstra wt0 analytic tg ldlf ldlm ldls ///
    nonhdl treated diabetes dm smoker egfr bpq020 ascvd hf mcq160e mcq160f s_complete ///
    prevent_cat rx_statin rx_ezet rx_pcsk9 rx_hiint ldlm_own exact tg_nat ldlf_nat ldlm_nat ldls_nat
tempfile pc
save `pc'

*=======================================================================*
* 2. 2021-2023 (values exactly as in the main analysis)
*=======================================================================*
use "`out'/rev3_analysis_data.dta", clear
gen byte cycle = 2
gen double wt0 = wtsaf2yr
rename (lbxtlg lbdldl lbdldlm lbdldln) (tg ldlf ldlm ldls)
gen goalgrp_rev3 = goalgrp
keep seqn cycle ridageyr ridexprg female ridreth3 sdmvpsu sdmvstra wt0 analytic tg ldlf ldlm ldls ///
    nonhdl treated diabetes dm smoker egfr bpq020 ascvd hf mcq160e mcq160f s_complete ///
    prevent_cat goalgrp_rev3
append using `pc'

*=======================================================================*
* 3. Goals (same rules as ldl_rev3_person_goals.do), both cycles
*=======================================================================*
gen n_major = (mcq160e == 1) + (mcq160f == 1)
gen n_hrc = (ridageyr >= 65) + (smoker == 1) + (dm == 1) + (hf == 1) ///
          + (bpq020 == 1) + (egfr < 60) + (treated == 1 & ldlm >= 100 & ldlm < .)
gen vhr = (n_major >= 2 | (n_major == 1 & n_hrc >= 2)) if ascvd == 1
gen goalgrp = .
replace goalgrp = 1 if ascvd == 0 & s_complete == 1 & prevent_cat == 1
replace goalgrp = 2 if ascvd == 0 & s_complete == 1 & inlist(prevent_cat, 2, 3)
replace goalgrp = 3 if ascvd == 0 & s_complete == 1 & prevent_cat == 4
replace goalgrp = 4 if ascvd == 1 & vhr == 0
replace goalgrp = 5 if ascvd == 1 & vhr == 1
assert goalgrp == goalgrp_rev3 if cycle == 2
gen goal = 160 if goalgrp == 1
replace goal = 100 if goalgrp == 2
replace goal = 70  if inlist(goalgrp, 3, 4)
replace goal = 55  if goalgrp == 5
gen s30 = (analytic == 1 & ridageyr >= 30 & ridageyr <= 79)
gen byte fA = (ldlf >= goal)        if s30 == 1 & goal < .
gen byte mA = (ldlm >= goal)        if s30 == 1 & goal < .
gen byte sA = (ldls >= goal)        if s30 == 1 & goal < .
gen byte nA = (nonhdl >= goal + 30) if s30 == 1 & goal < .
gen byte tg150 = (tg >= 150) if tg < .
gen tgcat = 1 if tg < 100
replace tgcat = 2 if tg >= 100 & tg < 150
replace tgcat = 3 if tg >= 150 & tg < 200
replace tgcat = 4 if tg >= 200 & tg < 400

*--- pooled design
quietly summarize sdmvstra if cycle == 1
local p_lo = r(min)
local p_hi = r(max)
quietly summarize sdmvstra if cycle == 2
assert `p_hi' < r(min) | `p_lo' > r(max)
gen double wt = wt0 * cond(cycle == 1, 3.2/5.2, 2.0/5.2)
svyset sdmvpsu [pweight = wt], strata(sdmvstra) vce(linearized) singleunit(centered)
svydescribe, single
quietly svy: mean nonhdl
display "POOLED DESIGN df = " e(df_r) "  strata = " e(N_strata) "  PSUs = " e(N_psu)

gen byte d_def   = (s30 == 1 & inrange(goalgrp, 2, 5))
gen byte d_prim  = (d_def == 1 & inlist(goalgrp, 2, 3))
gen byte d_sec   = (d_def == 1 & inlist(goalgrp, 4, 5))
gen byte d_t150  = (d_def == 1 & treated == 1 & tg150 == 1)
gen byte d_tlow  = (d_def == 1 & treated == 1 & tg150 == 0)
gen byte d_untr  = (d_def == 1 & treated == 0)
forvalues c = 1/4 {
    gen byte d_tg`c' = (d_def == 1 & tgcat == `c')
}
forvalues y = 1/2 {
    gen byte d_def_c`y'  = (d_def == 1 & cycle == `y')
    gen byte d_t150_c`y' = (d_t150 == 1 & cycle == `y')
    gen byte d_sec_c`y'  = (d_sec == 1 & cycle == `y')
}
display _n "=== Unweighted n with a defined goal, by cycle ==="
tab cycle if d_def == 1
tab cycle if d_t150 == 1
display "Exact-cell lookups, 2017-2020 adults with a defined goal (1 = cell observed in NHANES):"
tab exact if d_def_c1 == 1, missing

tempname T
postfile `T' str24 block str40 label n pF lF uF pE lE uE pM lM uM pS lS uS ///
    gap gap_lb gap_ub n_up n_dn rec rec_lb rec_ub n_rec res res_lb res_ub nr_up nr_dn ///
    str32 flagF str32 flagM using "`out'/rev6_pooled_table.dta", replace
rowstats `T' "P. Pooled 2017-2023" "All adults with a defined goal" d_def fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "Primary prevention, borderline-high" d_prim fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "ASCVD" d_sec fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "Treated, TG >=150" d_t150 fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "Treated, TG <150" d_tlow fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "Not treated" d_untr fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "TG <100" d_tg1 fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "TG 100-149" d_tg2 fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "TG 150-199" d_tg3 fA mA sA nA
rowstats `T' "P. Pooled 2017-2023" "TG 200-399" d_tg4 fA mA sA nA
rowstats `T' "Q. 2017-2020 alone" "All adults with a defined goal" d_def_c1 fA mA sA nA
rowstats `T' "Q. 2017-2020 alone" "ASCVD" d_sec_c1 fA mA sA nA
rowstats `T' "Q. 2017-2020 alone" "Treated, TG >=150" d_t150_c1 fA mA sA nA
rowstats `T' "R. 2021-2023 (check = main)" "All adults with a defined goal" d_def_c2 fA mA sA nA
rowstats `T' "R. 2021-2023 (check = main)" "Treated, TG >=150" d_t150_c2 fA mA sA nA
* own-factor approximation (revision 5) for comparison
gen byte mA_own = (cond(cycle == 1, ldlm_own, ldlm) >= goal) if s30 == 1 & goal < .
rowstats `T' "O. Pooled, own-factor MH" "All adults with a defined goal" d_def fA mA_own sA nA
rowstats `T' "O. Pooled, own-factor MH" "Treated, TG >=150" d_t150 fA mA_own sA nA
* 2017-March 2020 on its native assay, as reported by NCHS (no bridging)
gen byte fA_n = (ldlf_nat >= goal) if s30 == 1 & goal < . & cycle == 1
gen byte mA_n = (ldlm_nat >= goal) if s30 == 1 & goal < . & cycle == 1
gen byte sA_n = (ldls_nat >= goal) if s30 == 1 & goal < . & cycle == 1
gen byte d_t150_nat = (d_def_c1 == 1 & treated == 1 & tg_nat >= 150 & tg_nat < .)
rowstats `T' "N. 2017-2020 native assay" "All adults with a defined goal" d_def_c1 fA_n mA_n sA_n nA
rowstats `T' "N. 2017-2020 native assay" "Treated, TG >=150 (native)" d_t150_nat fA_n mA_n sA_n nA
postclose `T'

*=======================================================================*
* 4. Satterthwaite intervals
*=======================================================================*
gen dif_mh = mA - fA
gen byte eF = (fA == 1 | nA == 1) if fA < . & nA < .
gen resid = mA - eF
tempname Q
postfile `Q' str14 domain str8 comp n est se se_man df_design lb_t ub_t df_satt lb_s ub_s top3 ///
    using "`out'/rev6_pooled_satt.dta", replace
foreach d in d_def d_t150 d_sec d_def_c1 d_def_c2 {
    foreach y in dif_mh resid {
        quietly count if `d' == 1 & `y' < .
        local n = r(N)
        satt `y', sub(`d') wt(wt)
        post `Q' ("`d'") ("`y'") (`n') (r(est)) (r(se)) (r(se_man)) (r(df_design)) (r(lb_t)) (r(ub_t)) ///
            (r(df_satt)) (r(lb_s)) (r(ub_s)) (r(top3))
    }
}
postclose `Q'

*=======================================================================*
* 5. Person-level recovery, pooled
*=======================================================================*
gen byte up  = (mA == 1 & fA == 0) if mA < .
gen byte dn  = (mA == 0 & fA == 1) if mA < .
gen byte rup = (mA == 1 & eF == 0) if mA < .
gen byte rdn = (mA == 0 & eF == 1) if mA < .
gen byte gross_before = (up == 1 | dn == 1)   if mA < .
gen byte gross_after  = (rup == 1 | rdn == 1) if mA < .
gen byte caught = nA
gen byte add_mh = (mA == 0 & nA == 1) if mA < .
gen byte d_up     = (d_def == 1 & up == 1)
gen byte d_up_tg1 = (d_up == 1 & tg150 == 1)
gen byte d_up_tg0 = (d_up == 1 & tg150 == 0)
tempname A
postfile `A' str16 domain str14 outcome n x p lb ub neff df str32 flag ///
    using "`out'/rev6_pooled_person.dta", replace
foreach y in up dn gross_before gross_after rup rdn add_mh {
    kgprop `y', sub(d_def)
    post `A' ("d_def") ("`y'") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) (r(neff)) (r(df)) ("`r(flag)'")
}
foreach d in d_up d_up_tg1 d_up_tg0 {
    kgprop caught, sub(`d')
    post `A' ("`d'") ("caught") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) (r(neff)) (r(df)) ("`r(flag)'")
}
kgprop tg150, sub(d_up)
post `A' ("d_up") ("tg150") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) (r(neff)) (r(df)) ("`r(flag)'")
postclose `A'
gen gap_up = ldlm - goal if d_up == 1
_pctile gap_up [pweight = wt], percentiles(25 50 75)
display "Pooled: weighted median excess among the 'up' group: " r(r2) " (IQR " r(r1) "-" r(r3) ")"

*=======================================================================*
* 6. Drug names, 2017-2020 treated adults with a defined goal
*=======================================================================*
gen byte d_trt_c1    = (d_def_c1 == 1 & treated == 1)
gen byte d_trt150_c1 = (d_t150_c1 == 1)
gen byte rx_any = (rx_statin == 1 | rx_ezet == 1 | rx_pcsk9 == 1) if rx_statin < .
gen byte rx_none = 1 - rx_any
tempname R
postfile `R' str14 domain str10 drug n x p lb ub str32 flag using "`out'/rev6_pooled_rx.dta", replace
foreach d in d_trt_c1 d_trt150_c1 {
    foreach v in rx_statin rx_hiint rx_ezet rx_pcsk9 rx_none {
        kgprop `v', sub(`d')
        if "`r(flag)'" == "too few" post `R' ("`d'") ("`v'") (r(n)) (r(x)) (.) (.) (.) ("too few")
        else post `R' ("`d'") ("`v'") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) ("`r(flag)'")
    }
}
postclose `R'
save "`out'/rev6_pooled_data.dta", replace

foreach f in rev6_pooled_table rev6_pooled_rx rev6_pooled_satt rev6_pooled_person {
    use "`out'/`f'.dta", clear
    list, noobs clean abbreviate(12)
    export delimited using "`out'/`f'.csv", replace
}
log close
display "DONE"
