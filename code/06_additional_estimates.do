*=======================================================================*
* 06_additional_estimates.do
*   - very-high-risk LDL-C criterion assigned by Sampson
*   - what the paired non-HDL-C goal adds when LDL-C comes from
*     Martin-Hopkins vs Friedewald
*   - exploratory survey-weighted linear model of the paired difference
*     (TG category, goal strictness, treatment)
* Uses _kgprop.do and _moverkg.do.
* Outputs: results/rev5_nonhdl_added.csv, rev5_vhr_sampson.csv,
*          rev5_model.csv
*=======================================================================*
version 17
clear all
set more off
set varabbrev off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
local out  "$RESULTS"
local code "$CODE"
capture log close
log using "`out'/rev5_nonhdl_added.log", text replace
do "`code'/_kgprop.do"
do "`code'/_moverkg.do"
use "`out'/rev3_analysis_data.dta", clear
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)

gen byte d_def      = (d_bm_all == 1 & formal_goal == 1)
gen byte d_def_sec  = (d_def == 1 & inlist(goalgrp, 4, 5))
gen byte d_def_t150 = (d_def == 1 & treated == 1 & tg150 == 1)
gen byte fA = f_above
gen byte mA = above_mh
gen byte sA = above_s
gen byte nA = nh_miss

*--- (51) what the non-HDL-C goal adds to each equation
gen byte add_mh = (mA == 0 & nA == 1) if mA < . & nA < .
gen byte add_f  = (fA == 0 & nA == 1) if fA < . & nA < .
tempname E
postfile `E' str14 domain str8 outcome n x p lb ub str32 flag using "`out'/rev5_nonhdl_added.dta", replace
foreach d in d_def d_def_sec d_def_t150 {
    foreach y in add_mh add_f {
        kgprop `y', sub(`d')
        post `E' ("`d'") ("`y'") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) ("`r(flag)'")
    }
}
postclose `E'

*--- (29) very-high-risk criterion by Sampson
gen n_hrc_s = (ridageyr >= 65) + (smoker == 1) + (dm == 1) + (hf == 1) ///
    + (bpq020 == 1) + (egfr < 60) + (treated == 1 & lbdldln >= 100 & lbdldln < .)
gen vhr_s = (n_major >= 2 | (n_major == 1 & n_hrc_s >= 2)) if ascvd == 1
display _n "=== VHR: Martin-Hopkins criterion (rows) vs Sampson criterion (cols), ASCVD ==="
tab vhr vhr_s if d_bm_all == 1 & ascvd == 1, missing
gen goalgrp_s = goalgrp
replace goalgrp_s = 4 if goalgrp < . & ascvd == 1 & vhr_s == 0
replace goalgrp_s = 5 if goalgrp < . & ascvd == 1 & vhr_s == 1
gen goal_s = goal
replace goal_s = 70 if goalgrp_s == 4
replace goal_s = 55 if goalgrp_s == 5
gen byte fA_s = (lbdldl  >= goal_s)      if d_bm_all == 1 & goal_s < .
gen byte mA_s = (lbdldlm >= goal_s)      if d_bm_all == 1 & goal_s < .
gen byte sA_s = (lbdldln >= goal_s)      if d_bm_all == 1 & goal_s < .
gen byte nA_s = (nonhdl  >= goal_s + 30) if d_bm_all == 1 & goal_s < .
gen byte d_v     = (d_bm_all == 1 & inrange(goalgrp_s, 2, 5))
gen byte d_v_sec = (d_v == 1 & inlist(goalgrp_s, 4, 5))
gen byte d_v_g5  = (d_v == 1 & goalgrp_s == 5)
quietly count if d_bm_all == 1 & ascvd == 1 & vhr != vhr_s
display "ASCVD participants whose very-high-risk status changes with the Sampson criterion: " r(N)
local nchg = r(N)
tempname T
postfile `T' str24 block str40 label n pF lF uF pE lE uE pM lM uM pS lS uS ///
    gap gap_lb gap_ub n_up n_dn rec rec_lb rec_ub n_rec res res_lb res_ub nr_up nr_dn ///
    str32 flagF str32 flagM using "`out'/rev5_vhr_sampson.dta", replace
rowstats `T' "S. VHR by Sampson" "All adults with a defined goal" d_v fA_s mA_s sA_s nA_s
rowstats `T' "S. VHR by Sampson" "ASCVD" d_v_sec fA_s mA_s sA_s nA_s
rowstats `T' "S. VHR by Sampson" "  ASCVD, very high risk (<55)" d_v_g5 fA_s mA_s sA_s nA_s
postclose `T'

*--- (36) one exploratory model, as requested: TG category, strict goal, treatment
gen strict = inlist(goalgrp, 3, 4, 5) if goalgrp < .
gen y = 100 * dif_mh
svy, subpop(d_def): regress y i.tgcat i.strict i.treated
tempname A
postfile `A' str10 term coef lb ub p n_sub df using "`out'/rev5_model.dta", replace
local t = invttail(e(df_r), 0.025)
foreach x in 2.tgcat 3.tgcat 4.tgcat 1.strict 1.treated {
    quietly test `x'
    post `A' ("`x'") (_b[`x']) (_b[`x'] - `t'*_se[`x']) (_b[`x'] + `t'*_se[`x']) (r(p)) (e(N_sub)) (e(df_r))
}
postclose `A'

foreach f in rev5_nonhdl_added rev5_vhr_sampson rev5_model {
    use "`out'/`f'.dta", clear
    list, noobs clean abbreviate(10)
    export delimited using "`out'/`f'.csv", replace
}
display "VHR status changes with Sampson criterion: `nchg'"
log close
