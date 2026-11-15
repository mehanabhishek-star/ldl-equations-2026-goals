*=======================================================================*
* 08_participant_level_and_assay.do   (NHANES August 2021-August 2023)
*   A. Participant-level recovery: of adults above goal by Martin-Hopkins
*      but at goal by Friedewald ("up"), the weighted % above the paired
*      non-HDL-C goal; the same within TG >=150; participants classified
*      differently from Martin-Hopkins before and after adding the
*      non-HDL-C goal; weighted % in each direction; weighted median excess
*      among the "up" group.
*   B. Sampson minus Friedewald paired difference (MOVER-KG).
*   C. TG converted to the previous (glycerol-blanked) assay:
*      TG(old) = round(-12.19 + 0.9785 x TG); Friedewald and Sampson
*      recomputed and rounded as NCHS reports them; Martin-Hopkins =
*      non-HDL-C - v(band, TG) from results/mh_vldl_lookup.csv (step 07).
* Outputs: results/rev6_person.csv, rev6_sf.csv, rev6_bridge_rows.csv
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
log using "`out'/rev6_main.log", text replace
do "`code'/_kgprop.do"
do "`code'/_moverkg.do"

import delimited using "`out'/mh_vldl_lookup.csv", clear
rename tg tgi
tempfile lk
save `lk'

use "`out'/rev3_analysis_data.dta", clear
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)

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
gen byte fA = f_above
gen byte mA = above_mh
gen byte sA = above_s
gen byte nA = nh_miss
gen byte eF = (fA == 1 | nA == 1) if fA < . & nA < .
gen byte up  = (mA == 1 & fA == 0) if mA < .
gen byte dn  = (mA == 0 & fA == 1) if mA < .
gen byte rup = (mA == 1 & eF == 0) if mA < .
gen byte rdn = (mA == 0 & eF == 1) if mA < .
gen byte gross_before = (up == 1 | dn == 1)   if mA < .
gen byte gross_after  = (rup == 1 | rdn == 1) if mA < .

*=======================================================================*
* A. Person-level recovery
*=======================================================================*
tempname A
postfile `A' str16 domain str14 outcome n x p lb ub neff df str32 flag ///
    using "`out'/rev6_person.dta", replace
capture program drop postk
program define postk
    args h dom y
    kgprop `y', sub(`dom')
    if "`r(flag)'" == "too few" post `h' ("`dom'") ("`y'") (r(n)) (r(x)) (.) (.) (.) (.) (.) ("too few")
    else post `h' ("`dom'") ("`y'") (r(n)) (r(x)) (r(p)) (r(lb)) (r(ub)) (r(neff)) (r(df)) ("`r(flag)'")
end
foreach y in up dn gross_before gross_after rup rdn {
    postk `A' d_def `y'
}
gen byte d_up      = (d_def == 1 & up == 1)
gen byte d_up_tg1  = (d_up == 1 & tg150 == 1)
gen byte d_up_tg0  = (d_up == 1 & tg150 == 0)
gen byte caught    = nA
postk `A' d_up caught
postk `A' d_up_tg1 caught
postk `A' d_up_tg0 caught
postk `A' d_up tg150
gen byte trt150 = (treated == 1 & tg150 == 1) if treated < . & tg150 < .
postk `A' d_up trt150
postclose `A'

gen gap_up = lbdldlm - goal if d_up == 1
_pctile gap_up [pweight = wtsaf2yr], percentiles(25 50 75)
display "Weighted median excess above goal among the 'up' group with a defined goal: " r(r2) " (IQR " r(r1) "-" r(r3) ")"
quietly count if d_up == 1 & gap_up < .
display "   n = " r(N)

*=======================================================================*
* B. Sampson - Friedewald (MOVER-KG)
*=======================================================================*
gen byte sup = (sA == 1 & fA == 0) if sA < .
gen byte sdn = (sA == 0 & fA == 1) if sA < .
tempname S
postfile `S' str14 domain n est lb ub n_up n_dn using "`out'/rev6_sf.dta", replace
foreach d in d_def d_def_prim d_def_g2 d_def_g3 d_def_sec d_def_g4 d_def_g5 d_def_t150 d_def_tlow d_def_untr ///
             d_def_tg1 d_def_tg2 d_def_tg3 d_def_tg4 {
    quietly count if `d' == 1
    local n = r(N)
    moverkg sup sdn `d' 0
    post `S' ("`d'") (`n') (r(est)) (r(lb)) (r(ub)) (r(x1)) (r(x2))
}
postclose `S'

*=======================================================================*
* C. Previous TG assay, exact Martin-Hopkins lookup
*=======================================================================*
gen tgi = round(-12.19 + 0.9785 * lbxtlg) if d_bm_all == 1
gen byte band = 0 if nonhdl <= 99
replace band = 1 if nonhdl >= 100 & nonhdl <= 129
replace band = 2 if nonhdl >= 130 & nonhdl <= 159
replace band = 3 if nonhdl >= 160 & nonhdl <= 189
replace band = 4 if nonhdl >= 190 & nonhdl <= 219
replace band = 5 if nonhdl >= 220 & nonhdl < .
merge m:1 band tgi using `lk', keep(master match) nogen
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)
* check: the lookup reproduces the reported Martin-Hopkins value on the ORIGINAL TG
preserve
    keep if d_bm_all == 1
    keep seqn band lbxtlg nonhdl lbdldlm
    rename lbxtlg tgi
    merge m:1 band tgi using `lk', keep(master match) nogen
    gen chk = nonhdl - v
    assert chk == lbdldlm
    display "Lookup reproduces LBDLDLM exactly for all " _N " participants on the original TG"
restore
gen double f_x  = round(lbxtc - lbdhdd - tgi/5) if d_bm_all == 1
gen double s_x  = round(lbxtc/0.948 - lbdhdd/0.971 - (tgi/8.56 + tgi*nonhdl/2140 - tgi^2/16100) - 9.44) if d_bm_all == 1
gen double mh_x = nonhdl - v if d_bm_all == 1 & v < .
display "Exact-cell lookups among adults with a defined goal (1 = cell observed in NHANES):"
tab exact if d_def == 1, missing
gen byte fA_x = (f_x  >= goal) if d_bm_all == 1 & mh_x < .
gen byte mA_x = (mh_x >= goal) if d_bm_all == 1 & mh_x < .
gen byte sA_x = (s_x  >= goal) if d_bm_all == 1 & mh_x < .
gen byte nA_x = nA if mh_x < .
gen byte d_x      = (d_def == 1 & mh_x < .)
gen byte d_x_t150 = (d_x == 1 & treated == 1 & tgi >= 150)
gen byte d_x_sec  = (d_x == 1 & inlist(goalgrp, 4, 5))
tempname T
postfile `T' str24 block str40 label n pF lF uF pE lE uE pM lM uM pS lS uS ///
    gap gap_lb gap_ub n_up n_dn rec rec_lb rec_ub n_rec res res_lb res_ub nr_up nr_dn ///
    str32 flagF str32 flagM using "`out'/rev6_bridge_rows.dta", replace
rowstats `T' "X. Previous TG assay, exact" "All adults with a defined goal" d_x fA_x mA_x sA_x nA_x
rowstats `T' "X. Previous TG assay, exact" "ASCVD" d_x_sec fA_x mA_x sA_x nA_x
rowstats `T' "X. Previous TG assay, exact" "Treated, TG >=150 (previous assay)" d_x_t150 fA_x mA_x sA_x nA_x
postclose `T'
display "Mean change in Martin-Hopkins LDL-C on the previous assay (exact - reported), mg/dL:"
gen dmh = mh_x - lbdldlm if d_def == 1
summarize dmh
gen dmh_own = mh_x - round(mh_b) if d_def == 1
display "Exact lookup minus own-factor approximation (rev5), mg/dL:"
tab dmh_own

foreach f in rev6_person rev6_sf rev6_bridge_rows {
    use "`out'/`f'.dta", clear
    list, noobs clean abbreviate(12)
    export delimited using "`out'/`f'.csv", replace
}
log close
display "DONE"
