*=======================================================================*
* 03_goals_and_analysis_data.do
* Person-specific 2026 ACC/AHA LDL-C goals, following Abohashem et al,
* JAMA 2026;336:605-8:
*   Primary prevention (no ASCVD, no heart failure), PREVENT-ASCVD 10-year:
*       low <3%                  -> 160 mg/dL treatment threshold (no goal)
*       borderline/intermediate  -> <100
*       high >=10%               -> <70
*   ASCVD (CHD, angina, MI, stroke; MCQ160C/D/E/F), not very high risk -> <70
*   ASCVD at very high risk (VHR)                                      -> <55
*       VHR = MI and stroke, or (MI or stroke) + >=2 of: age >=65, current
*             smoking, diabetes, heart failure, hypertension, eGFR <60,
*             Martin-Hopkins LDL-C >=100 on lipid-lowering therapy
*   Paired non-HDL-C goal = LDL-C goal + 30 mg/dL.
* Goals are assigned once per person and applied to every equation.
* Inputs: ldl_reclass_data.dta (step 01), nhanes_prevent.dta (step 02)
* Output: results/rev3_analysis_data.dta (analysis dataset for steps 04-09)
*=======================================================================*
version 17
clear all
set more off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
cd "$DATA2123"
local out "$RESULTS"
capture log close
log using "`out'/rev3_person_goals.log", text replace


*=======================================================================*
* PART 1. KORN-GRAUBARD PROPORTION WITH NCHS RELIABILITY FLAGS
*   kgprop y, sub(s)   y = 0/1 outcome, s = 0/1 domain indicator
*   Returns r(n) r(x) r(p) r(lb) r(ub) r(neff) r(df) r(mil) r(flag)
*=======================================================================*
capture program drop kgprop
program define kgprop, rclass
    syntax varname, SUB(varname)
    tempvar dom tags tagp
    quietly gen byte `dom' = (`sub' == 1 & `varlist' < .)
    quietly count if `dom'
    local n = r(N)
    quietly count if `dom' & `varlist' == 1
    local x = r(N)
    if `n' < 2 {
        return scalar n = `n'
        return scalar x = `x'
        return local flag "too few"
        exit
    }
    * domain degrees of freedom = PSUs - strata containing the domain
    quietly egen `tags' = tag(sdmvstra)         if `dom'
    quietly egen `tagp' = tag(sdmvstra sdmvpsu) if `dom'
    quietly count if `tags' == 1
    local nstr = r(N)
    quietly count if `tagp' == 1
    local npsu = r(N)
    local df = `npsu' - `nstr'

    quietly svy, subpop(`dom'): mean `varlist'
    local p  = _b[`varlist']
    local se = _se[`varlist']
    quietly svy, subpop(`dom'): total `varlist'
    local mil = _b[`varlist'] / 1e6

    * effective sample size (capped at n; = n when p is 0 or 1)
    if `p' > 0 & `p' < 1 & `se' > 0 {
        local neff = min(`p' * (1 - `p') / `se'^2, `n')
    }
    else local neff = `n'
    * degrees-of-freedom adjustment
    if `df' >= 1 {
        local ndf = `neff' * (invttail(`n' - 1, 0.025) / invttail(`df', 0.025))^2
    }
    else local ndf = `neff'
    local xx = `p' * `ndf'
    if `xx' <= 0 {
        local lb = 0
        local ub = invibeta(1, `ndf', 0.975)
    }
    else if `xx' >= `ndf' {
        local lb = invibeta(`ndf', 1, 0.025)
        local ub = 1
    }
    else {
        local lb = invibeta(`xx', `ndf' - `xx' + 1, 0.025)
        local ub = invibeta(`xx' + 1, `ndf' - `xx', 0.975)
    }

    * NCHS reliability flags
    local w  = 100 * (`ub' - `lb')
    local flag "ok"
    if `p' == 0 | `p' == 1 local flag "review: 0/100%"
    if `df' < 8 local flag "review: df<8"
    if `w' > 5 & `w' < 30 & `p' > 0 {
        if `w' / (100 * `p') > 1.30 local flag "SUPPRESS: rel CI width>130%"
    }
    if `w' >= 30   local flag "SUPPRESS: CI width>=30"
    if `neff' < 30 local flag "SUPPRESS: n_eff<30"

    return scalar n    = `n'
    return scalar x    = `x'
    return scalar p    = 100 * `p'
    return scalar lb   = 100 * `lb'
    return scalar ub   = 100 * `ub'
    return scalar neff = `neff'
    return scalar df   = `df'
    return scalar mil  = `mil'
    return local  flag "`flag'"
end


*=======================================================================*
* PART 2. DATA: merge PREVENT categories and ASCVD history
*=======================================================================*
use ldl_reclass_data, clear
assert _N == 11933
merge 1:1 seqn using nhanes_prevent, keep(master match) nogen ///
    keepusing(mcq160b mcq160c mcq160d mcq160e mcq160f bpq020 smoker dm ///
              egfr prevent_cat s_complete r10_ascvd)
assert _N == 11933

svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)

* ASCVD = CHD, angina, MI or stroke (JAMA definition; HF alone is not ASCVD)
gen ascvd = .
replace ascvd = 0 if mcq160c == 2 & mcq160d == 2 & mcq160e == 2 & mcq160f == 2
replace ascvd = 1 if mcq160c == 1 | mcq160d == 1 | mcq160e == 1 | mcq160f == 1
gen hf = (mcq160b == 1) if inlist(mcq160b, 1, 2)

* Very high risk (JAMA eTable 2; LDL criterion uses Martin-Hopkins)
gen n_major = (mcq160e == 1) + (mcq160f == 1)
gen n_hrc = (ridageyr >= 65) + (smoker == 1) + (dm == 1) + (hf == 1) ///
          + (bpq020 == 1) + (egfr < 60) + (treated == 1 & lbdldlm >= 100 & lbdldlm < .)
gen vhr = (n_major >= 2 | (n_major == 1 & n_hrc >= 2)) if ascvd == 1

* Goal group and personal LDL-C goal
gen goalgrp = .
replace goalgrp = 1 if ascvd == 0 & s_complete == 1 & prevent_cat == 1
replace goalgrp = 2 if ascvd == 0 & s_complete == 1 & inlist(prevent_cat, 2, 3)
replace goalgrp = 3 if ascvd == 0 & s_complete == 1 & prevent_cat == 4
replace goalgrp = 4 if ascvd == 1 & vhr == 0
replace goalgrp = 5 if ascvd == 1 & vhr == 1
label define ggl 1 "Low <3% (<160, threshold)" 2 "Borderline/intermediate (<100)" ///
    3 "High >=10% (<70)" 4 "ASCVD, not VHR (<70)" 5 "ASCVD, VHR (<55)", replace
label values goalgrp ggl

gen goal = 160 if goalgrp == 1
replace goal = 100 if goalgrp == 2
replace goal = 70  if inlist(goalgrp, 3, 4)
replace goal = 55  if goalgrp == 5
gen goal_nh = goal + 30
gen formal_goal = inrange(goalgrp, 2, 5) if goalgrp < .

* Sample: fasting analytic sample (TG <400, all 3 LDL + non-HDL), 30-79
gen s30 = (analytic == 1 & ridageyr >= 30 & ridageyr <= 79)

display _n "=== SAMPLE FLOW (rev3) ==="
count if analytic == 1
count if s30 == 1
display "  ages 30-79 with an assigned goal group:"
count if s30 == 1 & goalgrp < .
display "  no goal group (missing PREVENT inputs / HF without ASCVD / missing history):"
count if s30 == 1 & goalgrp == .
tab goalgrp if s30 == 1, missing
display "  treated, TG >=150, by goal group:"
tab goalgrp if s30 == 1 & treated == 1 & tg150 == 1, missing

* Status by each equation against the personal goal
gen f_at    = (lbdldl  < goal)    if s30 == 1 & goal < .
gen mh_miss = (lbdldlm >= goal)   if s30 == 1 & goal < .
gen s_miss  = (lbdldln >= goal)   if s30 == 1 & goal < .
gen nh_miss = (nonhdl  >= goal_nh) if s30 == 1 & goal < .
gen any_miss  = (mh_miss == 1 | s_miss == 1 | nh_miss == 1) if mh_miss < .
gen conc_miss = (mh_miss == 1 & nh_miss == 1) if mh_miss < .
gen mh_m5   = (lbdldlm >= goal + 5)  if mh_miss < .
gen mh_m10  = (lbdldlm >= goal + 10) if mh_miss < .
gen nh_m10  = (nonhdl  >= goal_nh + 10) if mh_miss < .
gen f_above  = 1 - f_at

* Domains
* PRIMARY: treated, TG >=150, formal goal, at personal goal by Friedewald
gen d_primary = (s30 == 1 & treated == 1 & tg150 == 1 & formal_goal == 1 & f_at == 1)
* Comparison: same, TG <150
gen d_lowtg   = (s30 == 1 & treated == 1 & tg150 == 0 & formal_goal == 1 & f_at == 1)
* All treated at personal goal by Friedewald (any TG)
gen d_alltrt  = (s30 == 1 & treated == 1 & formal_goal == 1 & f_at == 1)
* Sensitivity: primary + low-risk with the <160 threshold
gen d_prim_low = (s30 == 1 & treated == 1 & tg150 == 1 & goal < . & f_at == 1)
* Untreated, TG >=150, formal goal, at goal by Friedewald
gen d_untr    = (s30 == 1 & treated == 0 & tg150 == 1 & formal_goal == 1 & f_at == 1)

display _n "Unweighted n: primary " _c
count if d_primary


*=======================================================================*
* PART 3. ESTIMATES -> postfile
*=======================================================================*
tempname H
postfile `H' str40 block str44 domain str28 outcome n x p lb ub neff df mil ///
    str32 flag using "`out'/rev3_person_goals.dta", replace

capture program drop postkg
program define postkg
    args handle block domain outcome y d
    kgprop `y', sub(`d')
    if "`r(flag)'" == "too few" {
        post `handle' ("`block'") ("`domain'") ("`outcome'") (r(n)) (r(x)) ///
            (.) (.) (.) (.) (.) (.) ("too few")
    }
    else {
        post `handle' ("`block'") ("`domain'") ("`outcome'") (r(n)) (r(x)) ///
            (r(p)) (r(lb)) (r(ub)) (r(neff)) (r(df)) (r(mil)) ("`r(flag)'")
    }
end

*--- A. PRIMARY and comparison domains
foreach d in d_primary d_lowtg d_alltrt d_prim_low d_untr {
    foreach y in mh_miss s_miss nh_miss conc_miss any_miss mh_m5 mh_m10 nh_m10 {
        postkg `H' "A. Discordance at personal goal" "`d'" "`y'" `y' `d'
    }
}

*--- B. Primary domain split by goal group (supplement; expect suppression)
forvalues k = 2/5 {
    gen d_gg`k' = (d_primary == 1 & goalgrp == `k')
    local lab : label ggl `k'
    foreach y in mh_miss s_miss nh_miss {
        postkg `H' "B. Primary, by goal group" "`lab'" "`y'" `y' d_gg`k'
    }
}

*--- C. Benchmark check: replicate JAMA (Martin-Hopkins, all adults 30-79)
*       then show the same quantity with Friedewald and Sampson
*       JAMA: primary prevention 32.2% above goal; ASCVD 79.6%
gen above_mh = mh_miss
gen above_s  = s_miss
gen d_bm_prim = (s30 == 1 & inrange(goalgrp, 1, 3))
gen d_bm_sec  = (s30 == 1 & inrange(goalgrp, 4, 5))
gen d_bm_all  = (s30 == 1 & goalgrp < .)
foreach d in d_bm_prim d_bm_sec d_bm_all {
    foreach y in above_mh f_above above_s {
        postkg `H' "C. Above personal goal (JAMA check)" "`d'" "`y'" `y' `d'
    }
}
* People above goal by MH who would look at goal with Friedewald
gen hidden = (mh_miss == 1 & f_at == 1) if mh_miss < .
foreach d in d_bm_all d_bm_prim d_bm_sec {
    postkg `H' "C2. Above by MH but at goal by F" "`d'" "hidden (all)" hidden `d'
    gen `d'_abv = (`d' == 1 & mh_miss == 1)
    postkg `H' "C2. Above by MH but at goal by F" "`d' & above by MH" "hidden (share)" hidden `d'_abv
}
* Same, restricted to treated adults with TG >=150
gen d_bm_trt150 = (d_bm_all == 1 & treated == 1 & tg150 == 1)
postkg `H' "C2. Above by MH but at goal by F" "treated TG>=150" "hidden (all)" hidden d_bm_trt150
foreach y in above_mh f_above above_s {
    postkg `H' "C. Above personal goal (JAMA check)" "d_bm_trt150" "`y'" `y' d_bm_trt150
}

*--- D. Context: who is at their own goal at all? (treated, formal goal)
gen at_mh  = 1 - mh_miss
gen at_s   = 1 - s_miss
gen at_nh  = 1 - nh_miss
gen at_all = (f_at == 1 & mh_miss == 0 & s_miss == 0 & nh_miss == 0) if mh_miss < .
gen d_trt150_formal = (s30 == 1 & treated == 1 & tg150 == 1 & formal_goal == 1)
gen d_trtlow_formal = (s30 == 1 & treated == 1 & tg150 == 0 & formal_goal == 1)
foreach d in d_trt150_formal d_trtlow_formal {
    foreach y in f_at at_mh at_s at_nh at_all {
        postkg `H' "D. At own goal (treated, formal goal)" "`d'" "`y'" `y' `d'
    }
}

*--- E. Triglyceride assay sensitivity: NCHS backward bridging equation
*   TRIGLY_L documentation: LBXTR(old) = -12.19 + 0.9785 * LBXTLG(new)
*   Friedewald and Sampson recomputed with bridged TG; Martin-Hopkins
*   recomputed with each person's own NHANES adjustable factor
*   (factor = TG / (non-HDL - MH LDL-C)), an approximation because the
*   360-cell table is not public in NHANES files.
gen double tg_b = -12.19 + 0.9785 * lbxtlg if analytic == 1
* check the Sampson formula reproduces NHANES LBDLDLN before using it
gen double s_chk = lbxtc/0.948 - lbdhdd/0.971 - (lbxtlg/8.56 + lbxtlg*nonhdl/2140 ///
    - lbxtlg^2/16100) - 9.44 if analytic == 1
gen double s_err = abs(s_chk - lbdldln)
summarize s_err if analytic == 1
assert s_err < 1.01 if analytic == 1
gen double mhfac = lbxtlg / (nonhdl - lbdldlm) if analytic == 1 & nonhdl - lbdldlm > 0
gen double f_b  = lbxtc - lbdhdd - tg_b/5 if analytic == 1
gen double mh_b = nonhdl - tg_b/mhfac   if mhfac < .
gen double s_b  = lbxtc/0.948 - lbdhdd/0.971 - (tg_b/8.56 + tg_b*nonhdl/2140 ///
    - tg_b^2/16100) - 9.44 if analytic == 1
gen tg150_b    = (tg_b >= 150) if tg_b < .
gen f_at_b     = (f_b  < goal)  if s30 == 1 & goal < . & f_b < .
gen mh_miss_b  = (mh_b >= goal) if s30 == 1 & goal < . & mh_b < .
gen s_miss_b   = (s_b  >= goal) if s30 == 1 & goal < . & s_b < .
gen hidden_b   = (mh_miss_b == 1 & f_at_b == 1) if mh_miss_b < . & f_at_b < .
gen d_primary_b = (s30 == 1 & treated == 1 & tg150_b == 1 & formal_goal == 1 & f_at_b == 1)
gen d_lowtg_b   = (s30 == 1 & treated == 1 & tg150_b == 0 & formal_goal == 1 & f_at_b == 1)
foreach y in mh_miss_b s_miss_b {
    postkg `H' "E. TG bridged (backward equation)" "d_primary_b" "`y'" `y' d_primary_b
    postkg `H' "E. TG bridged (backward equation)" "d_lowtg_b" "`y'" `y' d_lowtg_b
}
gen d_bm_all_abv_b = (d_bm_all == 1 & mh_miss_b == 1)
postkg `H' "E. TG bridged (backward equation)" "d_bm_all" "hidden_b (all)" hidden_b d_bm_all
postkg `H' "E. TG bridged (backward equation)" "d_bm_all & above by MH" "hidden_b (share)" hidden_b d_bm_all_abv_b

postclose `H'

*--- F. Paired differences in % above own goal: MH (or Sampson) minus Friedewald
gen dif_mh = above_mh - f_above
gen dif_s  = above_s  - f_above
tempname D
postfile `D' str24 domain str8 comp n pct lb ub mil mil_lb mil_ub ///
    using "`out'/rev3_paired_diff.dta", replace
foreach d in d_bm_all d_bm_prim d_bm_sec d_bm_trt150 {
    foreach v in dif_mh dif_s hidden {
        quietly count if `d' == 1 & `v' < .
        local n = r(N)
        quietly svy, subpop(`d'): mean `v'
        local t = invttail(e(df_r), 0.025)
        local b = _b[`v']
        local s = _se[`v']
        quietly svy, subpop(`d'): total `v'
        local tb = _b[`v']
        local ts = _se[`v']
        post `D' ("`d'") ("`v'") (`n') (100*`b') (100*(`b' - `t'*`s')) (100*(`b' + `t'*`s')) ///
            (`tb'/1e6) ((`tb' - `t'*`ts')/1e6) ((`tb' + `t'*`ts')/1e6)
    }
}
postclose `D'


*=======================================================================*
* PART 4. TESTS AND MAGNITUDE
*=======================================================================*
display _n "=== Rao-Scott: TG >=150 vs <150 among treated at personal goal (Friedewald) ==="
foreach y in mh_miss s_miss nh_miss {
    svy, subpop(d_alltrt): tab tg150 `y', row percent
    display "   `y': Rao-Scott p = " %6.4f e(p_Pear)
}

display _n "=== Treated vs untreated (TG >=150, at personal goal by Friedewald) ==="
gen d_tg150_atgoal = (d_primary == 1 | d_untr == 1)
svy, subpop(d_tg150_atgoal): tab treated mh_miss, row percent
display "   Rao-Scott p = " %6.4f e(p_Pear)

display _n "=== How far above goal (primary domain, discordant by MH) ==="
gen gap_mh = lbdldlm - goal if d_primary == 1 & mh_miss == 1
gen gap_nh = nonhdl - goal_nh if d_primary == 1 & nh_miss == 1
foreach v in gap_mh gap_nh {
    quietly count if `v' < .
    if r(N) > 0 {
        _pctile `v' [pweight = wtsaf2yr], percentiles(25 50 75)
        display "   `v': median " %4.1f r(r2) " (IQR " %4.1f r(r1) " to " %4.1f r(r3) ") mg/dL"
    }
}
gen d_mh_f = lbdldlm - lbdldl if d_primary == 1
_pctile d_mh_f [pweight = wtsaf2yr], percentiles(25 50 75)
display "   MH - F in primary domain: median " %4.1f r(r2) " (IQR " %4.1f r(r1) " to " %4.1f r(r3) ")"

display _n "=== Base population of the primary domain ==="
gen one = 1
svy, subpop(d_primary): total one
display "   weighted millions at personal goal by Friedewald: " %5.2f _b[one]/1e6
svy, subpop(d_primary): tab goalgrp, percent


*=======================================================================*
* PART 5. SAVE AND EXPORT
*=======================================================================*
save "`out'/rev3_analysis_data.dta", replace

use "`out'/rev3_person_goals.dta", clear
format p lb ub %5.1f
format neff mil %6.2f
list block domain outcome n x p lb ub neff df mil flag, noobs sepby(block domain) abbreviate(12)
export delimited using "`out'/rev3_person_goals.csv", replace

use "`out'/rev3_paired_diff.dta", clear
format pct lb ub mil mil_lb mil_ub %6.2f
list, noobs clean
export delimited using "`out'/rev3_paired_diff.csv", replace

display _n "Compare block C with JAMA: primary prevention 32.2% (29.6-34.9), ASCVD 79.6% (73.5-84.8)"
display "Differences expected: JAMA used DIQ010 only for diabetes; this file uses the"
display "prevent_build.do diabetes definition (also HbA1c >=6.5, insulin, pills) and requires"
display "Friedewald and Sampson to be available as well as Martin-Hopkins."

log close
display "DONE - see `out'/rev3_person_goals.csv"
