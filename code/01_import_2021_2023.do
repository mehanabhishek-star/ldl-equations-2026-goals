*=======================================================================*
* 01_import_2021_2023.do
* NHANES August 2021-August 2023 (fasting subsample): imports and merges the
* lipid, laboratory, examination and questionnaire files and derives the
* variables used later. Saves ldl_reclass_data.dta in the 2021-2023 data
* folder, which step 03 uses as its base.
* This step also writes descriptive tables and figures (log, CSVs and
* graphs_ldl_reclass/ in the data folder); these are intermediate outputs
* and are not reported in the paper.
* Inputs (.xpt, data/nhanes_2021_2023): DEMO_L TRIGLY_L TCHOL_L HDL_L BPQ_L
*   BMX_L DIQ_L GHB_L INS_L GLU_L LUX_L
*=======================================================================*

version 16
clear all
set more off

* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
cd "$DATA2123"

capture log close
log using "ldl_reclass_log.smcl", replace

capture mkdir "graphs_ldl_reclass"


*=======================================================================*
* PART A. IMPORT AND MERGE
*=======================================================================*
foreach f in DEMO_L TRIGLY_L TCHOL_L HDL_L BPQ_L BMX_L DIQ_L GHB_L ///
             INS_L GLU_L LUX_L {
    capture confirm file "`f'.xpt"
    if _rc {
        display as error "`f'.xpt not found in this folder"
        log close
        exit 601
    }
}

import sasxport5 "TRIGLY_L.xpt", clear
keep seqn wtsaf2yr lbxtlg lbdldl lbdldlm lbdldln
save tmp_trig, replace

import sasxport5 "TCHOL_L.xpt", clear
keep seqn lbxtc
save tmp_tc, replace

import sasxport5 "HDL_L.xpt", clear
keep seqn lbdhdd
save tmp_hdl, replace

import sasxport5 "BPQ_L.xpt", clear
keep seqn bpq101d
save tmp_bpq, replace

import sasxport5 "BMX_L.xpt", clear
keep seqn bmxbmi bmxwaist
save tmp_bmx, replace

import sasxport5 "DIQ_L.xpt", clear
keep seqn diq010 diq050
save tmp_diq, replace

import sasxport5 "GHB_L.xpt", clear
keep seqn lbxgh
save tmp_ghb, replace

import sasxport5 "INS_L.xpt", clear
keep seqn lbxin
save tmp_ins, replace

import sasxport5 "GLU_L.xpt", clear
keep seqn lbxglu
save tmp_glu, replace

import sasxport5 "LUX_L.xpt", clear
keep seqn luaxstat luxcapm luxsmed luxsiqr
save tmp_lux, replace

import sasxport5 "DEMO_L.xpt", clear
keep seqn ridageyr riagendr ridreth3 ridexprg indfmpir sdmvpsu sdmvstra
count
assert _N == 11933          // stops here if the wrong file is loaded

foreach f in tmp_trig tmp_tc tmp_hdl tmp_bpq tmp_bmx tmp_diq tmp_ghb ///
             tmp_ins tmp_glu tmp_lux {
    merge 1:1 seqn using `f', gen(m_`f')
    drop if m_`f' == 2
    drop m_`f'
}
assert _N == 11933

* ApoB is not measured in 2021-2023 - expect no apoB variable here
lookfor apo


*=======================================================================*
* PART B. VARIABLES
*=======================================================================*
label define yn 0 "No" 1 "Yes", replace

gen female = (riagendr == 2) if riagendr < .
label define sexl 0 "Male" 1 "Female", replace
label values female sexl

recode ridreth3 (1=1 "Mexican American") (2=2 "Other Hispanic") (3=3 "NH White") ///
    (4=4 "NH Black") (6=5 "NH Asian") (7=6 "Other/Multiracial"), gen(race6)

* Non-HDL cholesterol
gen nonhdl = lbxtc - lbdhdd if lbxtc < . & lbdhdd < .

* Cholesterol-lowering medication (BPQ_L - RXQ_RX_L has no drug names,
* so this is "any cholesterol medication", not statin-specific)
gen treated = (bpq101d == 1) if inlist(bpq101d, 1, 2)
label values treated yn

* Diabetes: self-report or HbA1c >=6.5
gen diabetes = (diq010 == 1) if inlist(diq010, 1, 2, 3)
replace diabetes = 1 if lbxgh >= 6.5 & lbxgh < .
label values diabetes yn

* HOMA-IR = fasting glucose (mg/dL) x insulin (uU/mL) / 405
* Not valid in people taking insulin -> set to missing
gen homa_ir = (lbxglu * lbxin) / 405 if lbxglu < . & lbxin < .
replace homa_ir = . if diq050 == 1

* High waist circumference (Obesity Society/ADA: men >=102 cm, women >=88 cm)
gen hi_waist = .
replace hi_waist = (bmxwaist >= 102) if female == 0 & bmxwaist < .
replace hi_waist = (bmxwaist >= 88)  if female == 1 & bmxwaist < .
label values hi_waist yn

* Liver elastography (valid = complete exam and IQR/median <0.30)
gen valid_te = (luaxstat == 1 & luxcapm < . & luxsmed < .)
replace valid_te = 0 if luxsiqr < . & luxsmed > 0 & luxsmed < . ///
    & (luxsiqr / luxsmed) >= 0.30
gen steatosis = (luxcapm >= 274) if valid_te == 1
gen fibrosis  = (luxsmed >= 8)   if valid_te == 1
gen cap = luxcapm if valid_te == 1
gen lsm = luxsmed if valid_te == 1
label values steatosis fibrosis yn

* Triglycerides
gen tg150 = (lbxtlg >= 150) if lbxtlg < .
label define tg150l 0 "TG <150" 1 "TG >=150", replace
label values tg150 tg150l

gen tgcat = .
replace tgcat = 1 if lbxtlg <  100
replace tgcat = 2 if lbxtlg >= 100 & lbxtlg < 150
replace tgcat = 3 if lbxtlg >= 150 & lbxtlg < 200
replace tgcat = 4 if lbxtlg >= 200 & lbxtlg < 400
label define tgcatl 1 "<100" 2 "100-149" 3 "150-199" 4 "200-399", replace
label values tgcat tgcatl
label variable tgcat "Triglycerides (mg/dL)"

* Difference from Friedewald (positive = Friedewald underestimates)
gen d_mh = lbdldlm - lbdldl
gen d_s  = lbdldln - lbdldl
label variable d_mh "Martin-Hopkins minus Friedewald LDL-C (mg/dL)"
label variable d_s  "Sampson minus Friedewald LDL-C (mg/dL)"

* How many have TG >=400 (Friedewald and Martin-Hopkins not reported)?
display "Adults 20+ in fasting subsample with TG >=400 (excluded from main analysis):"
count if ridageyr >= 20 & ridexprg != 1 & lbxtlg >= 400 & lbxtlg < .

* Analytic sample: adults 20+, not pregnant, fasting subsample,
* all three LDL estimates and non-HDL available (in practice TG <400)
gen analytic = (ridageyr >= 20 & ridexprg != 1 & lbdldl < . & lbdldlm < . ///
    & lbdldln < . & nonhdl < . & wtsaf2yr > 0 & wtsaf2yr < .)
tab analytic
tab treated tg150 if analytic == 1, missing

* Goal flags. For each LDL goal g:
*   f_goal = at goal by Friedewald
*   mh_miss / s_miss = NOT at goal by Martin-Hopkins / Sampson
*   nh_miss = NOT at non-HDL goal (g + 30)
*   any_miss = missed by at least one of the three
foreach g in 100 70 55 {
    local nh = `g' + 30
    gen f_goal`g'   = (lbdldl  <  `g') if analytic == 1
    gen mh_miss`g'  = (lbdldlm >= `g') if analytic == 1
    gen s_miss`g'   = (lbdldln >= `g') if analytic == 1
    gen nh_miss`g'  = (nonhdl  >= `nh') if analytic == 1
    gen any_miss`g' = (mh_miss`g' == 1 | s_miss`g' == 1 | nh_miss`g' == 1) if analytic == 1

    * Main subpopulation: treated, TG >=150, at goal by Friedewald
    gen sub_main`g'  = (analytic == 1 & treated == 1 & tg150 == 1 & f_goal`g' == 1)
    * Comparison: treated, TG <150, at goal by Friedewald
    gen sub_lowtg`g' = (analytic == 1 & treated == 1 & tg150 == 0 & f_goal`g' == 1)
}

* Secondary subpopulations
gen sub_treat   = (analytic == 1 & treated == 1)
gen sub_main_dm = (analytic == 1 & treated == 1 & tg150 == 1 & diabetes == 1 & f_goal70 == 1)
gen sub_untr100 = (analytic == 1 & treated == 0 & tg150 == 1 & f_goal100 == 1)
forvalues c = 1/4 {
    gen sub_tg`c' = (analytic == 1 & tgcat == `c')
}

* Phenotype analysis (Part H): ALL treated adults at Friedewald goal,
* any TG level (TG >=150 alone is too small), split into
*   0 = truly at goal (agrees with Martin-Hopkins, Sampson and non-HDL)
*   1 = falsely reassured (missed by at least one of them)
label define recll 0 "Truly at goal" 1 "Falsely reassured", replace
foreach g in 100 70 {
    gen sub_atgoal`g' = (analytic == 1 & treated == 1 & f_goal`g' == 1)
    gen recl`g' = any_miss`g' if sub_atgoal`g' == 1
    label values recl`g' recll
}

display "Unweighted n in main subpopulations (treated, TG >=150, at Friedewald goal):"
foreach g in 100 70 55 {
    quietly count if sub_main`g' == 1
    display "  LDL <`g': n = " r(N)
}

save ldl_reclass_data, replace


*=======================================================================*
* PART C. SURVEY DESIGN (fasting weight)
*=======================================================================*
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)


*=======================================================================*
* PART D. TABLE 1 - TREATED ADULTS BY TRIGLYCERIDE STATUS
*=======================================================================*
tab tg150 if sub_treat == 1

svy, subpop(sub_treat): mean ridageyr bmxbmi bmxwaist lbxtlg lbdldl lbdldlm ///
    lbdldln nonhdl lbdhdd, over(tg150)
svy, subpop(sub_treat): tab female tg150, col percent
svy, subpop(sub_treat): tab race6 tg150, col percent
svy, subpop(sub_treat): tab diabetes tg150, col percent


*=======================================================================*
* PART E. HOW FAR DOES FRIEDEWALD UNDERESTIMATE? (all analytic adults)
*=======================================================================*
svy, subpop(analytic): mean d_mh d_s, over(tgcat)

* Weighted median (IQR) difference in the main group (treated, TG >=150)
gen main_any = (analytic == 1 & treated == 1 & tg150 == 1)
_pctile d_mh [pweight = wtsaf2yr] if main_any == 1, percentiles(25 50 75)
display "Martin-Hopkins minus Friedewald, treated TG>=150: median " %4.1f r(r2) ///
    " (IQR " %4.1f r(r1) " to " %4.1f r(r3) ") mg/dL"
_pctile d_s [pweight = wtsaf2yr] if main_any == 1, percentiles(25 50 75)
display "Sampson minus Friedewald, treated TG>=150: median " %4.1f r(r2) ///
    " (IQR " %4.1f r(r1) " to " %4.1f r(r3) ") mg/dL"

* Post mean differences by TG category for Figure 2
tempname pd
postfile `pd' tgcat eq est lo hi using "tmp_diff.dta", replace
forvalues c = 1/4 {
    local e 0
    foreach v in d_mh d_s {
        local ++e
        quietly svy, subpop(sub_tg`c'): mean `v'
        local t = invttail(e(df_r), 0.025)
        post `pd' (`c') (`e') (_b[`v']) (_b[`v'] - `t' * _se[`v']) (_b[`v'] + `t' * _se[`v'])
    }
}
postclose `pd'


*=======================================================================*
* PART F. MAIN RESULT - RECLASSIFICATION AMONG ADULTS "AT GOAL"
*   Weighted % (logit-based 95% CI) and weighted number (millions)
*=======================================================================*
capture program drop post_prop
program define post_prop
    * Usage: post_prop postname "group" "label" thr eq var subpopvar
    args pname grp label thr eq var sub
    quietly count if `sub' == 1
    local n = r(N)
    quietly count if `sub' == 1 & `var' == 1
    local k = r(N)
    local p = .
    local lo = .
    local hi = .
    local pop = .
    if `k' > 0 & `k' < `n' {
        capture quietly svy, subpop(`sub'): logit `var'
        if !_rc {
            local t  = invttail(e(df_r), 0.025)
            local p  = 100 * invlogit(_b[_cons])
            local lo = 100 * invlogit(_b[_cons] - `t' * _se[_cons])
            local hi = 100 * invlogit(_b[_cons] + `t' * _se[_cons])
        }
        capture quietly svy, subpop(`sub'): total `var'
        if !_rc local pop = _b[`var'] / 1e6
    }
    else if `n' > 0 & `k' == 0 local p = 0
    else if `n' > 0 & `k' == `n' local p = 100
    if `n' < 30 display as text "Note: n = `n' (<30) for `grp' / `label' - interpret with caution"
    post `pname' ("`grp'") ("`label'") (`thr') (`eq') (`n') (`k') (`p') (`lo') (`hi') (`pop')
end

tempname ph
postfile `ph' str40 group str50 label thr eq n k pct lo hi pop_m ///
    using "ldl_reclass_results.dta", replace

foreach g in 100 70 55 {
    local nh = `g' + 30
    * Main: treated, TG >=150
    post_prop `ph' "Treated, TG>=150" "LDL-F <`g' but Martin-Hopkins >=`g'" `g' 1 mh_miss`g' sub_main`g'
    post_prop `ph' "Treated, TG>=150" "LDL-F <`g' but Sampson >=`g'"        `g' 2 s_miss`g'  sub_main`g'
    post_prop `ph' "Treated, TG>=150" "LDL-F <`g' but non-HDL >=`nh'"        `g' 3 nh_miss`g' sub_main`g'
    post_prop `ph' "Treated, TG>=150" "LDL-F <`g' but missed by any"         `g' 4 any_miss`g' sub_main`g'
    * Comparison: treated, TG <150
    post_prop `ph' "Treated, TG<150"  "LDL-F <`g' but Martin-Hopkins >=`g'" `g' 1 mh_miss`g' sub_lowtg`g'
    post_prop `ph' "Treated, TG<150"  "LDL-F <`g' but non-HDL >=`nh'"        `g' 3 nh_miss`g' sub_lowtg`g'
}

* Secondary: treated with diabetes and TG >=150 (high-risk goal <70)
post_prop `ph' "Treated, TG>=150, diabetes" "LDL-F <70 but Martin-Hopkins >=70" 70 1 mh_miss70 sub_main_dm
post_prop `ph' "Treated, TG>=150, diabetes" "LDL-F <70 but non-HDL >=100"       70 3 nh_miss70 sub_main_dm

* Secondary: untreated adults with TG >=150 (primary prevention, <100)
post_prop `ph' "Untreated, TG>=150" "LDL-F <100 but Martin-Hopkins >=100" 100 1 mh_miss100 sub_untr100
post_prop `ph' "Untreated, TG>=150" "LDL-F <100 but non-HDL >=130"        100 3 nh_miss100 sub_untr100

postclose `ph'


*=======================================================================*
* PART G. HOW FAR ABOVE GOAL ARE THE RECLASSIFIED? (noise check)
*   G1. Weighted median (IQR) distance above goal
*   G2. Reclassification only counted if >=5 or >=10 mg/dL above goal
*       (eq 5 = MH >= goal+5, eq 6 = MH >= goal+10, eq 7 = non-HDL >= goal+10)
*=======================================================================*

*--- G1. Distance above goal among the reclassified (treated, TG >=150)
foreach g in 100 70 55 {
    local nh = `g' + 30
    gen gap_mh`g' = lbdldlm - `g'  if sub_main`g' == 1 & mh_miss`g' == 1
    gen gap_nh`g' = nonhdl  - `nh' if sub_main`g' == 1 & nh_miss`g' == 1

    quietly count if gap_mh`g' < .
    if r(N) >= 2 {
        _pctile gap_mh`g' [pweight = wtsaf2yr], percentiles(25 50 75)
        display "Goal <`g': reclassified by Martin-Hopkins are a median " %4.1f r(r2) ///
            " mg/dL above goal (IQR " %4.1f r(r1) " to " %4.1f r(r3) ")"
    }
    quietly count if gap_nh`g' < .
    if r(N) >= 2 {
        _pctile gap_nh`g' [pweight = wtsaf2yr], percentiles(25 50 75)
        display "Goal <`g': reclassified by non-HDL are a median " %4.1f r(r2) ///
            " mg/dL above the non-HDL goal (IQR " %4.1f r(r1) " to " %4.1f r(r3) ")"
    }

    gen mh_m5_`g'  = (lbdldlm >= `g' + 5)   if analytic == 1
    gen mh_m10_`g' = (lbdldlm >= `g' + 10)  if analytic == 1
    gen nh_m10_`g' = (nonhdl  >= `nh' + 10) if analytic == 1
}

*--- G2. Reclassification with a safety margin
tempname pm
postfile `pm' str40 group str50 label thr eq n k pct lo hi pop_m ///
    using "tmp_margin.dta", replace
foreach g in 100 70 55 {
    local g5   = `g' + 5
    local g10  = `g' + 10
    local nh10 = `g' + 40
    post_prop `pm' "Treated, TG>=150" "LDL-F <`g' but Martin-Hopkins >=`g5'"  `g' 5 mh_m5_`g'  sub_main`g'
    post_prop `pm' "Treated, TG>=150" "LDL-F <`g' but Martin-Hopkins >=`g10'" `g' 6 mh_m10_`g' sub_main`g'
    post_prop `pm' "Treated, TG>=150" "LDL-F <`g' but non-HDL >=`nh10'"       `g' 7 nh_m10_`g' sub_main`g'
}
postclose `pm'

*--- Combine Parts F and G into one results table
preserve
use "ldl_reclass_results.dta", clear
append using "tmp_margin.dta"
sort group thr eq
format pct lo hi %5.1f
format pop_m %5.2f
save "ldl_reclass_results.dta", replace
list group label n k pct lo hi pop_m, noobs sepby(group thr) abbreviate(20)
export delimited using "ldl_reclass_results.csv", replace
restore


*=======================================================================*
* PART H. WHO IS FALSELY REASSURED? METABOLIC AND LIVER PHENOTYPE
*   All treated adults at Friedewald goal (<70 primary, <100 secondary),
*   falsely reassured vs truly at goal. Descriptive only - small groups.
*   Note: TG is part of how reclassification happens, so a TG difference
*   is expected; the question is what ELSE travels with it.
*=======================================================================*
foreach g in 70 100 {
    display _n "=== Goal <`g': unweighted counts ==="
    tab recl`g' if sub_atgoal`g' == 1
    tab recl`g' valid_te if sub_atgoal`g' == 1
}

capture program drop post_cmp
program define post_cmp
    * Usage: post_cmp postname thr "label" var type subpopvar groupvar
    *   type: bin  = % (CI clipped to 0-100)
    *         cont = mean
    *         log  = geometric mean (skewed variables)
    args pname thr label var type sub grp
    quietly count if `sub' == 1 & `grp' == 0 & `var' < .
    local n0 = r(N)
    quietly count if `sub' == 1 & `grp' == 1 & `var' < .
    local n1 = r(N)
    forvalues j = 0/1 {
        local est`j' = .
        local lo`j'  = .
        local hi`j'  = .
    }
    local p = .
    tempvar y
    if "`type'" == "bin"  quietly gen `y' = 100 * `var'
    if "`type'" == "cont" quietly gen `y' = `var'
    if "`type'" == "log"  quietly gen `y' = ln(`var') if `var' > 0 & `var' < .
    if `n0' >= 2 & `n1' >= 2 {
        capture quietly svy, subpop(`sub'): regress `y' ibn.`grp', noconstant
        if !_rc {
            local t = invttail(e(df_r), 0.025)
            forvalues j = 0/1 {
                local b = _b[`j'.`grp']
                local s = _se[`j'.`grp']
                if "`type'" == "log" {
                    local est`j' = exp(`b')
                    local lo`j'  = exp(`b' - `t' * `s')
                    local hi`j'  = exp(`b' + `t' * `s')
                }
                else {
                    local est`j' = `b'
                    local lo`j'  = `b' - `t' * `s'
                    local hi`j'  = `b' + `t' * `s'
                }
                if "`type'" == "bin" {
                    local lo`j' = max(0, `lo`j'')
                    local hi`j' = min(100, `hi`j'')
                }
            }
            capture quietly test 0.`grp' = 1.`grp'
            if !_rc local p = r(p)
        }
    }
    post `pname' (`thr') ("`label'") ("`type'") (`n0') (`est0') (`lo0') (`hi0') ///
        (`n1') (`est1') (`lo1') (`hi1') (`p')
end

tempname pc
postfile `pc' thr str40 label str4 type n0 est0 lo0 hi0 n1 est1 lo1 hi1 p ///
    using "ldl_reclass_phenotype.dta", replace
foreach g in 70 100 {
    local s sub_atgoal`g'
    local r recl`g'
    post_cmp `pc' `g' "TG >=150 mg/dL (%)"           tg150     bin  `s' `r'
    post_cmp `pc' `g' "Triglycerides (geo. mean)"    lbxtlg    log  `s' `r'
    post_cmp `pc' `g' "HDL-C (mg/dL)"                lbdhdd    cont `s' `r'
    post_cmp `pc' `g' "Age (years)"                  ridageyr  cont `s' `r'
    post_cmp `pc' `g' "Female (%)"                   female    bin  `s' `r'
    post_cmp `pc' `g' "Diabetes (%)"                 diabetes  bin  `s' `r'
    post_cmp `pc' `g' "HbA1c (%)"                    lbxgh     cont `s' `r'
    post_cmp `pc' `g' "HOMA-IR (geo. mean)"          homa_ir   log  `s' `r'
    post_cmp `pc' `g' "BMI (kg/m2)"                  bmxbmi    cont `s' `r'
    post_cmp `pc' `g' "Waist (cm)"                   bmxwaist  cont `s' `r'
    post_cmp `pc' `g' "High waist (%)"               hi_waist  bin  `s' `r'
    post_cmp `pc' `g' "CAP (dB/m)"                   cap       cont `s' `r'
    post_cmp `pc' `g' "Steatosis, CAP >=274 (%)"     steatosis bin  `s' `r'
    post_cmp `pc' `g' "LSM (kPa, geo. mean)"         lsm       log  `s' `r'
    post_cmp `pc' `g' "Fibrosis, LSM >=8 kPa (%)"    fibrosis  bin  `s' `r'
}
postclose `pc'

preserve
use "ldl_reclass_phenotype.dta", clear
format est0 lo0 hi0 est1 lo1 hi1 %6.1f
format p %5.3f
display _n "Phenotype: 0 = truly at goal, 1 = falsely reassured (weighted)"
list thr label n0 est0 lo0 hi0 n1 est1 lo1 hi1 p, noobs sepby(thr) abbreviate(12)
export delimited using "ldl_reclass_phenotype.csv", replace
restore


*=======================================================================*
* PART I. FIGURES
*=======================================================================*

*--- Figure 1: % reclassified, treated adults with TG >=150
preserve
use "ldl_reclass_results.dta", clear
keep if group == "Treated, TG>=150" & eq <= 3
gen thr_pos = cond(thr == 100, 1, cond(thr == 70, 2, 3))
gen x = thr_pos + (eq - 2) * 0.24
twoway (bar pct x if eq == 1, barwidth(0.22) color(navy)) ///
       (bar pct x if eq == 2, barwidth(0.22) color(maroon)) ///
       (bar pct x if eq == 3, barwidth(0.22) color(forest_green)) ///
       (rcap lo hi x, lcolor(black)), ///
    xlabel(1 "LDL-C goal <100" 2 "<70" 3 "<55", noticks) xtitle("") ///
    yscale(range(0 100)) ylabel(0(20)100, angle(0)) ///
    ytitle("% NOT at goal by alternative measure") ///
    title("Adults on cholesterol medication with TG >=150 mg/dL" ///
          "who are at goal by Friedewald LDL-C", size(medsmall)) ///
    legend(order(1 "Martin-Hopkins LDL-C" 2 "Sampson LDL-C" 3 "Non-HDL-C (goal + 30)") ///
           rows(1) size(vsmall)) ///
    note("NHANES 2021-2023 fasting subsample, survey-weighted; bars = 95% CI.", size(vsmall)) ///
    graphregion(color(white))
graph export "graphs_ldl_reclass/fig1_reclassification.png", width(2400) replace
restore

*--- Figure 2: Friedewald underestimation by TG category
preserve
use "tmp_diff.dta", clear
gen x = tgcat + (eq - 1.5) * 0.25
twoway (rcap lo hi x if eq == 1, lcolor(navy)) ///
       (scatter est x if eq == 1, msymbol(D) mcolor(navy)) ///
       (rcap lo hi x if eq == 2, lcolor(maroon)) ///
       (scatter est x if eq == 2, msymbol(O) mcolor(maroon)), ///
    yline(0, lpattern(dash) lcolor(gs8)) ///
    xlabel(1 "<100" 2 "100-149" 3 "150-199" 4 "200-399", noticks) ///
    xtitle("Fasting triglycerides (mg/dL)") ///
    ytitle("Mean difference vs Friedewald (mg/dL)") ylabel(, angle(0)) ///
    title("How much Friedewald underestimates LDL-C", size(medsmall)) ///
    legend(order(2 "Martin-Hopkins minus Friedewald" 4 "Sampson minus Friedewald") ///
           rows(1) size(small)) ///
    note("Adults 20+, NHANES 2021-2023 fasting subsample, survey-weighted; 95% CI.", size(vsmall)) ///
    graphregion(color(white))
graph export "graphs_ldl_reclass/fig2_difference_by_tg.png", width(2400) replace
restore

*--- Figure 3: Friedewald vs Martin-Hopkins, treated adults
twoway (scatter lbdldlm lbdldl if sub_treat == 1 & tg150 == 0, ///
            msize(vsmall) mcolor(gs10%50)) ///
       (scatter lbdldlm lbdldl if sub_treat == 1 & tg150 == 1, ///
            msize(vsmall) mcolor(maroon%60)) ///
       (function y = x, range(20 250) lpattern(dash) lcolor(black)), ///
    xline(70, lpattern(dot) lcolor(gs6)) yline(70, lpattern(dot) lcolor(gs6)) ///
    xtitle("Friedewald LDL-C (mg/dL)") ytitle("Martin-Hopkins LDL-C (mg/dL)") ///
    ylabel(, angle(0)) ///
    title("Adults on cholesterol medication", size(medsmall)) ///
    legend(order(1 "TG <150" 2 "TG >=150" 3 "Line of equality") rows(1) size(small)) ///
    note("Points left of x = 70 and above y = 70 are at goal (<70) by Friedewald" ///
         "but not by Martin-Hopkins. Unweighted.", size(vsmall)) ///
    graphregion(color(white))
graph export "graphs_ldl_reclass/fig3_scatter_f_vs_mh.png", width(2400) replace

*--- Figure 4: phenotype of the falsely reassured (goal <70, % variables)
preserve
use "ldl_reclass_phenotype.dta", clear
keep if thr == 70 & type == "bin" & label != "Female (%)"
gen lab2 = subinstr(label, " (%)", "", .)
local N = _N
gen pos = `N' - _n + 1
gen y0 = pos + 0.15
gen y1 = pos - 0.15
local ylab ""
forvalues i = 1/`N' {
    local ps = pos[`i']
    local l  = lab2[`i']
    local ylab `"`ylab' `ps' `"`l'"'"'
}
twoway (rcap lo0 hi0 y0, horizontal lcolor(navy)) ///
       (scatter y0 est0, msymbol(O) mcolor(navy)) ///
       (rcap lo1 hi1 y1, horizontal lcolor(maroon)) ///
       (scatter y1 est1, msymbol(D) mcolor(maroon)), ///
    ylabel(`ylab', angle(0) noticks labsize(small) nogrid) ///
    yscale(range(0.5 `=`N' + 0.5')) ///
    xlabel(0(20)100) xscale(range(0 100)) xtitle("Weighted %") ytitle("") ///
    title("Who is falsely reassured by Friedewald LDL-C <70?", size(medsmall)) ///
    legend(order(2 "Truly at goal" 4 "Falsely reassured") rows(1) size(small)) ///
    note("Adults on cholesterol medication with Friedewald LDL-C <70 mg/dL." ///
         "Falsely reassured = not at goal by Martin-Hopkins, Sampson or non-HDL-C." ///
         "NHANES 2021-2023, survey-weighted; 95% CI. Liver measures: valid elastography only.", ///
         size(vsmall)) ///
    graphregion(color(white))
graph export "graphs_ldl_reclass/fig4_phenotype.png", width(2400) replace
restore

*--- Figure 5: how far above goal (Martin-Hopkins), treated at Friedewald <70
gen gap_all70 = lbdldlm - 70 if sub_atgoal70 == 1 & mh_miss70 == 1
quietly count if gap_all70 < .
if r(N) >= 2 {
    histogram gap_all70, width(5) start(0) frequency color(navy%70) lcolor(white) ///
        xtitle("Martin-Hopkins LDL-C minus 70 (mg/dL)") ytitle("Number of participants") ///
        ylabel(, angle(0)) ///
        title("How far above the <70 goal are the reclassified?", size(medsmall)) ///
        note("Adults on cholesterol medication at goal (<70) by Friedewald but not by" ///
             "Martin-Hopkins. Unweighted counts.", size(vsmall)) ///
        graphregion(color(white))
    graph export "graphs_ldl_reclass/fig5_gap_above_goal.png", width(2400) replace
}


*=======================================================================*
* CLEAN UP
*=======================================================================*
save ldl_reclass_data, replace

foreach f in tmp_trig tmp_tc tmp_hdl tmp_bpq tmp_bmx tmp_diq tmp_ghb ///
             tmp_ins tmp_glu tmp_lux tmp_diff tmp_margin {
    capture erase "`f'.dta"
}

log close
display "DONE - see ldl_reclass_log.smcl, ldl_reclass_results.csv,"
display "       ldl_reclass_phenotype.csv and graphs_ldl_reclass/"
*=======================================================================*
* END
*=======================================================================*
