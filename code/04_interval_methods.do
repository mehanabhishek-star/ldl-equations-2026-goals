*=======================================================================*
* 04_interval_methods.do
* 95% CIs for every paired difference (adults with a defined goal), by:
*   A. z interval             est +/- 1.96 SE
*   B. t, design df            est +/- t(15) SE (15 strata, 30 PSUs)
*   C. t, domain df            PSUs - strata containing the domain
*   D. t, Satterthwaite df     effective df from the stratum variance
*                              contributions
*   E. MOVER-KG                difference of the two directions of
*                              reclassification (up - down), combining their
*                              Korn-Graubard limits with the design-based
*                              correlation (Newcombe 1998; Zou & Donner 2008)
* The non-HDL-C gain has no 'down' direction by construction and is exactly
* 0 where TG <150 mg/dL. A manual linearization SE is checked against
* Stata's SE before D is used.
* Input: results/rev3_analysis_data.dta   Output: results/rev5_p2_ci.csv
*=======================================================================*
version 17
clear all
set more off
set varabbrev off
* Paths are set by config.do; run from the repository root (see README).
if "$RESULTS" == "" do "config.do"
local out "$RESULTS"
capture log close
log using "`out'/rev5_p2_ci.log", text replace
do "$CODE/_kgprop.do"

use "`out'/rev3_analysis_data.dta", clear
svyset sdmvpsu [pweight = wtsaf2yr], strata(sdmvstra) vce(linearized) singleunit(centered)

*--- domains
gen byte d_def      = (d_bm_all == 1 & formal_goal == 1)
gen byte d_def_prim = (d_def == 1 & inlist(goalgrp, 2, 3))
gen byte d_def_sec  = (d_def == 1 & inlist(goalgrp, 4, 5))
gen byte d_def_t150 = (d_def == 1 & treated == 1 & tg150 == 1)
gen byte d_def_tlow = (d_def == 1 & treated == 1 & tg150 == 0)
gen byte d_def_untr = (d_def == 1 & treated == 0)
forvalues c = 1/4 {
    gen byte d_def_tg`c' = (d_def == 1 & tgcat == `c')
}

*--- differences and their two directions (each difference = up - down)
gen byte either_f = (f_above == 1 | nh_miss == 1) if mh_miss < .
gen gain_nh = either_f - f_above
gen resid   = above_mh - either_f
* dif_mh: MH vs F
gen byte up_mh = (above_mh == 1 & f_at == 1)   if mh_miss < .
gen byte dn_mh = (above_mh == 0 & f_at == 0)   if mh_miss < .
* dif_s: Sampson vs F
gen byte up_s  = (above_s == 1 & f_at == 1)    if mh_miss < .
gen byte dn_s  = (above_s == 0 & f_at == 0)    if mh_miss < .
* gain_nh: non-HDL-C flags someone F calls at goal (never negative)
gen byte up_g  = (f_at == 1 & nh_miss == 1)    if mh_miss < .
gen byte dn_g  = 0                             if mh_miss < .
* resid: above by MH but at goal by F and non-HDL-C, minus the reverse
gen byte up_r  = (above_mh == 1 & either_f == 0) if mh_miss < .
gen byte dn_r  = (above_mh == 0 & either_f == 1) if mh_miss < .
foreach p in mh s g r {
    if "`p'" == "mh" local v dif_mh
    if "`p'" == "s"  local v dif_s
    if "`p'" == "g"  local v gain_nh
    if "`p'" == "r"  local v resid
    capture assert `v' == up_`p' - dn_`p' if mh_miss < .
    if _rc display as error "component check failed for `v'"
    assert `v' == up_`p' - dn_`p' if mh_miss < .
}

*-----------------------------------------------------------------------*
* Satterthwaite effective df for a domain mean (with SE check)
*-----------------------------------------------------------------------*
capture program drop satt
program define satt, rclass
    syntax varname, SUB(varname)
    quietly svy, subpop(`sub'): mean `varlist'
    local b  = _b[`varlist']
    local se = _se[`varlist']
    preserve
    quietly keep if wtsaf2yr < .
    tempvar dom z
    quietly gen byte `dom' = (`sub' == 1 & `varlist' < .)
    quietly summarize wtsaf2yr if `dom', meanonly
    local W = r(sum)
    quietly gen double `z' = cond(`dom', wtsaf2yr * (`varlist' - `b') / `W', 0)
    quietly collapse (sum) zz = `z', by(sdmvstra sdmvpsu)
    quietly bysort sdmvstra: gen nh = _N
    quietly by sdmvstra: egen double zbar = mean(zz)
    quietly gen double dev2 = (zz - zbar)^2
    quietly collapse (sum) dev2 (first) nh, by(sdmvstra)
    quietly count if nh < 2
    local nsing = r(N)
    quietly gen double vh  = nh / (nh - 1) * dev2
    quietly gen double vh2 = vh^2 / (nh - 1)
    quietly summarize vh, meanonly
    local V = r(sum)
    quietly summarize vh2, meanonly
    local den = r(sum)
    * how concentrated is the variance? share carried by the top 3 strata
    gsort -vh
    quietly gen double cum = sum(vh) / `V'
    local top3 = cum[3]
    restore
    return scalar est   = `b'
    return scalar se    = `se'
    return scalar se_man = sqrt(`V')
    return scalar df_s  = cond(`den' > 0, `V'^2 / `den', .)
    return scalar top3  = `top3'
    return scalar nsing = `nsing'
end

*-----------------------------------------------------------------------*
* MOVER with Korn-Graubard limits for d = p(up) - p(down)
*-----------------------------------------------------------------------*
capture program drop moverkg
program define moverkg, rclass
    args up dn sub struct
    * struct = 1: the down direction is zero BY CONSTRUCTION (gain_nh can
    * never be negative), so it gets no sampling interval
    kgprop `up', sub(`sub')
    local p1 = r(p)/100
    local l1 = r(lb)/100
    local u1 = r(ub)/100
    local x1 = r(x)
    kgprop `dn', sub(`sub')
    local p2 = r(p)/100
    local l2 = r(lb)/100
    local u2 = r(ub)/100
    local x2 = r(x)
    if "`struct'" == "1" {
        local p2 = 0
        local l2 = 0
        local u2 = 0
    }
    local rho = 0
    if `x1' > 0 & `x2' > 0 {
        quietly svy, subpop(`sub'): mean `up' `dn'
        matrix V = e(V)
        if V[1,1] > 0 & V[2,2] > 0 local rho = V[1,2] / sqrt(V[1,1] * V[2,2])
    }
    local d = `p1' - `p2'
    local L = `d' - sqrt((`p1' - `l1')^2 + (`u2' - `p2')^2 - 2*`rho'*(`p1' - `l1')*(`u2' - `p2'))
    local U = `d' + sqrt((`u1' - `p1')^2 + (`p2' - `l2')^2 - 2*`rho'*(`u1' - `p1')*(`p2' - `l2'))
    return scalar lb  = 100 * `L'
    return scalar ub  = 100 * `U'
    return scalar rho = `rho'
    return scalar x1  = `x1'
    return scalar x2  = `x2'
end

*-----------------------------------------------------------------------*
* Domain df (as in _kgprop.do)
*-----------------------------------------------------------------------*
capture program drop domdf
program define domdf, rclass
    syntax varname, SUB(varname)
    tempvar dom tags tagp
    quietly gen byte `dom' = (`sub' == 1 & `varlist' < .)
    quietly egen `tags' = tag(sdmvstra)         if `dom'
    quietly egen `tagp' = tag(sdmvstra sdmvpsu) if `dom'
    quietly count if `tags' == 1
    local nstr = r(N)
    quietly count if `tagp' == 1
    return scalar df = r(N) - `nstr'
    * df counted only over PSUs that contain a reclassified person
    quietly count if `dom'
    return scalar n = r(N)
end

*-----------------------------------------------------------------------*
* Run
*-----------------------------------------------------------------------*
tempname P
postfile `P' str12 domain str8 comp n n_up n_dn est se se_man ///
    z_lb z_ub t15_lb t15_ub df_dom tdom_lb tdom_ub df_satt tsat_lb tsat_ub ///
    mover_lb mover_ub rho top3 nsing using "`out'/rev5_p2_ci.dta", replace

foreach d in d_bm_all d_def d_def_prim d_def_sec d_def_t150 d_def_tlow d_def_untr ///
             d_def_tg1 d_def_tg2 d_def_tg3 d_def_tg4 {
    foreach p in mh s g r {
        if "`p'" == "mh" local v dif_mh
        if "`p'" == "s"  local v dif_s
        if "`p'" == "g"  local v gain_nh
        if "`p'" == "r"  local v resid
        satt `v', sub(`d')
        local b = r(est)
        local s = r(se)
        local sm = r(se_man)
        local dfs = r(df_s)
        local top3 = r(top3)
        local nsing = r(nsing)
        domdf `v', sub(`d')
        local dfd = r(df)
        local n = r(n)
        local st = cond("`p'" == "g", 1, 0)
        moverkg up_`p' dn_`p' `d' `st'
        local ml = r(lb)
        local mu = r(ub)
        local rho = r(rho)
        local x1 = r(x1)
        local x2 = r(x2)
        local t15  = invttail(15, 0.025)
        local tdom = invttail(max(`dfd', 1), 0.025)
        if `dfs' < . & `dfs' >= 1 local tsat = invttail(`dfs', 0.025)
        else local tsat = .
        post `P' ("`d'") ("`v'") (`n') (`x1') (`x2') (100*`b') (100*`s') (100*`sm') ///
            (100*(`b' - 1.96*`s')) (100*(`b' + 1.96*`s')) ///
            (100*(`b' - `t15'*`s')) (100*(`b' + `t15'*`s')) ///
            (`dfd') (100*(`b' - `tdom'*`s')) (100*(`b' + `tdom'*`s')) ///
            (`dfs') (100*(`b' - `tsat'*`s')) (100*(`b' + `tsat'*`s')) ///
            (`ml') (`mu') (`rho') (`top3') (`nsing')
    }
}
postclose `P'

use "`out'/rev5_p2_ci.dta", clear
* the manual linearization SE must reproduce Stata's before df_satt is trusted
gen double se_err = abs(se - se_man)
summarize se_err
assert se_err < 0.01 if se > 0
format est se se_man z_lb-tsat_ub mover_lb mover_ub %6.2f
format df_satt rho top3 %5.2f
list domain comp n n_up n_dn est t15_lb t15_ub df_dom df_satt tsat_lb tsat_ub ///
    mover_lb mover_ub top3, noobs sepby(domain) abbreviate(8)
export delimited using "`out'/rev5_p2_ci.csv", replace
log close
display "DONE"
