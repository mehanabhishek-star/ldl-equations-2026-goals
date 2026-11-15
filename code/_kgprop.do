*=======================================================================*
* _kgprop.do - program kgprop
*   kgprop var, sub(domain)
*   Survey-weighted proportion with a Korn-Graubard 95% CI (degrees-of-
*   freedom-adjusted effective sample size) and NCHS presentation-standard
*   flags (Parker et al 2017). Returns r(n) r(x) r(p) r(lb) r(ub) r(flag).
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
