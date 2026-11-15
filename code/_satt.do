*=======================================================================*
* _satt.do - program satt (Satterthwaite effective degrees of freedom)
*   satt y, sub(domain) wt(weightvar)
*   Survey-weighted domain mean of y with linearization SE, plus the
*   Satterthwaite effective df from the stratum variance contributions,
*   the share of variance carried by the 3 largest strata, and t intervals
*   with design df and effective df.
*   Returns r(est) r(se) r(se_man) r(df_design) r(df_satt) r(top3)
*           r(lb_t) r(ub_t) r(lb_s) r(ub_s)   (est, se, limits x100)
*=======================================================================*
capture program drop satt
program define satt, rclass
    syntax varname, SUB(varname) WT(varname)
    quietly svy, subpop(`sub'): mean `varlist'
    local b  = _b[`varlist']
    local se = _se[`varlist']
    local dfd = e(df_r)
    preserve
    quietly keep if `wt' < . & `wt' > 0
    tempvar dom z
    quietly gen byte `dom' = (`sub' == 1 & `varlist' < .)
    quietly summarize `wt' if `dom', meanonly
    local W = r(sum)
    quietly gen double `z' = cond(`dom', `wt' * (`varlist' - `b') / `W', 0)
    quietly collapse (sum) zz = `z', by(sdmvstra sdmvpsu)
    quietly bysort sdmvstra: gen nh = _N
    quietly by sdmvstra: egen double zbar = mean(zz)
    quietly gen double dev2 = (zz - zbar)^2
    quietly collapse (sum) dev2 (first) nh, by(sdmvstra)
    quietly gen double vh  = nh / (nh - 1) * dev2
    quietly gen double vh2 = vh^2 / (nh - 1)
    quietly summarize vh, meanonly
    local V = r(sum)
    quietly summarize vh2, meanonly
    local den = r(sum)
    gsort -vh
    quietly gen double cum = sum(vh) / `V'
    local top3 = cum[3]
    restore
    local dfs = cond(`den' > 0, `V'^2 / `den', .)
    local td = invttail(`dfd', 0.025)
    local ts = invttail(`dfs', 0.025)
    return scalar est = 100 * `b'
    return scalar se = 100 * `se'
    return scalar se_man = 100 * sqrt(`V')
    return scalar df_design = `dfd'
    return scalar df_satt = `dfs'
    return scalar top3 = `top3'
    return scalar lb_t = 100 * (`b' - `td' * `se')
    return scalar ub_t = 100 * (`b' + `td' * `se')
    return scalar lb_s = 100 * (`b' - `ts' * `se')
    return scalar ub_s = 100 * (`b' + `ts' * `se')
end
