*=======================================================================*
* _moverkg.do - programs moverkg and rowstats (needs _kgprop.do)
*
* moverkg up dn dom [struct]
*   95% CI for d = p(up) - p(down) by MOVER (Newcombe 1998; Zou & Donner
*   2008), combining Korn-Graubard limits for each direction with their
*   design-based correlation. struct = 1: 'down' is zero by construction.
*   Returns r(est) r(lb) r(ub) (percentage points), r(x1) r(x2)
*   (participants moving up / down), r(rho).
*
* rowstats handle block label dom fA mA sA nA
*   One row of the main table for domain dom, given 0/1 "above goal"
*   indicators by Friedewald (fA), Martin-Hopkins (mA), Sampson (sA) and
*   the paired non-HDL-C goal (nA).
*=======================================================================*
capture program drop moverkg
program define moverkg, rclass
    args up dn sub struct
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
        matrix V__ = e(V)
        if V__[1,1] > 0 & V__[2,2] > 0 local rho = V__[1,2] / sqrt(V__[1,1] * V__[2,2])
    }
    local d = `p1' - `p2'
    local L = `d' - sqrt((`p1' - `l1')^2 + (`u2' - `p2')^2 - 2*`rho'*(`p1' - `l1')*(`u2' - `p2'))
    local U = `d' + sqrt((`u1' - `p1')^2 + (`p2' - `l2')^2 - 2*`rho'*(`u1' - `p1')*(`p2' - `l2'))
    return scalar est = 100 * `d'
    return scalar lb  = 100 * `L'
    return scalar ub  = 100 * `U'
    return scalar rho = `rho'
    return scalar x1  = `x1'
    return scalar x2  = `x2'
end

capture program drop rowstats
program define rowstats
    args h block label dom fA mA sA nA
    tempvar up dn eF gup rup rdn
    quietly gen byte `up'  = (`mA' == 1 & `fA' == 0) if `mA' < . & `fA' < .
    quietly gen byte `dn'  = (`mA' == 0 & `fA' == 1) if `mA' < . & `fA' < .
    quietly gen byte `eF'  = (`fA' == 1 | `nA' == 1) if `fA' < . & `nA' < .
    quietly gen byte `gup' = (`fA' == 0 & `nA' == 1) if `fA' < . & `nA' < .
    quietly gen byte `rup' = (`mA' == 1 & `eF' == 0) if `mA' < . & `eF' < .
    quietly gen byte `rdn' = (`mA' == 0 & `eF' == 1) if `mA' < . & `eF' < .
    foreach k in F E M S {
        if "`k'" == "F" local v `fA'
        if "`k'" == "E" local v `eF'
        if "`k'" == "M" local v `mA'
        if "`k'" == "S" local v `sA'
        kgprop `v', sub(`dom')
        local n    = r(n)
        local p`k' = r(p)
        local l`k' = r(lb)
        local u`k' = r(ub)
        local f`k' "`r(flag)'"
    }
    moverkg `up' `dn' `dom' 0
    local g  = r(est)
    local gl = r(lb)
    local gu = r(ub)
    local gx1 = r(x1)
    local gx2 = r(x2)
    kgprop `gup', sub(`dom')
    local c  = r(p)
    local cl = r(lb)
    local cu = r(ub)
    local cx = r(x)
    moverkg `rup' `rdn' `dom' 0
    local q  = r(est)
    local ql = r(lb)
    local qu = r(ub)
    local qx1 = r(x1)
    local qx2 = r(x2)
    post `h' ("`block'") ("`label'") (`n') (`pF') (`lF') (`uF') (`pE') (`lE') (`uE') ///
        (`pM') (`lM') (`uM') (`pS') (`lS') (`uS') (`g') (`gl') (`gu') (`gx1') (`gx2') ///
        (`c') (`cl') (`cu') (`cx') (`q') (`ql') (`qu') (`qx1') (`qx2') ///
        ("`fF'") ("`fM'")
end
