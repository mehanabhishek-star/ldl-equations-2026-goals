*=======================================================================*
* _prevent_programs.do - PREVENT calculators (same equations as
* 02_prevent_scores.do), used for the 2017-March 2020 cycle in step 09.
*=======================================================================*

capture program drop prevent_calc
program define prevent_calc
    capture drop t_* lp_* r10_ascvd r10_cvd r30_ascvd r10_ascvd_full

    * --- centred / spline terms (Khan 2024, Table S12) ---
    gen double t_age      = (p_age - 55) / 10                if p_age < .
    gen double t_age2     = t_age^2
    gen double t_nhdl     = (p_tc - p_hdl) * 0.02586 - 3.5   if p_tc < . & p_hdl < .
    gen double t_hdl      = (p_hdl * 0.02586 - 1.3) / 0.3    if p_hdl < .
    gen double t_sbp_lo   = (min(p_sbp, 110) - 110) / 20     if p_sbp < .
    gen double t_sbp_hi   = (max(p_sbp, 110) - 130) / 20     if p_sbp < .
    gen double t_dm       = p_dm
    gen double t_smk      = p_smk
    gen double t_bmi_lo   = (min(p_bmi, 30) - 25) / 5        if p_bmi < .
    gen double t_bmi_hi   = (max(p_bmi, 30) - 30) / 5        if p_bmi < .
    gen double t_egfr_lo  = (min(p_egfr, 60) - 60) / -15     if p_egfr < .
    gen double t_egfr_hi  = (max(p_egfr, 60) - 90) / -15     if p_egfr < .
    gen double t_bptx     = p_bptx
    gen double t_statin   = p_statin
    gen double t_bptx_sbp    = t_bptx * t_sbp_hi
    gen double t_statin_nhdl = t_statin * t_nhdl
    gen double t_age_nhdl = t_age * t_nhdl
    gen double t_age_hdl  = t_age * t_hdl
    gen double t_age_sbp  = t_age * t_sbp_hi
    gen double t_age_dm   = t_age * t_dm
    gen double t_age_smk  = t_age * t_smk
    gen double t_age_bmi  = t_age * t_bmi_hi
    gen double t_age_egfr = t_age * t_egfr_lo
    * full-model extras (SDI needs ZIP code -> always "missing" in NHANES)
    gen double t_sdi46    = 0
    gen double t_sdi710   = 0
    gen double t_sdimiss  = 1
    gen double t_lnuacr   = cond(p_uacr < ., ln(p_uacr), 0)
    gen double t_uacrmiss = (p_uacr >= .)
    gen double t_a1c_dm   = cond(p_a1c < . & p_dm == 1, p_a1c - 5.3, 0)
    gen double t_a1c_nodm = cond(p_a1c < . & p_dm == 0, p_a1c - 5.3, 0)
    gen double t_a1cmiss  = (p_a1c >= .)

    * --- linear predictors (coefficients generated from preventr) ---
    * base_10yr - ascvd
    gen double lp_b10_ascvd = .
    replace lp_b10_ascvd = -3.819975 + (0.719883*t_age) + (0.1176967*t_nhdl) ///
        + (-0.151185*t_hdl) + (-0.0835358*t_sbp_lo) + (0.3592852*t_sbp_hi) ///
        + (0.8348585*t_dm) + (0.4831078*t_smk) + (0.4864619*t_egfr_lo) ///
        + (0.0397779*t_egfr_hi) + (0.2265309*t_bptx) + (-0.0592374*t_statin) ///
        + (-0.0395762*t_bptx_sbp) + (0.0844423*t_statin_nhdl) ///
        + (-0.0567839*t_age_nhdl) + (0.0325692*t_age_hdl) + (-0.1035985*t_age_sbp) ///
        + (-0.2417542*t_age_dm) + (-0.0791142*t_age_smk) + (-0.1671492*t_age_egfr) if female == 1
    replace lp_b10_ascvd = -3.500655 + (0.7099847*t_age) + (0.1658663*t_nhdl) ///
        + (-0.1144285*t_hdl) + (-0.2837212*t_sbp_lo) + (0.3239977*t_sbp_hi) ///
        + (0.7189597*t_dm) + (0.3956973*t_smk) + (0.3690075*t_egfr_lo) ///
        + (0.0203619*t_egfr_hi) + (0.2036522*t_bptx) + (-0.0865581*t_statin) ///
        + (-0.0322916*t_bptx_sbp) + (0.114563*t_statin_nhdl) ///
        + (-0.0300005*t_age_nhdl) + (0.0232747*t_age_hdl) + (-0.0927024*t_age_sbp) ///
        + (-0.2018525*t_age_dm) + (-0.0970527*t_age_smk) + (-0.1217081*t_age_egfr) if female == 0

    * base_10yr - total_cvd
    gen double lp_b10_cvd = .
    replace lp_b10_cvd = -3.307728 + (0.7939329*t_age) + (0.0305239*t_nhdl) ///
        + (-0.1606857*t_hdl) + (-0.2394003*t_sbp_lo) + (0.3600781*t_sbp_hi) ///
        + (0.8667604*t_dm) + (0.5360739*t_smk) + (0.6045917*t_egfr_lo) ///
        + (0.0433769*t_egfr_hi) + (0.3151672*t_bptx) + (-0.1477655*t_statin) ///
        + (-0.0663612*t_bptx_sbp) + (0.1197879*t_statin_nhdl) ///
        + (-0.0819715*t_age_nhdl) + (0.0306769*t_age_hdl) + (-0.0946348*t_age_sbp) ///
        + (-0.27057*t_age_dm) + (-0.078715*t_age_smk) + (-0.1637806*t_age_egfr) if female == 1
    replace lp_b10_cvd = -3.031168 + (0.7688528*t_age) + (0.0736174*t_nhdl) ///
        + (-0.0954431*t_hdl) + (-0.4347345*t_sbp_lo) + (0.3362658*t_sbp_hi) ///
        + (0.7692857*t_dm) + (0.4386871*t_smk) + (0.5378979*t_egfr_lo) ///
        + (0.0164827*t_egfr_hi) + (0.288879*t_bptx) + (-0.1337349*t_statin) ///
        + (-0.0475924*t_bptx_sbp) + (0.150273*t_statin_nhdl) ///
        + (-0.0517874*t_age_nhdl) + (0.0191169*t_age_hdl) + (-0.1049477*t_age_sbp) ///
        + (-0.2251948*t_age_dm) + (-0.0895067*t_age_smk) + (-0.1543702*t_age_egfr) if female == 0

    * base_30yr - ascvd
    gen double lp_b30_ascvd = .
    replace lp_b30_ascvd = -1.974074 + (0.4669202*t_age) + (-0.0893118*t_age2) ///
        + (0.1256901*t_nhdl) + (-0.1542255*t_hdl) + (-0.0018093*t_sbp_lo) ///
        + (0.322949*t_sbp_hi) + (0.6296707*t_dm) + (0.268292*t_smk) ///
        + (0.100106*t_egfr_lo) + (0.0499663*t_egfr_hi) + (0.1875292*t_bptx) ///
        + (0.0152476*t_statin) + (-0.0276123*t_bptx_sbp) ///
        + (0.0736147*t_statin_nhdl) + (-0.0521962*t_age_nhdl) ///
        + (0.0316918*t_age_hdl) + (-0.1046101*t_age_sbp) + (-0.2727793*t_age_dm) ///
        + (-0.1530907*t_age_smk) + (-0.1299149*t_age_egfr) if female == 1
    replace lp_b30_ascvd = -1.736444 + (0.3994099*t_age) + (-0.0937484*t_age2) ///
        + (0.1744643*t_nhdl) + (-0.120203*t_hdl) + (-0.0665117*t_sbp_lo) ///
        + (0.2753037*t_sbp_hi) + (0.4790257*t_dm) + (0.1782635*t_smk) ///
        + (-0.0218789*t_egfr_lo) + (0.0602553*t_egfr_hi) + (0.1421182*t_bptx) ///
        + (0.0135996*t_statin) + (-0.0218265*t_bptx_sbp) ///
        + (0.1013148*t_statin_nhdl) + (-0.0312619*t_age_nhdl) ///
        + (0.020673*t_age_hdl) + (-0.0920935*t_age_sbp) + (-0.2159947*t_age_dm) ///
        + (-0.1548811*t_age_smk) + (-0.0712547*t_age_egfr) if female == 0

    * full_10yr - ascvd
    gen double lp_f10_ascvd = .
    replace lp_f10_ascvd = -4.291503 + (0.7023067*t_age) + (0.0898765*t_nhdl) ///
        + (-0.1407316*t_hdl) + (-0.0256648*t_sbp_lo) + (0.314511*t_sbp_hi) ///
        + (0.4799217*t_dm) + (0.4062049*t_smk) + (0.3847744*t_egfr_lo) ///
        + (0.0495174*t_egfr_hi) + (0.2133861*t_bptx) + (-0.0678552*t_statin) ///
        + (-0.0451416*t_bptx_sbp) + (0.0788187*t_statin_nhdl) ///
        + (-0.0535985*t_age_nhdl) + (0.0291762*t_age_hdl) + (-0.0961839*t_age_sbp) ///
        + (-0.2001466*t_age_dm) + (-0.0586472*t_age_smk) + (-0.1537791*t_age_egfr) ///
        + (0.1413965*t_sdi46) + (0.228136*t_sdi710) + (0.1588908*t_sdimiss) ///
        + (0.1371824*t_lnuacr) + (0.0061613*t_uacrmiss) + (0.123192*t_a1c_dm) ///
        + (0.1410572*t_a1c_nodm) + (0.005866*t_a1cmiss) if female == 1
    replace lp_f10_ascvd = -3.969788 + (0.7128741*t_age) + (0.1465201*t_nhdl) ///
        + (-0.1125794*t_hdl) + (-0.3387216*t_sbp_lo) + (0.2980252*t_sbp_hi) ///
        + (0.399583*t_dm) + (0.3379111*t_smk) + (0.2582604*t_egfr_lo) ///
        + (0.0147769*t_egfr_hi) + (0.1686621*t_bptx) + (-0.1073619*t_statin) ///
        + (-0.0381038*t_bptx_sbp) + (0.1034169*t_statin_nhdl) ///
        + (-0.0228755*t_age_nhdl) + (0.0267453*t_age_hdl) + (-0.0897449*t_age_sbp) ///
        + (-0.1497464*t_age_dm) + (-0.077206*t_age_smk) + (-0.1198368*t_age_egfr) ///
        + (0.0651121*t_sdi46) + (0.2676683*t_sdi710) + (0.1388492*t_sdimiss) ///
        + (0.1375837*t_lnuacr) + (0.0652944*t_uacrmiss) + (0.101282*t_a1c_dm) ///
        + (0.1092726*t_a1c_nodm) + (-0.0112852*t_a1cmiss) if female == 0

    * --- risks in % ---
    gen double r10_ascvd      = 100 * invlogit(lp_b10_ascvd)
    gen double r10_cvd        = 100 * invlogit(lp_b10_cvd)
    gen double r30_ascvd      = 100 * invlogit(lp_b30_ascvd)
    gen double r10_ascvd_full = 100 * invlogit(lp_f10_ascvd)
    drop t_* lp_*
end

capture program drop prevent_hf
program define prevent_hf
    syntax, GENerate(name)
    tempvar age slo shi blo bhi elo ehi lp
    quietly {
        gen double `age' = (h_age - 55) / 10
        gen double `slo' = (min(h_sbp, 110) - 110) / 20
        gen double `shi' = (max(h_sbp, 110) - 130) / 20
        gen double `blo' = (min(h_bmi, 30) - 25) / 5
        gen double `bhi' = (max(h_bmi, 30) - 30) / 5
        gen double `elo' = (min(h_egfr, 60) - 60) / -15
        gen double `ehi' = (max(h_egfr, 60) - 90) / -15
        gen double `lp' = .
        replace `lp' = -4.310409 + 0.8998235*`age' - 0.4559771*`slo' + 0.3576505*`shi' ///
            + 1.038346*h_dm + 0.583916*h_smk - 0.0072294*`blo' + 0.2997706*`bhi' ///
            + 0.7451638*`elo' + 0.0557087*`ehi' + 0.3534442*h_bptx ///
            - 0.0981511*h_bptx*`shi' - 0.0946663*`age'*`shi' - 0.3581041*`age'*h_dm ///
            - 0.1159453*`age'*h_smk - 0.003878*`age'*`bhi' - 0.1884289*`age'*`elo' ///
            if female == 1
        replace `lp' = -3.946391 + 0.8972642*`age' - 0.6811466*`slo' + 0.3634461*`shi' ///
            + 0.923776*h_dm + 0.5023736*h_smk - 0.0485841*`blo' + 0.3726929*`bhi' ///
            + 0.6926917*`elo' + 0.0251827*`ehi' + 0.2980922*h_bptx ///
            - 0.0497731*h_bptx*`shi' - 0.1289201*`age'*`shi' - 0.3040924*`age'*h_dm ///
            - 0.1401688*`age'*h_smk + 0.0068126*`age'*`bhi' - 0.1797778*`age'*`elo' ///
            if female == 0
        gen double `generate' = 100 * invlogit(`lp')
    }
end
