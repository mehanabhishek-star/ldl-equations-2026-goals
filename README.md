# LDL-C equations and attainment of 2026 guideline goals (NHANES)

Analysis code for:

> Mehan A, Patel D, Sairam K, Ganatra S, Dani SS. **Friedewald vs Martin–Hopkins LDL Cholesterol and Attainment of 2026 Guideline Goals in US Adults.** Manuscript under review.

The code assigns each NHANES participant their own 2026 ACC/AHA LDL-C goal (from PREVENT-ASCVD risk or clinical ASCVD), compares goal attainment by Friedewald, Martin–Hopkins and Sampson LDL-C, and measures how much of the difference the guideline's paired non-HDL-C goal (LDL-C goal + 30 mg/dL) recovers. The main analysis uses NHANES August 2021–August 2023. Sensitivity analyses pool 2017–March 2020 and convert triglycerides between the two NHANES assays.

## Requirements

- **Stata 17** (uses `svy`, `import sasxport5`). No user-written packages are needed.
- **Python 3.10+** with the packages in `requirements.txt`:
  ```
  pip install -r requirements.txt
  ```

## Data

All data are public-use NHANES files from the National Center for Health Statistics (NCHS). They are not stored in this repository. Download them into `data/` with:

```
python code/00_download_nhanes.py
```

| Cycle | Files |
|---|---|
| August 2021–August 2023 (`data/nhanes_2021_2023`) | DEMO_L, TRIGLY_L, TCHOL_L, HDL_L, BPQ_L, BMX_L, DIQ_L, GHB_L, INS_L, GLU_L, LUX_L, BPXO_L, BIOPRO_L, ALB_CR_L, SMQ_L, MCQ_L |
| 2017–March 2020 (`data/nhanes_2017_2020`) | P_DEMO, P_TRIGLY, P_TCHOL, P_HDL, P_BPXO, P_BPQ, P_DIQ, P_GHB, P_BIOPRO, P_SMQ, P_MCQ, P_BMX, P_RXQ_RX |

Users of NHANES data agree to the NCHS data-use restrictions, including making no attempt to identify participants.

## Running the analysis

From the repository root:

```
python run_all.py --stata "C:/Program Files/Stata17/StataMP-64.exe" --download
```

On macOS or Linux, pass the console Stata executable instead, for example `--stata /usr/local/stata17/stata-mp`. Leave out `--download` if the data are already in `data/`, and use `--from-step 05` to restart from a given step. The runner stops at the first Stata error and names the log to read.

To run the steps by hand, start Stata in the repository root and run each do-file in order. Each do-file loads `config.do`, which sets the folder paths. Run the Python steps with `python code/<file>`.

| Step | File | What it does | Main outputs |
|---|---|---|---|
| 00 | `00_download_nhanes.py` | Downloads the NHANES files | `data/` |
| 01 | `01_import_2021_2023.do` | Imports and merges 2021–2023 lipid, laboratory and questionnaire files | `ldl_reclass_data.dta` (data folder) |
| 02 | `02_prevent_scores.do` | PREVENT 10-year ASCVD risk (Khan et al, 2024), checked against the published worked example | `nhanes_prevent.dta` (data folder) |
| 03 | `03_goals_and_analysis_data.do` | Person-specific 2026 goals (Abohashem et al, JAMA 2026); analysis dataset | `results/rev3_analysis_data.dta` |
| 04 | `04_interval_methods.do` | Interval methods for paired differences (t, domain df, Satterthwaite, MOVER) | `results/rev5_p2_ci.csv` |
| 05 | `05_main_table.do` | Main table; margin rule; very-high-risk by Friedewald; previous TG assay (approximate); characteristics; sample flow | `results/rev5_main_table.csv`, `rev5_margin.csv`, `rev5_table1.csv`, `rev5_flow.csv` |
| 06 | `06_additional_estimates.do` | Very-high-risk by Sampson; non-HDL-C added to each equation; exploratory model | `results/rev5_nonhdl_added.csv`, `rev5_vhr_sampson.csv`, `rev5_model.csv` |
| 07 | `07_mh_lookup_build.py` | Recovers the exact Martin–Hopkins VLDL-C values from the NCHS files (see below) | `results/mh_vldl_lookup.csv` |
| 08 | `08_participant_level_and_assay.do` | Participant-level recovery; Sampson − Friedewald; previous TG assay with exact Martin–Hopkins | `results/rev6_person.csv`, `rev6_sf.csv`, `rev6_bridge_rows.csv` |
| 09 | `09_pooled_2017_2023.do` | Pooled 2017–March 2020 + 2021–2023; Satterthwaite intervals; drug names (2017–2020) | `results/rev6_pooled_table.csv`, `rev6_pooled_satt.csv`, `rev6_pooled_person.csv`, `rev6_pooled_rx.csv` |
| 10 | `10_figure.py` | Main figure | `figures/Figure_LDL_equations.tiff/.pdf/.png` |
| 11 | `11_build_tables.py` | Table and eTables 1–5 as a Word file, read from the CSVs | `results/Tables_and_eTables.docx` |

Helper programs: `_kgprop.do` (Korn–Graubard intervals and NCHS presentation flags), `_moverkg.do` (MOVER intervals for paired differences and main-table rows), `_satt.do` (Satterthwaite effective degrees of freedom), `_prevent_programs.do` (PREVENT equations for 2017–March 2020).

The `rev5_`/`rev6_` prefixes on output files mark analysis versions; the files listed above are the ones used in the paper.

## Where each result comes from

| In the paper | File |
|---|---|
| Table; eTable 3 (risk groups) | `results/rev5_main_table.csv`; Sampson − Friedewald in `rev6_sf.csv` |
| eTable 1 (characteristics) | `results/rev5_table1.csv` |
| eTable 2 (sensitivity analyses) | `rev5_main_table.csv`, `rev5_margin.csv`, `rev5_p2_ci.csv`, `rev6_bridge_rows.csv`, `rev6_pooled_table.csv`, `rev6_pooled_satt.csv`, `rev6_pooled_rx.csv` |
| eTable 4 (exploratory model) | `results/rev5_model.csv` |
| eTable 5 (participant level) | `rev6_person.csv`, `rev6_pooled_person.csv`, `rev5_nonhdl_added.csv`, `rev6_pooled_table.csv` |
| Sample flow (eMethods) | `results/rev5_flow.csv` |
| Figure | `figures/Figure_LDL_equations.*` (from `rev5_main_table.csv`) |

## Methods notes

- **Survey design.** All estimates use the NHANES fasting-subsample weights with Taylor-series linearization. Proportions have Korn–Graubard 95% CIs and are checked against NCHS presentation standards. Paired differences (percentage reclassified upward minus downward) have MOVER intervals that combine the Korn–Graubard limits of each direction with their design-based correlation.
- **Exact Martin–Hopkins at a converted TG.** NCHS reports TG, non-HDL-C and Martin–Hopkins LDL-C as integers. So non-HDL-C minus Martin–Hopkins LDL-C (the VLDL-C estimate) is identical for every participant with the same TG within a non-HDL-C stratum (<100, 100–129, 130–159, 160–189, 190–219, ≥220 mg/dL). Step 07 verifies this with no exceptions across both cycles and tabulates the values. Steps 08–09 then compute Martin–Hopkins LDL-C at a converted TG exactly as NCHS would.
- **Triglyceride assays.** The 2017–March 2020 assay was glycerol blanked; the 2021–2023 assay is not. TG were converted with the NCHS equations (TRIGLY_L documentation): forward TG(new) = 12.45 + 1.022 × TG(old); backward TG(old) = −12.19 + 0.9785 × TG(new).
- **Pooling.** Weights were scaled by years covered (3.2/5.2 and 2.0/5.2), following the NCHS caution on combining 2021–2023 with earlier cycles.
- **PREVENT.** Coefficients are from the `preventr` R package (Khan SS et al, Circulation 2024;149:430–449); step 02 stops if the published worked example is not reproduced.

## Citation

Please cite the article once published. Citation details and an archived DOI for this code will be added on acceptance.

## License

MIT (see `LICENSE`). NHANES data are produced by NCHS and are in the public domain.
