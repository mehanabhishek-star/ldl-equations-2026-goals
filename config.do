*=======================================================================*
* config.do - folder locations. Run Stata from the repository root.
* run_all.py calls this automatically; each analysis do-file also calls it
* if the paths are not yet set.
*=======================================================================*
global ROOT     "`c(pwd)'"
global CODE     "$ROOT/code"
global RESULTS  "$ROOT/results"
global FIGURES  "$ROOT/figures"
global DATA2123 "$ROOT/data/nhanes_2021_2023"
global DATA1720 "$ROOT/data/nhanes_2017_2020"
capture mkdir "$RESULTS"
capture mkdir "$FIGURES"
