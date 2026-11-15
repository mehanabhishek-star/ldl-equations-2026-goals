"""
run_all.py - runs the whole analysis from the repository root.
    python run_all.py --stata "C:/Program Files/Stata17/StataMP-64.exe"
    python run_all.py --stata /usr/local/stata17/stata-mp          (macOS/Linux)
Add --download to fetch the NHANES files first (code/00_download_nhanes.py).
Stops at the first Stata error and names the log to read.
"""
import argparse, os, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = ["01_import_2021_2023.do", "02_prevent_scores.do", "03_goals_and_analysis_data.do",
         "04_interval_methods.do", "05_main_table.do", "06_additional_estimates.do",
         "07_mh_lookup_build.py", "08_participant_level_and_assay.do", "09_pooled_2017_2023.do",
         "10_figure.py", "11_build_tables.py"]

def run_stata(stata, step):
    do = ROOT / "code" / step
    cmd = [stata, "/e", "do", str(do)] if os.name == "nt" else [stata, "-b", "do", str(do)]
    subprocess.run(cmd, cwd=ROOT, check=True)
    # Stata writes the batch log next to the repository root; keep it in logs/ (not part of the repository)
    (ROOT / "logs").mkdir(exist_ok=True)
    src, log = ROOT / (Path(step).stem + ".log"), ROOT / "logs" / (Path(step).stem + ".log")
    if src.exists():
        src.replace(log)
    text = log.read_text(errors="replace") if log.exists() else ""
    err = re.search(r"^r\((\d+)\);", text, flags=re.M)
    if err or not text:
        sys.exit(f"Stata stopped in {step} (r({err.group(1) if err else '?'})); see {log}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stata", required=True, help="path to the Stata executable")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--from-step", default="01")
    a = ap.parse_args()
    if a.download:
        subprocess.run([sys.executable, str(ROOT / "code" / "00_download_nhanes.py")], check=True)
    for step in STEPS:
        if step[:2] < a.from_step:
            continue
        print("==>", step, flush=True)
        if step.endswith(".do"):
            run_stata(a.stata, step)
        else:
            subprocess.run([sys.executable, str(ROOT / "code" / step)], cwd=ROOT, check=True)
    print("Done. Results in results/, figure in figures/.")

if __name__ == "__main__":
    main()
