"""
00_download_nhanes.py
Downloads the public NHANES files used by this analysis from the CDC website into
data/nhanes_2021_2023 and data/nhanes_2017_2020. Files already present are skipped.
"""
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
FILES = {
    ("2021", "nhanes_2021_2023"): ["DEMO_L", "TRIGLY_L", "TCHOL_L", "HDL_L", "BPQ_L", "BMX_L", "DIQ_L", "GHB_L",
                                   "INS_L", "GLU_L", "LUX_L", "BPXO_L", "BIOPRO_L", "ALB_CR_L", "SMQ_L", "MCQ_L"],
    ("2017", "nhanes_2017_2020"): ["P_DEMO", "P_TRIGLY", "P_TCHOL", "P_HDL", "P_BPXO", "P_BPQ", "P_DIQ", "P_GHB",
                                   "P_BIOPRO", "P_SMQ", "P_MCQ", "P_BMX", "P_RXQ_RX"],
}
for (year, folder), names in FILES.items():
    dest = ROOT / "data" / folder
    dest.mkdir(parents=True, exist_ok=True)
    for n in names:
        f = dest / f"{n}.xpt"
        if f.exists() and f.stat().st_size > 0:
            print("present ", f.name)
            continue
        url = f"{BASE}/{year}/DataFiles/{n}.xpt"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research data download)"})
        data = urllib.request.urlopen(req, timeout=300).read()
        if not data.startswith(b"HEADER RECORD"):
            raise SystemExit(f"{url} did not return an XPT file")
        f.write_bytes(data)
        print(f"saved    {f.name} ({len(data) / 1e6:.1f} MB)")
