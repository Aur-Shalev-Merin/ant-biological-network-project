"""
Pull a small set of files from the Stroeymeyt et al. 2018 Zenodo archive
(data.zip, 24.7 GB) without downloading the whole thing.

remotezip reads the zip's directory and individual files over HTTP range
requests, so you only download what you ask for.

Usage:
    pip install remotezip
    python fetch_ant_files.py           # download the selected files
    python fetch_ant_files.py --list    # just save every path in the zip to file_list.txt
"""

import os
import sys
import shutil
from itertools import islice

from remotezip import RemoteZip

# If this URL fails, copy the download link for data.zip from
# https://zenodo.org/records/1480665 and paste it here.
URL = "https://zenodo.org/records/1480665/files/data.zip?download=1"

OUT_DIR = "ant_data"
TRACKING_LINES = 7200  # 2 frames/sec, so 7200 lines is about 1 hour

# Matched by suffix, so it doesn't matter what top-level folder the zip uses.
WANTED = [
    # 1. colony info
    "main_experiment/original_data/info.txt",
    # 2. task groups and treated workers
    "main_experiment/original_data/task_groups.txt",
    "main_experiment/original_data/treated_worker_list.txt",
    # 3. their observed vs random networks (novelty check)
    "topology_comparison/colony276_pathogen_PreTreatment_TH-27_TD9_interactions_observed.txt",
    "topology_comparison/colony276_pathogen_PreTreatment_TH-27_TD9_interactions_random.txt",
    "topology_comparison/colony276_pathogen_PostTreatment_TH21_TD9_interactions_observed.txt",
    # 4. a full post-treatment interaction list for colony 276
    "node_color_f_property/colony276_pathogen_PostTreatment_interactions_observed.txt",
    # 5. body geometry for colony 276
    "main_experiment/original_data/tag_files/colony276_pathogen.tags",
    # optional extras
    "main_experiment/original_data/qPCR/qPCR_results.txt",
    "main_experiment/original_data/plume_files/colony276_pathogen.plume",
]

# This one is ~330 MB, so we only stream the first TRACKING_LINES lines.
TRACKING = "main_experiment/original_data/tracking/colony276_pathogen_PreTreatment.csv"


def find(names, suffix):
    """Return the zip path ending with suffix (warn if 0 or several match)."""
    hits = [n for n in names if n.endswith(suffix)]
    if not hits:
        print(f"  NOT FOUND: {suffix}")
        return None
    if len(hits) > 1:
        print(f"  {len(hits)} matches for {suffix}, using the first:")
        for h in hits:
            print(f"    {h}")
    return hits[0]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Opening remote zip (reads only the file directory)...")

    with RemoteZip(URL) as z:
        names = z.namelist()
        print(f"Zip contains {len(names)} entries.\n")

        if "--list" in sys.argv:
            with open("file_list.txt", "w") as f:
                for info in z.infolist():
                    f.write(f"{info.file_size:>14,}  {info.filename}\n")
            print("Saved full listing to file_list.txt")
            return

        # Small files: download whole
        for suffix in WANTED:
            name = find(names, suffix)
            if name is None:
                continue
            dest = os.path.join(OUT_DIR, os.path.basename(name))
            with z.open(name) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"  saved {dest}")

        # Big tracking file: stream only the first N lines
        name = find(names, TRACKING)
        if name is not None:
            base = os.path.basename(name).replace(".csv", f"_first{TRACKING_LINES}.csv")
            dest = os.path.join(OUT_DIR, base)
            with z.open(name) as src, open(dest, "wb") as dst:
                for line in islice(src, TRACKING_LINES):
                    dst.write(line)
            print(f"  saved {dest} (first {TRACKING_LINES} lines)")

    print(f"\nDone. Files are in ./{OUT_DIR}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nError: {e}")
        print("If this is a range-request or 4xx error, try copying the exact "
              "data.zip download link from the Zenodo page into URL.")
        raise