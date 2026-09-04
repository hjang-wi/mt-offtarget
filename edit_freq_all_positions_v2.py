import pandas as pd
import glob
import ast
import re
import os

# ======================
# Paths
# ======================
# mito_offtarget_single_sample_v2.py writes each sample's REDItools output to:
#   {sample}_output_revised/reditools_MT_revised/denovo_<id>/outTable_<id>
# Point DATA_DIR at the parent folder containing all the *_output_revised dirs
# (e.g. the DdCBE folder), not at a flat outTable directory.
# NOTE: replace with your own path before running (same fastq_dir/output
# parent used in mito_offtarget_single_sample_v2.py).
DATA_DIR = "<YOUR_DATA_DIR>"                     # e.g. "/path/to/mtDNA_WGS/DdCBE"
OUTPUT_DIR = os.path.join(DATA_DIR, "processed")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================
# Function to compute C->T or G->A edit %
# ======================
def compute_edit_percent(ref, A, C, G, T):
    ref = ref.upper()
    if ref == "C":
        total = A + C + G+ T
        return (T / total * 100)
    elif ref == "G":
        total = A + C + G+ T
        return (A / total * 100)
    else:
        return 0

def get_sample_id(outtable_path):
    """Derive the sample_id from the enclosing folder structure, e.g.:
        ND1-DdCBE-TFL-rep1_S7_L001_output_revised/reditools_MT_revised/denovo_123/outTable_123
        -> "ND1-DdCBE-TFL-rep1_S7_L001"
    (outTable's own filename is just REDItools' internal job ID and carries
    no sample information - the sample_id lives in the folder name.)"""
    output_dir = os.path.dirname(os.path.dirname(os.path.dirname(outtable_path)))
    return os.path.basename(output_dir).replace("_output_revised", "")

# ======================
# Process all outTable files
# ======================
file_list = glob.glob(os.path.join(
    DATA_DIR, "*_output_revised", "reditools_MT_revised", "denovo_*", "outTable_*"
))

for filepath in file_list:
    sample_id = get_sample_id(filepath)
    print(f"Processing {filepath} (sample: {sample_id})...")
    processed_rows = []

    with open(filepath, "r", encoding="utf-8") as f:
        header = f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue

            # Split first 6 columns; rest contains BaseCount etc
            parts = line.strip().split(None, 6)
            if len(parts) < 7:
                continue
            region, pos, ref, strand, coverage, meanq, rest = parts

            # Only process positions with Reference C or G
            if ref.upper() not in ["C", "G"]:
                continue

            # Extract BaseCount[A,C,G,T] string
            match = re.search(r"\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]", rest)
            if match:
                basecount_str = "[" + ",".join(match.groups()) + "]"
                A = int(match.group(1))
                C = int(match.group(2))
                G = int(match.group(3))
                T = int(match.group(4))
            else:
                basecount_str = "[0,0,0,0]"
                A = C = G = T = 0

            edit_percent = compute_edit_percent(ref, A, C, G, T)

            processed_rows.append([
                region, pos, ref, strand, coverage, meanq,
                basecount_str, A, C, G, T, edit_percent
            ])

    if processed_rows:
        df_out = pd.DataFrame(processed_rows, columns=[
            "Region","Position","Reference","Strand","Coverage-q30","MeanQ",
            "BaseCount[A,C,G,T]", "A_count","C_count","G_count","T_count","edit_percent"
        ])
        outname = f"outTable_{sample_id}_processed.csv"
        df_out.to_csv(os.path.join(OUTPUT_DIR, outname), index=False)
        print(f"Saved {outname}")
    else:
        print(f"No C/G positions found in {filepath} (sample: {sample_id})")

print("All files processed successfully!")
