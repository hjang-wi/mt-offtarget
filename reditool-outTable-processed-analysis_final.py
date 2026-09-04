# -*- coding: utf-8 -*-

import re
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt

# Ensure text stays as editable text in vector output
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['axes.linewidth'] = 0.75
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['ytick.major.width'] = 0.75

# colorblind-safe palette, avoiding red-green pairings
COLOR_BACKGROUND = "#999999"   # gray
COLOR_OFFTARGET  = "#0072B2"   # blue
COLOR_ONTARGET   = "#D55E00"   # vermillion

# Figure size 
FIG_WIDTH_MM = 180
FIG_HEIGHT_MM = 70
MM_TO_IN = 1 / 25.4
FIG_SIZE_IN = (FIG_WIDTH_MM * MM_TO_IN, FIG_HEIGHT_MM * MM_TO_IN)

MARKER_SIZE = 6  
ERRORBAR_CAPSIZE = 2
ERRORBAR_LINEWIDTH = 0.75

PLOTS_DIR = "plots"
os.makedirs(PLOTS_DIR, exist_ok=True)

# =========================================
# User-defined settings
# =========================================

# Vector output format(s) to save alongside/instead of PNG.
VECTOR_FORMATS = ["svg", "pdf"]
SAVE_PNG_TOO = False  # PNG preview disabled — only vector (SVG/PDF) output

DATA_DIR = "./processed"  # folder containing outTable files

# On-target sites for each DdCBE
ONTARGET_ND1 = [3654, 3655, 3657, 3659]
ONTARGET_ND4 = [11696, 11698]

MITO_GENOME_SIZE = 16569  # human mtDNA length

# =========================================
# Helper functions
# =========================================

def parse_sample_id(sample_id):

    m = re.match(
        r'^(?:(?P<target>ND1|ND4)-)?'
        r'(?P<group_token>DdCBE-TFL|DdCBE-TFH|VLP|UT)'
        r'-rep(?P<rep>\d+)_S(?P<snum>\d+)(?:_L\d+)?$',
        sample_id
    )
    if not m:
        return None
    target = m.group('target')  # None for UT
    token = m.group('group_token')
    rep = int(m.group('rep'))
    if token == 'UT':
        group, dose = 'Control', None
    elif token.startswith('DdCBE'):
        group, dose = 'DdCBE', token.split('-')[1]  # TFL or TFH
    else:
        group, dose = token, None  # VLP
    return {'target': target, 'group': group, 'dose': dose, 'rep': rep}

def get_sample_id_from_processed_filename(fname):
    """outTable_{sample_id}_processed.csv -> sample_id, as written by the
    updated edit_freq_all_positions.py."""
    key = re.sub(r'\.csv$', '', fname, flags=re.IGNORECASE)
    key = re.sub(r'^outTable_', '', key, flags=re.IGNORECASE)
    key = re.sub(r'_processed$', '', key, flags=re.IGNORECASE)
    return key

# Manual overrides - only needed for sample_ids that don't match the
SAMPLE_GROUP_OVERRIDE = {}
SAMPLE_TARGET_REGION_OVERRIDE = {}
REPLICATE_GROUP_OVERRIDE = {}

def detect_group(fname):
    sample_id = get_sample_id_from_processed_filename(fname)
    parsed = parse_sample_id(sample_id)
    return SAMPLE_GROUP_OVERRIDE.get(sample_id, parsed['group'] if parsed else 'Unknown')

def detect_target_region(fname):
    sample_id = get_sample_id_from_processed_filename(fname)
    parsed = parse_sample_id(sample_id)
    return SAMPLE_TARGET_REGION_OVERRIDE.get(sample_id, parsed['target'] if parsed else None)
    
    
def read_file(path):
    """Read Excel/CSV files with 'outTable_' prefix and convert edit_percent to fraction."""
    df = pd.read_csv(path, encoding='utf-8')
    df["Frequency"] = df["edit_percent"] / 100  # convert percent to fraction
    return df

def get_condition_key(fname):

    sample_id = get_sample_id_from_processed_filename(fname)
    if sample_id in REPLICATE_GROUP_OVERRIDE:
        return REPLICATE_GROUP_OVERRIDE[sample_id]
    parsed = parse_sample_id(sample_id)
    if parsed is None:
        return sample_id  
    if parsed['group'] == 'Control':
        return 'Control'
    parts = [parsed['target'] or 'none', parsed['group']]
    if parsed['dose']:
        parts.append(parsed['dose'])
    return "_".join(parts)

# =========================================
# Step 1 - Load all sample files
# =========================================

all_samples = {}

for path in glob.glob(os.path.join(DATA_DIR, "outTable_*.csv")):
    fname = os.path.basename(path)
    df = read_file(path)
    df["sample"] = fname
    df["group"] = detect_group(fname)
    df["target_region"] = detect_target_region(fname)
    all_samples[fname] = df

print(f"Loaded {len(all_samples)} samples.")

# =========================================
# Step 2 - Determine germline SNV sites
#          (>=5% conversion in BOTH a treated AND an untreated sample)
# =========================================

SNV_THRESHOLD = 0.05  # 5%

treated_dfs = [df for df in all_samples.values() if df["group"].iloc[0] != "Control"]
untreated_dfs = [df for df in all_samples.values() if df["group"].iloc[0] == "Control"]

if len(untreated_dfs) == 0:
    raise ValueError("No untreated controls found!")
if len(treated_dfs) == 0:
    raise ValueError("No treated samples found!")

treated_concat = pd.concat(treated_dfs, ignore_index=True)
untreated_concat = pd.concat(untreated_dfs, ignore_index=True)

positions_ge5_treated = set(
    treated_concat.loc[treated_concat["Frequency"] >= SNV_THRESHOLD, "Position"]
)
positions_ge5_untreated = set(
    untreated_concat.loc[untreated_concat["Frequency"] >= SNV_THRESHOLD, "Position"]
)

# A position is a germline SNV only if it reaches >=5% in a treated sample
# AND in an untreated sample - present regardless of treatment.
background_sites = sorted(positions_ge5_treated & positions_ge5_untreated)

print(f"Positions >=5% in >=1 treated sample:   {len(positions_ge5_treated)}")
print(f"Positions >=5% in >=1 untreated sample: {len(positions_ge5_untreated)}")
print(f"SNV sites excluded (>=5% in both):      {len(background_sites)}")


# =========================================
# Step 3 - mito off-target summary
# =========================================

OFFTARGET_COUNT_THRESHOLD = 0.001  # 0.1% - only sites above this count as "edited"

# Representative on-target position used for the ND1/ND4 on-target frequency columns
ND1_ONTARGET_POS = 3654
ND4_ONTARGET_POS = 11698

def get_position_freq_pct(df, position):
    """Return the Frequency (as a percent) at a given Position in df, or
    None if that position isn't present in this sample's data."""
    match = df.loc[df["Position"] == position, "Frequency"]
    if match.empty:
        return None
    return match.iloc[0] * 100

combined_summary = []

for fname, df in all_samples.items():
    target = df["target_region"].iloc[0]
    group = df["group"].iloc[0]

    # On-target frequency at the fixed ND1/ND4 representative positions,
    nd1_ontarget_freq = get_position_freq_pct(df, ND1_ONTARGET_POS)
    nd4_ontarget_freq = get_position_freq_pct(df, ND4_ONTARGET_POS)

    # Determine on-target sites
    if target == "ND1":
        on_targets = ONTARGET_ND1
    elif target == "ND4":
        on_targets = ONTARGET_ND4
    else:
        on_targets = []

    df_offtarget = df.copy()

    # Exclude on-target positions
    df_offtarget = df_offtarget[~df_offtarget["Position"].isin(on_targets)]

    # Exclude germline SNV positions (>=5% in both treated & untreated)
    df_offtarget = df_offtarget[~df_offtarget["Position"].isin(background_sites)]

    # --- n_offtarget_edited_sites (editing frequency > 0.1% count as "edited" sites)
    edited_df = df_offtarget[df_offtarget["Frequency"] > OFFTARGET_COUNT_THRESHOLD]
    n_edited_sites = edited_df.shape[0]

    # --- mito_offtarget_frequency(%) ---
    # Sum of edited reads (C->T or G->A) over ALL off-target positions,
    # divided by the sum of total coverage over ALL off-target positions.
    df_offtarget = df_offtarget.copy()
    df_offtarget["edited_reads"] = df_offtarget.apply(
        lambda row: row["T_count"] if row["Reference"] == "C"
                    else (row["A_count"] if row["Reference"] == "G" else 0),
        axis=1
    )
    df_offtarget["total_coverage"] = df_offtarget["Coverage-q30"]

    total_edits = df_offtarget["edited_reads"].sum()
    total_cov = df_offtarget["total_coverage"].sum()
    mito_offtarget_freq = (total_edits / total_cov * 100) if total_cov > 0 else 0

    combined_summary.append({
        "sample": fname,
        "group": group,
        "target_region": target,
        "n_offtarget_edited_sites": n_edited_sites,
        "mito_offtarget_frequency(%)": mito_offtarget_freq,
        "ND1_ontarget_freq(%)": nd1_ontarget_freq,
        "ND4_ontarget_freq(%)": nd4_ontarget_freq
    })

combined_df = pd.DataFrame(combined_summary)
combined_df.to_csv(os.path.join(PLOTS_DIR, "combined_editing_summary1.csv"), index=False, encoding='utf-8')
print("\nCombined genome-wide and mito off-target summary1:")
print(combined_df)



# =========================================
# Step 4 - Volcano-style plots (replicates grouped by condition, mean +/- SD per position)
# =========================================

EDITING_FREQ_THRESHOLD = 0.01  # 1%

# Group loaded samples by condition 
from collections import defaultdict
condition_groups = defaultdict(list)
for fname, df in all_samples.items():
    condition_groups[get_condition_key(fname)].append(df)

print(f"\nGrouped {len(all_samples)} sample files into {len(condition_groups)} conditions:")
for key, dfs in condition_groups.items():
    print(f"  {key}: {len(dfs)} replicate(s)")

for condition_key, df_list in condition_groups.items():

    group = df_list[0]["group"].iloc[0]
    target = df_list[0]["target_region"].iloc[0]
    n_reps = len(df_list)

    # Combine all replicates and compute mean +/- SD per position.
    combined = pd.concat(df_list, ignore_index=True)
    agg = (
        combined.groupby("Position")["Frequency"]
        .agg(["mean", "std"])
        .reset_index()
        .rename(columns={"mean": "Frequency"})
    )
    agg["std"] = agg["std"].fillna(0) 
    # Reference base is the same across replicates at a given position
    ref_map = combined.drop_duplicates("Position").set_index("Position")["Reference"]
    agg["Reference"] = agg["Position"].map(ref_map)
    df = agg

    plt.figure(figsize=FIG_SIZE_IN)

    # Determine on-target sites
    if target == "ND1":
        on_targets = ONTARGET_ND1
    elif target == "ND4":
        on_targets = ONTARGET_ND4
    else:
        on_targets = []

    # Identify off-target CT/GA positions (exclude background & on-target)
    df_filtered = df[~df["Position"].isin(background_sites)]
    df_filtered = df_filtered.copy()
    df_filtered["AllSubs"] = df_filtered.apply(lambda row: "CT" if row["Reference"]=="C" else ("GA" if row["Reference"]=="G" else ""), axis=1)
    off_targets = df_filtered[df_filtered["AllSubs"].isin(["CT", "GA"])]["Position"].tolist()
    off_targets = [pos for pos in off_targets if pos not in on_targets]



    # Blue: off-target CT/GA above threshold
    if off_targets:
        ot_cg_df = df[df["Position"].isin(off_targets) & (df["Frequency"] > EDITING_FREQ_THRESHOLD)]
        if not ot_cg_df.empty:
            plt.errorbar(
                ot_cg_df["Position"], ot_cg_df["Frequency"],
                yerr=ot_cg_df["std"], fmt="o",
                color=COLOR_OFFTARGET, ecolor=COLOR_OFFTARGET,
                markersize=MARKER_SIZE ** 0.5, linewidth=0,
                elinewidth=ERRORBAR_LINEWIDTH, capsize=ERRORBAR_CAPSIZE,
                label="Off-target CT/GA"
            )

    # Vermillion: on-target above threshold
    if on_targets:
        ot_df = df[df["Position"].isin(on_targets) & (df["Frequency"] > EDITING_FREQ_THRESHOLD)]
        if not ot_df.empty:
            plt.errorbar(
                ot_df["Position"], ot_df["Frequency"],
                yerr=ot_df["std"], fmt="o",
                color=COLOR_ONTARGET, ecolor=COLOR_ONTARGET,
                markersize=MARKER_SIZE ** 0.5, linewidth=0,
                elinewidth=ERRORBAR_LINEWIDTH, capsize=ERRORBAR_CAPSIZE,
                label="On-target"
            )

    plt.xlabel("Position")
    plt.ylabel("Editing Frequency")
    # No in-figure title: Nature-style figures carry the title/description
    # in the figure caption, not inside the plot itself.
    plt.xlim(0, MITO_GENOME_SIZE)
    plt.xticks(
        [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000],
        ["0", "2000", "4000", "6000", "8000", "10000", "12000", "14000", "16000"]
    )
    plt.ylim(0, 1.0)
    plt.yticks(
        [0, 0.2, 0.4, 0.6, 0.8, 1.0],
        ["0%", "20%", "40%", "60%", "80%", "100%"]
    )
    # Build custom legend handles with plain dot markers 
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="",
               color=COLOR_OFFTARGET,
               markersize=MARKER_SIZE ** 0.5, label="Off-target CT/GA"),
        Line2D([0], [0], marker="o", linestyle="",
               color=COLOR_ONTARGET,
               markersize=MARKER_SIZE ** 0.5, label="On-target"),
    ]
    plt.legend(handles=legend_handles)
    plt.tight_layout()

    base_name = os.path.join(PLOTS_DIR, f"volcano_{condition_key}")

    # Save vector formats
    for ext in VECTOR_FORMATS:
        plt.savefig(f"{base_name}.{ext}")

    # Optionally still keep a PNG for quick previewing / slides
    if SAVE_PNG_TOO:
        plt.savefig(f"{base_name}.png", dpi=200)

    plt.close()

print("\nVolcano plots (mean +/- SD across replicates, shown as dots with error bars) "
      "with blue CT/GA off-target and vermillion on-target points "
      "saved successfully.")
