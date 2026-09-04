# mtDNA Off-Target Editing Analysis Pipeline

Pipeline for computing mtDNA-wide C→T / G→A off-target editing
frequency from sequencing data. It excludes on-target sites and background
germline SNVs, then generates volcano-style visualizations.

## Pipeline Overview

```
samples.txt
    │
    ▼
run_mito_samplev2.sh  (SLURM array job)
    │  runs the script below for each sample_id
    ▼
mito_offtarget_single_sample_v2.py   (Python 2.7)
    - BWA-MEM alignment (paired-end, MT reference)
    - samtools sort / fixmate / index
    - Runs REDItoolDenovo.py (via node-local /tmp scratch)
    - Only produces outTable_* (no parsing/filtering here)
    │
    ▼
edit_freq_all_positions_v2.py        (Python 3)
    - Batch-reads every sample's outTable_*
    - Parses BaseCount[A,C,G,T], keeps only C/G reference positions,
      computes edit_percent, saves outTable_{sample_id}_processed.csv
    │
    ▼
reditool-outTable-processed-analysis_final.py   (Python 3)
    - Loads all processed CSVs
    - Excludes germline SNVs (positions ≥5% in BOTH a treated and untreated sample)
    - Excludes on-target positions (ND1: 3654/3655/3657/3659, ND4: 11696/11698)
    - Computes mtDNA-wide off-target editing frequency (sum of edited reads / sum of coverage)
    - Saves combined_editing_summary1.csv
    - Groups replicates by condition (target/group/dose) and generates mean ± SD volcano plots (SVG/PDF)
```

## Requirements

This pipeline uses two separate Python environments.

| Step | Script | Python | Key dependencies |
|---|---|---|---|
| 1 | `mito_offtarget_single_sample_v2.py` | 2.7 | `bwa`, `samtools`, `REDItools 1.2.1` (`REDItoolDenovo.py`) |
| 2 | `edit_freq_all_positions_v2.py` | 3.x | `pandas` |
| 3 | `reditool-outTable-processed-analysis_final.py` | 3.x | `pandas`, `matplotlib` |

```bash
# Step 1 environment (conda)
conda create -n reditools_py2 python=2.7
conda activate reditools_py2
conda install -c bioconda bwa samtools

# Step 2 & 3 environment
pip install pandas matplotlib
```

REDItools 1.2.1 must be installed/cloned separately:
https://github.com/BioinfoUNIBA/REDItools

## Directory Structure (generated at runtime)

```
{DATA_DIR}/
├── {sample}_output_revised/
│   ├── aligned_MT.sorted.bam(.bai)
│   └── reditools_MT_revised/
│       └── denovo_<id>/
│           └── outTable_<id>
├── processed/
│   └── outTable_{sample_id}_processed.csv
└── plots/
    ├── combined_editing_summary1.csv
    └── volcano_{condition}.svg / .pdf
```

## Configuration

This repo ships with lab-specific paths and personal info removed and
replaced with placeholders. **You must edit these before running anything.**

| File | Variable | Placeholder | Example value |
|---|---|---|---|
| `mito_offtarget_single_sample_v2.py` | `fastq_dir` | `<YOUR_DATA_DIR>` | `/path/to/mtDNA_WGS/DdCBE` |
| `mito_offtarget_single_sample_v2.py` | `mt_reference` | `<YOUR_REDITOOLS_DIR>/Homo_sapiens.GRCh38_MT.fa` | `/path/to/REDItools-1.2.1/Homo_sapiens.GRCh38_MT.fa` |
| `mito_offtarget_single_sample_v2.py` | `reditool` | `<YOUR_REDITOOLS_DIR>/reditools/REDItoolDenovo.py` | `/path/to/REDItools-1.2.1/reditools/REDItoolDenovo.py` |
| `edit_freq_all_positions_v2.py` | `DATA_DIR` | `<YOUR_DATA_DIR>` | same folder as `fastq_dir` above |
| `run_mito_samplev2.sh` | `--mail-user` | `<YOUR_EMAIL>` | your own email, for SLURM job notifications |

`reditool-outTable-processed-analysis_final.py` needs no changes — its
`DATA_DIR = "./processed"` is already a relative path.

You'll also likely want to adjust the `#SBATCH` options in
`run_mito_samplev2.sh` (partition name, `--array` indices, memory/time
limits) to match your own cluster.

## Usage

### 1. Prepare the sample list
List one sample_id per line in `samples.txt` — this is the portion of the FASTQ filename before `_R1_`/`_R2_`. 

### 2. Set your paths
Edit the hardcoded paths at the top of each script to match your
environment:
- `mito_offtarget_single_sample_v2.py`: `fastq_dir`, `mt_reference`, `reditool`
- `edit_freq_all_positions_v2.py`: `DATA_DIR`
- `reditool-outTable-processed-analysis_final.py`: `DATA_DIR`

See **Configuration** above for the exact placeholder names to replace.

### 3. Run alignment + REDItools on the SLURM cluster
```bash
mkdir -p logs
sbatch run_mito_samplev2.sh
```
The `--array` indices correspond to line numbers in `samples.txt`. Update the
`#SBATCH` options (partition, account, email, etc.) to match your cluster.

### 4. Batch-parse all sample outTables
```bash
conda activate <python3 env>
python edit_freq_all_positions_v2.py
```

### 5. Run off-target analysis and generate plots
```bash
python reditool-outTable-processed-analysis_final.py
```

## Output

- **`combined_editing_summary1.csv`**: per-sample off-target edited-site
  counts, mtDNA-wide off-target editing frequency (%), and ND1/ND4 on-target
  frequency (%)
- **`volcano_{condition}.svg / .pdf`**: mtDNA position (x-axis) vs. editing
  frequency (y-axis), replicate mean ± SD, on-target (vermillion) vs.
  off-target CT/GA (blue)


