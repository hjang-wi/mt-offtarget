#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
MITOCHONDRIAL OFF-TARGET EDITING PIPELINE (Python 2.7)
Performs alignment + REDItools only - stops once the REDItoolDenovo.py
outTable is produced. All parsing/filtering/analysis (BaseCount extraction,
C->T/G->A conversion, SNV removal, on-target exclusion, plotting) is done
downstream, once per batch, by
reditool-outTable-processed-analysis_final_replicates_combined_v2.py,
which reads the raw outTable files directly - mirroring how
edit_freq_all_positions_v2.py used to do this parsing in one batch step.

"""

import os, sys, subprocess, glob

# ---------- INPUT ARG ----------
if len(sys.argv) < 2:
    print("Usage: python mito_offtarget_single_sample_v2.py <sample_id>")
    sys.exit(1)

sample = sys.argv[1]
print("\n==============================")
print(" Processing sample: %s" % sample)
print("==============================")

# ---------- CONFIGURATION ----------
# NOTE: replace the placeholders below with your own paths before running.
fastq_dir = "<YOUR_DATA_DIR>"                     # e.g. "/path/to/mtDNA_WGS/DdCBE"
mt_reference = "<YOUR_REDITOOLS_DIR>/Homo_sapiens.GRCh38_MT.fa"
reditool = "<YOUR_REDITOOLS_DIR>/reditools/REDItoolDenovo.py"

mt_chrom = "MT"
bwa = "bwa"
samtools = "samtools"
threads = "7"

# ---------- FIND FASTQs ----------
# Files are gzipped .fastq.gz, e.g.:
#   ND1-DdCBE-TFL-rep1_S7_L001_R1_001.fastq.gz
#   ND4-DdCBE-TFH-rep1_S1_R1_001.fastq.gz          (no _L001_ token)
#   UT-rep1_S12_R1_001.fastq.gz
# `sample` (the SLURM array argument / samples.txt entry) is expected to be
# everything before "_R1_"/"_R2_", e.g. "ND1-DdCBE-TFL-rep1_S7_L001" or
# "UT-rep1_S12".
fastq_r1_list = glob.glob(os.path.join(fastq_dir, sample + "_R1_*.fastq.gz"))
fastq_r2_list = glob.glob(os.path.join(fastq_dir, sample + "_R2_*.fastq.gz"))

if len(fastq_r1_list) == 0 or len(fastq_r2_list) == 0:
    print("FASTQ files not found for %s" % sample)
    print("Checked path:", os.path.join(fastq_dir, sample + "_R1_*.fastq.gz"))
    sys.exit(1)

fastq_r1 = fastq_r1_list[0]
fastq_r2 = fastq_r2_list[0]

print("Using FASTQs:")
print("R1 -> " + fastq_r1)
print("R2 -> " + fastq_r2)

# ---------- OUTPUT DIR ----------
output_dir = sample + "_output_revised"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

bam_file = os.path.join(output_dir, "aligned_MT.sorted.bam")

def run_cmd(cmd, stdout=None):
    print(">>> " + " ".join(cmd))
    if stdout is not None:
        subprocess.check_call(cmd, stdout=stdout)
    else:
        subprocess.check_call(cmd)

# ---------- STEP 1: ALIGNMENT ----------
if not os.path.exists(bam_file):
    print("Running BWA-MEM alignment (paired-end, single call for R1+R2)...")

    raw_sam = os.path.join(output_dir, sample + "_aligned.sam")
    sam_out = open(raw_sam, "w")
    # Both mates aligned together in one bwa mem call, as in the paper's
    # Methods, instead of aligning R1 and R2 separately and merging BAMs.
    subprocess.check_call([bwa, "mem", "-t", threads, mt_reference, fastq_r1, fastq_r2], stdout=sam_out)
    sam_out.close()

    print("Converting SAM -> BAM...")
    raw_bam = os.path.join(output_dir, sample + "_aligned.bam")
    run_cmd([samtools, "view", "-bS", raw_sam, "-o", raw_bam])

    print("Name-sorting for fixmate...")
    namesorted_bam = os.path.join(output_dir, sample + "_namesorted.bam")
    run_cmd([samtools, "sort", "-n", "-@", threads, "-o", namesorted_bam, raw_bam])

    print("Running samtools fixmate (fixes mate-pair info/flags, per Methods)...")
    fixmate_bam = os.path.join(output_dir, sample + "_fixmate.bam")
    run_cmd([samtools, "fixmate", "-m", namesorted_bam, fixmate_bam])

    print("Coordinate-sorting and indexing final BAM...")
    run_cmd([samtools, "sort", "-@", threads, "-o", bam_file, fixmate_bam])
    run_cmd([samtools, "index", bam_file])

    for f in [raw_sam, raw_bam, namesorted_bam, fixmate_bam]:
        if os.path.exists(f):
            os.remove(f)
else:
    print("Skipping alignment (sorted BAM already exists).")

# ---------- STEP 2: ALIGNMENT STATS ----------
print("\nChecking MT mapping statistics...")
try:
    idxstats = subprocess.check_output([samtools, "idxstats", bam_file])
    for line in idxstats.splitlines():
        if line.startswith(mt_chrom):
            parts = line.split("\t")
            mapped = int(parts[2])
            length = parts[1]
            print("MT length: %s | Mapped reads: %s" % (length, mapped))
            if mapped == 0:
                print("No reads mapped to MT - skipping REDItools.")
                sys.exit(0)
except Exception as e:
    print("Could not obtain mapping statistics: %s" % str(e))
    sys.exit(1)

# ---------- STEP 3: REDItools ----------
# REDItoolDenovo.py does per-position pileup reads against the BAM, which
# means constant small I/O requests. Running this directly against a
# network filesystem (e.g. /lab/raguram_lab/...) is much slower than local
# disk, especially at high mtDNA coverage and/or with many concurrent array
# tasks contending for the same network storage. To avoid this bottleneck,
# the BAM/index/reference are copied to node-local scratch (/tmp) first,
# REDItools runs entirely on local disk, and only the (small) output table
# is copied back to the shared output_dir afterward.
reditools_out_dir = os.path.join(output_dir, "reditools_MT_revised")
if not os.path.exists(reditools_out_dir):
    print("Running REDItools Denovo on MT only (via local scratch)...")

    scratch_tag = os.environ.get("SLURM_JOB_ID", "nojob") + "_" + os.environ.get("SLURM_ARRAY_TASK_ID", "na")
    scratch_dir = os.path.join("/tmp", "reditools_scratch_%s_%s" % (sample, scratch_tag))
    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    local_bam = os.path.join(scratch_dir, "aligned_MT.sorted.bam")
    local_bai = local_bam + ".bai"
    local_ref = os.path.join(scratch_dir, os.path.basename(mt_reference))
    local_out = os.path.join(scratch_dir, "reditools_MT_revised")

    print("Copying BAM/index/reference to local scratch: " + scratch_dir)
    run_cmd(["cp", bam_file, local_bam])
    if os.path.exists(bam_file + ".bai"):
        run_cmd(["cp", bam_file + ".bai", local_bai])
    else:
        run_cmd([samtools, "index", local_bam])
    run_cmd(["cp", mt_reference, local_ref])
    # Copy the .fai index too if it already exists next to the reference;
    # otherwise let REDItools/samtools generate one locally on first use.
    if os.path.exists(mt_reference + ".fai"):
        run_cmd(["cp", mt_reference + ".fai", local_ref + ".fai"])

    cmd = [
        "python", reditool,
        "-i", local_bam,
        "-f", local_ref,
        "-o", local_out,
        "-c", "10", "-q", "30", "-s", "2",
        "-v", "1",
        "-n", "0.001"
    ]
    try:
        run_cmd(cmd)
        print("Copying REDItools results back to: " + reditools_out_dir)
        run_cmd(["cp", "-r", local_out, reditools_out_dir])
    finally:
        # Clean up local scratch regardless of success/failure, so /tmp
        # doesn't fill up across many array tasks on the same node.
        run_cmd(["rm", "-rf", scratch_dir])
else:
    print("Skipping REDItools (output exists).")

# ---------- STEP 4: FIND REDItools OUTPUT ----------
print("\nLocating REDItools output table...")
search_pattern = os.path.join(reditools_out_dir, "denovo_*", "outTable_*")
matches = glob.glob(search_pattern)
if len(matches) == 0:
    print("No REDItools output found - skipping sample.")
    sys.exit(0)
reditools_table = matches[0]
print("Found REDItools output table: %s" % reditools_table)
print("Completed sample: " + sample)

# NOTE: Parsing/filtering the REDItools outTable (extracting BaseCount[A,C,G,T],
# computing C->T / G->A conversion, applying the >0.1% threshold, removing
# SNV/on-target sites) is now done ALL AT ONCE across every sample by
# reditool-outTable-processed-analysis_final_replicates_combined_revised.py,
# which reads directly from this outTable file.

