#!/bin/bash
# ============================================================================
# PURPOSE
#   Aim to analyse and process long read nanopore sequencing data from Transgenic A.gambaie mosquitio strains 
#   to identify the location of the transgene in the host native genome. 
#
# LOGIC OF THE PIPELINE
#   Raw reads are aligned to the Ifakara reference so that specfic reads touching the
#   artifical transgene can be isolated and visualsied in IGV. 
#   This allows for visual conformation of the transgene presence in the Trasgeneic mosquito genome. 
#
# HOW TO USE THE PIPELINE Script. 
#   The use of this pipeline requires access to the shared High Performance Cluster at Imperial owned by the Crisanti Lab 
#   An analysis file needs to be created in Crisanti ephemeral directory *name as perferred*
#   Here the script will be run and the data stored to save disk space until it can be moved to a local storage drive.
# 
# The analysis file should contain: 
#   - Ref folder with host ifakara geneome. This is the reference genome for Anopholes Gambaie. 
#   - Long read Nanopore sequencing file in ziped format 
#   - Analysis_pipeline.sh script file 
#
#  ****REVIEW THE SCRIPT AND ADJUST THE FOLLOWING VARIABLES ****
#   - SAMPLE="" inside the quotes should be the sample of your FASTQ file you want to analyse 
#   - TARGET_CHROMOSOME=""inside the quotes should be the name of your artifical chromosome as named in your reference genome FASTA file located in ref folder.  
#   - BASE_DIR="" This should be the location in ephemeral where the analysis file with its required content is stored eg BASE_DIR="/rds/general/project/tmstorage/ephemeral/"
#   - REFERENCE="$BASE_DIR/ref/VectorBase-68_AgambiaeIfakara_Genome.fasta" 
#
# Script inputs and outputs
#   INPUT   : Sample Nanopore FASTQ zip file, Host genome reference FASTA file with artifial transgene chromosome. 
#   OUTPUT  : Assorted BAM files that map reads to the transgene and the refernce file. Summary statistics on the coverage of the transgene. 
# ============================================================
#PBS -N AnalysisPipeline
#PBS -l select=1:ncpus=16:mem=128gb
#PBS -l walltime=08:00:00 
#PBS -j oe 

# ============================================================
#  CONDA Enviroment 
source /rds/general/user/nkhan3/home/miniforge3/etc/profile.d/conda.sh
conda activate Minimap2_miniforge
# NOTE: this env lives in another user's home directory. It works while their
# permissions allow it, but your pipeline breaks silently if they rebuild or
# move it. Worth cloning into your own space for reproducibility:
#     conda env export -n Minimap2_miniforge > Minimap2_miniforge.yml
#     conda env create -f Minimap2_miniforge.yml
#
# Exit immediately if a command exits with a non-zero status,
# -e           : abort on any command failure
# -u           : abort on use of an undefined variable (catches typos in paths)
# -o pipefail  : abort if ANY command in a pipe fails, not just the last one.
#                Critical here: without it, `minimap2 | samtools sort` would
#                report success even if minimap2 crashed, leaving a truncated BAM.
set -euo pipefail  

# ============================================================================
# SAMPLE & RESOURCE CONFIGURATION
#
# TARGET_CHROMOSOME must match a sequence NAME inside the reference FASTA
# exactly, not the reference's filename. Verified in the preflight below.
# ============================================================================
SAMPLE="Sample1" # ADJUST TO THE NUMBER OF SAMPLES TO ANALYSE ("" "NAME SPECFIC") 
TARGET_CHROMOSOME="vasa2cas9_3xP3_GFP" # Named exactly as the artifical chromosome in the reference FASTA 
THREADS=16
echo "=========================================="
echo "Running sample  : $SAMPLE"
echo "Started         : $(date)"
echo "=========================================="


# ============================================================
# PATHS — update these if needed, Must First create a file contianed samples, and ref folder+content !
## All outputs are nested under BASE_DIR so a run is self-contained.
# WARNING: BASE_DIR is on RDS `ephemeral`, which is purged on a fixed schedule.
# Anything you need to keep must be copied to `live` or your home space before
# that window elapses. Confirm the current retention period against Imperial
# RCS documentation — do not assume.
#============================================================
# Correct path pointing directly to the project's system-managed ephemeral folder
BASE_DIR="/rds/general/user/nkhan3/projects/tmstorage/ephemeral/vasa2cas9_analysis_testing"


# Create your job directory (using -p ensures it succeeds safely)
mkdir -p "$BASE_DIR"
# ============================================================
# Input what the sample FASTA sequecing file name is between the ""
# ============================================================
REFERENCE="$BASE_DIR/ref/VectorBase-68_AgambiaeIfakara_Genome.fasta" # ADD the transgene sequence into the fasta file  
INDEX="$BASE_DIR/ref/VectorBase-68_AgambiaeIfakara_Genome.mmi" # will get re built if missing or outdated
READS="$BASE_DIR/${SAMPLE}.fastq.gz" # this the zip file of the reads for each sample, make sure the path is correct and the file exists.

# Output prefix and dirs are sample-specific so they don't overwrite each other
OUTPUT_PREFIX="${SAMPLE}_results" # Prefix for output files (e.g., BAM, coverage) to keep them organized by sample
RESULTS_DIR="$BASE_DIR/${SAMPLE}_results" # Directory to store intermediate results (e.g., BAM files, coverage files) for each sample
REPORT_DIR="$BASE_DIR/${SAMPLE}_report" # Directory to store final reports for each sample
CHROM_DIR="$REPORT_DIR/chromosome_specific" # Subdirectory within the report directory for chromosome-specific analyses (e.g., coverage plots, stats for target chromosomes)


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================
mkdir -p "$RESULTS_DIR" # creates a results firectory for each sample, where the output of minimap2 and samtools will be stored.
mkdir -p "$REPORT_DIR" # creates a report directory for each sample, where the final report and chromosome-specific analyses will be stored.
mkdir -p "$CHROM_DIR" # creates a subdirectory within the report directory for chromosome-specific analyses.


echo "=========================================="
echo "Running sample  : $SAMPLE"
echo "Date & Time     : $(date)"
echo "=========================================="
echo "Reference        : $REFERENCE"
echo "Reads            : $READS"
echo "Results dir      : $RESULTS_DIR"
echo "Report dir       : $REPORT_DIR"
echo "Target chromosome: $TARGET_CHROMOSOME"
echo "Threads          : $THREADS"
echo "=========================================="
# ============================================================

# ============================================================================
# PREFLIGHT: Double check if the tools required for the script to run is present in the enviroment. Also if the 
#
# Each check below corresponds to a way v1.1 could fail late and expensively.
# ============================================================================
echo "--- PREFLIGHT ---"

# (a) Every binary the script calls must exist in the activated env. Without
#     this, a missing tool surfaces only when its step is reached.
for tool in minimap2 samtools ; do
    command -v "$tool" >/dev/null 2>&1 \
        || { echo "ERROR: '$tool' not found in the active conda env." >&2; exit 1; }
done
echo "All required tools present."

# (b) Inputs must exist and be non-empty. -s tests both in one operator.
[[ -s "$READS"     ]] || { echo "ERROR: The Zip fastQ file not found or empty: $READS" >&2; exit 1; }
[[ -s "$REFERENCE" ]] || { echo "ERROR: reference not found or empty: $REFERENCE" >&2; exit 1; }

# (c) TARGET_CHROMOSOME must be a real sequence name in the reference.
#     If it is not, `samtools view <bam> <region>` returns zero reads WITHOUT
#     erroring, and every subsequent step processes empty files. This is the
#     single most likely silent failure in the whole pipeline.
if ! grep -q "^>${TARGET_CHROMOSOME}\([[:space:]]\|$\)" "$REFERENCE"; then
    echo "ERROR: '$TARGET_CHROMOSOME' is not a sequence name in the reference." >&2
    echo "Available sequence names:" >&2
    grep "^>" "$REFERENCE" >&2
    exit 1
fi
echo "Target contig '$TARGET_CHROMOSOME' confirmed present in reference."
echo "--- PREFLIGHT PASSED ---"


# ============================================================
# STEP 1: Build minimap2 index when the Index file is older than the refeerence file. 
# ============================================================
if [ ! -f "$INDEX" ] || [ "$REFERENCE" -nt "$INDEX" ]; then
    echo "Rebuilding minimap2 index..."
    minimap2 -d "$INDEX" "$REFERENCE"
else
    echo "Index already up-to-date."
fi

echo "Rebuilding samtools FASTA index..."
samtools faidx "$REFERENCE"

# ============================================================
# STEP 2 & 3 & 4: Stream Alignment directly to Sorted BAM (Saves 70% Disk Space)
# ============================================================
echo "Aligning and sorting $SAMPLE in a single stream..."
minimap2 -ax map-ont -t $THREADS $INDEX $READS | \
samtools view -@ $THREADS -b -q 21 - | \
samtools sort -@ $THREADS -o $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam -

# ============================================================
# STEP 5: Index sorted BAM
# ============================================================
echo "Indexing sorted BAM..."
samtools index $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam

# ============================================================
# STEP 6: Alignment statistics
# ============================================================
echo "Generating alignment statistics..."
echo "=== BASIC ALIGNMENT STATISTICS ===" > $REPORT_DIR/${OUTPUT_PREFIX}_stats.txt
samtools flagstat $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam >> $REPORT_DIR/${OUTPUT_PREFIX}_stats.txt
samtools stats $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam > $REPORT_DIR/${OUTPUT_PREFIX}_detailed_stats.txt

# ============================================================
# STEP 7: Per-chromosome mapping statistics
# ============================================================
echo "=== READS MAPPED PER CHROMOSOME ===" > $REPORT_DIR/${OUTPUT_PREFIX}_per_chromosome.txt
samtools idxstats $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam >> $REPORT_DIR/${OUTPUT_PREFIX}_per_chromosome.txt

# ============================================================
# STEP 8: Coverage statistics
# ============================================================
echo "=== COVERAGE PER CHROMOSOME ===" > $REPORT_DIR/${OUTPUT_PREFIX}_coverage.txt
samtools coverage $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam >> $REPORT_DIR/${OUTPUT_PREFIX}_coverage.txt

# ============================================================
# STEP 9: Mapping quality distribution
# ============================================================
samtools view -q 1 $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam | cut -f5 | sort -n | uniq -c > $REPORT_DIR/${OUTPUT_PREFIX}_mapq_distribution.txt

# ============================================================
# STEP 10: Focused analysis on target chromosomes
# ============================================================
echo "Focused analysis on: $TARGET_CHROMOSOME"

# Extract all alignments for target chromosome
samtools view -@ $THREADS -h $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam $TARGET_CHROMOSOME \
    > $CHROM_DIR/${TARGET_CHROMOSOME}_all_alignments.sam

# Extract primary alignments and convert to FASTA
samtools view -@ $THREADS -F 256 $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam $TARGET_CHROMOSOME \
    | awk '{print ">"$1"\n"$10}' > $CHROM_DIR/${TARGET_CHROMOSOME}_primary_reads.fasta

# Extract secondary alignments
samtools view -@ $THREADS -f 256 $RESULTS_DIR/${OUTPUT_PREFIX}_sorted.bam $TARGET_CHROMOSOME \
    > $CHROM_DIR/${TARGET_CHROMOSOME}_secondary_alignments.sam

echo "Job completed successfully for $SAMPLE!"
