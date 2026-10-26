
### POTENTIAL FOR BIOFILM FORMATION IS WIDESPREAD AND INFLUENCED BY PRECIPITATION VARIABILITY 
### AND ENVIRONMENTAL EXTREMES ACROSS U.S. SOILS 

# Authors: Eduardo M. Choreño-Parra, Nitish Malhotra, Alejandra Rodríguez-Verdugo, 
# Steven D. Allison, and Kathleen K. Treseder.

# ===================================================================================================================== #

# Description: This Bash pipeline processes 119 paired-end soil metagenomic samples from NEON. It performs quality filtering, contaminant removal, 
# sequence assembly, CDS prediction, functional annotation, read mapping, and quantification of raw read counts assigned to predicted CDSs.

# Inputs: 
# - metagenome_urls_2017.txt: URLs for downloading the 119 selected paired-end metagenomes (238 URLs).

# Outputs:
# - Summary statistics from the metagenomic processing pipeline present in TableS2_summary_pipeline_stats.csv (CSV version of Table S2)
# - all_samples.emapper.annotations: Functional annotations of predicted coding sequences (CDSs) across all metagenomic samples. TSV version of Table S3.
# - all_genes_long_all_samples.tsv: raw read counts assigned to all predicted CDSs (biofilm-associated and non-biofilm-associated) across all metagenomic samples.
#   TSV version of Table S4.

# Execution:
# Run this pipeline after section 1 in “1. NEON_metagenomic_environmental_data.R” and before “3. Biofilm_gene_analysis.R”. 
# The processing steps are designed for execution on an HPC cluster using SLURM. Submit each batch script in the indicated order, 
# ensuring that preceding jobs have completed before running dependent steps.




====================================================================================================================================
====================================================================================================================================
#################################################### 1. Downloading sequences ######################################################


# Define project directory
PROJECT_DIR="$HOME/metagenomic_analysis"
mkdir -p "$PROJECT_DIR"

cd "$PROJECT_DIR"
mkdir -p "1.raw_reads"
cd "$PROJECT_DIR/1.raw_reads" # Upload metagenome_urls_2017.txt to this folder

# Download manually the metagenomes from the URLs
wget -i metagenome_urls_2017.txt


###############################################


============================================
#### Summary: Counting reads per sample ####

# Paired read counts in each file
output_file="read_paired_reads.txt"

# Empty the file if it already exists
> "$output_file"

# Loop and write to file
for file in *.fastq.gz; do
  count=$(gzip -cd "$file" | wc -l)
  echo "$file: $((count / 4)) reads" >> "$output_file"
done

# This information is summarized in Table S2 (raw_paired_reads)




====================================================================================================================================
====================================================================================================================================
######################################################## 2. Run FastP ##############################################################


# Install FastP if necessary
# conda install -c bioconda fastp


===============================================
################# FastP loop ##################

# Create the fastP script and paste the SLURM batch job below
# touch 2.fastp_code.sh

# Submit the completed script using
# sbatch 2.fastp_code.sh


#!/bin/bash

#SBATCH --job-name=fastp_processing   
#SBATCH -A YOUR_ACCOUNT              
#SBATCH -p standard                  
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=18           
#SBATCH --mem-per-cpu=10G            
#SBATCH --error=fastp_error.txt         

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/1.raw_reads"
out_dir="$PROJECT_DIR/2.clean_reads_fastp"
mkdir -p "$out_dir"

module load fastp/0.23.2

for inp_read1 in "$inp_dir"/*_R1.fastq.gz; do
  base_name=$(basename "$inp_read1" _R1.fastq.gz)  
  inp_read2="$inp_dir/${base_name}_R2.fastq.gz"
 
  if [[ ! -f "$inp_read2" ]]; then
    echo "Skipping $base_name: Missing _R2 file!"
    continue
  fi

  out_read1="${base_name}.R1.fastq.gz"
  out_read2="${base_name}.R2.fastq.gz"
  sample_name="${base_name}"
  
fastp -i "$inp_read1" -I "$inp_read2" -o "$out_dir/$out_read1" -O "$out_dir/$out_read2" \
        --detect_adapter_for_pe -j "$out_dir/$sample_name.json" -h "$out_dir/$sample_name.html"

echo "Completed: $sample_name"
done

echo "Fastp processing finished for all samples."

###############################################


# NOTE: Output filenames use .R1.fastq.gz and .R2.fastq.gz instead of
# _R1.fastq.gz and _R2.fastq.gz for compatibility with Kneaddata.




====================================================================================================================================
====================================================================================================================================
############################################ 3. Run Kneaddata (remove human sequences) #############################################


===============================================
########## 3.1. Install human genome ########## 

# Install Kneaddata if necessary
# pip install kneaddata

# Create the human genome download script and paste the SLURM batch job below
# touch 3.humangenome_download.sh

# Submit the completed script using
# sbatch 3.humangenome_download.sh


#!/bin/bash

#SBATCH --job-name=human_genome_download
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=30
#SBATCH --time=6:00:00
#SBATCH --error=human_genome_error.txt

PROJECT_DIR="$HOME/metagenomic_analysis"

cd "$PROJECT_DIR"
mkdir -p "3.human_genome"

kneaddata_database --download human_genome bowtie2 "$PROJECT_DIR/3.human_genome/"

###############################################


# NOTE: Adapter trimming was performed previously using FastP.



===============================================
############## 3.2. Kneaddata array ###########

# Create the Kneaddata script and paste the SLURM batch job below
# touch 3.kneaddata_code.sh

# Submit the completed script using
# sbatch 3.kneaddata_code.sh


#!/bin/bash

#SBATCH --job-name=kneaddata_processing_array
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=16G
#SBATCH --error=kneaddata_error.txt
#SBATCH --array=1-119 # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/2.clean_reads_fastp"
out_dir="$PROJECT_DIR/3.kneaddata"
mkdir -p "$out_dir"

module load kneaddata

read1=`ls "$inp_dir"/*R1.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
  # Define the corresponding R2 file
  read2="${read1/R1.fastq.gz/R2.fastq.gz}"

  # Check if the R2 file exists
  if [[ -f "$read2" ]]; then
    # Extract the base name for output
    base=$(basename "$read1" .R1.fastq.gz)

    kneaddata \
    --input "$read1" --input "$read2" \
    --output "$out_dir" \
    -db "$PROJECT_DIR/3.human_genome/hg37dec_v0.1.1.bt2" \
    --bypass-trim \
    --output-prefix "$base" \
    --threads 8

 else
    echo "Paired file for $read1 not found. Skipping."
  fi

###############################################


# NOTE: Kneaddata processes files ending in .R1.fastq.gz and .R2.fastq.gz.

# Check number of Kneaddata paired-end output files
ls "$PROJECT_DIR"/3.kneaddata/*paired_1.fastq | wc -l




=======================================================================================================================================
=======================================================================================================================================
############################ 4. FastQC and MultiQC to check quality of each individual file ###########################################


===============================================
################ 4.1. FastQC array ############

# Create the FastQC script and paste the SLURM batch job below
# touch 4.fastqc_code.sh

# Submit the completed script using
# sbatch 4.fastqc_code.sh


#!/bin/bash

#SBATCH --job-name=fastqc_processing
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=16G
#SBATCH --error=fastqc_error.txt
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/3.kneaddata"
out_dir="$PROJECT_DIR/4.fast_multi_qc_files"
mkdir -p "$out_dir"

module load fastqc

# Select one file per SLURM job
read1=$(ls "$inp_dir"/*_paired_1.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# Define the corresponding R2 file
read2="${read1/_paired_1.fastq/_paired_2.fastq}"

# Check if R2 exists before running FastQC
if [[ -f "$read2" ]]; then
    fastqc -o "$out_dir" "$read1"
    fastqc -o "$out_dir" "$read2"
    echo "Processed: $(basename "$read1") and $(basename "$read2")"
else
    echo "Missing R2 pair for $(basename "$read1"). Skipping."
fi

###############################################


===============================================
################# 4.2 MultiQC ### #############

# Install MultiQC in a Conda environment if necessary
# conda create -n multiqc_env -c conda-forge -c bioconda multiqc

# Create the MultiQC script and paste the SLURM batch job below
# touch 4.2.multiqc_code.sh

# Submit the completed script using
# sbatch 4.2.multiqc_code.sh


#!/bin/bash

#SBATCH --job-name=multiqc_run
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --error=multiqc_error.txt

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/4.fast_multi_qc_files"
out_dir="$PROJECT_DIR/4.fast_multi_qc_files"

# Load MultiQC environment
module load anaconda3
source activate multiqc_env

# Run MultiQC after FastQC finishes
multiqc -o "$out_dir" "$inp_dir"

echo "MultiQC completed."

###############################################




=======================================================================================================================================
=======================================================================================================================================
########################################### 5. MEGAHIT to assemble clean sequences ####################################################


=================================================
#### 5.1. MEGAHIT assembly array (per sample) ###

# Install MEGAHIT if necessary
# conda install -c bioconda megahit

# Create the MEGAHIT script and paste the SLURM batch job below
# touch 5.megahit_code_single_array.sh

# Submit the completed script using
# sbatch 5.megahit_code_single_array.sh


#!/bin/bash

#SBATCH --job-name=megahit_single_array
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=35G
#SBATCH --error=megahit_single_array_error.txt
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/3.kneaddata"
out_dir="$PROJECT_DIR/5.contigs_megahit_single_array"
mkdir -p "$out_dir"

module load megahit/1.2.9

# Select one file per SLURM job
read1=$(ls "$inp_dir"/*_paired_1.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# Define the corresponding R2 file
read2="${read1/_paired_1.fastq/_paired_2.fastq}"

# Extract sample name
sample_name=$(basename "$read1" _paired_1.fastq)

# Run MEGAHIT per sample
megahit -1 "$read1" \
        -2 "$read2" \
        -o "$out_dir/${sample_name}_megahit_out" \
        --min-contig-len 1000 \
        -t 12

echo "MEGAHIT assembly completed for ${sample_name}"

###############################################



======================================
#### Summary: assembly statistics ####

PROJECT_DIR="$HOME/metagenomic_analysis"
cd "$PROJECT_DIR/5.contigs_megahit_single_array"

# Create a tab-delimited file with one line per sample
echo -e "sample\tstats" > megahit_contig_stats.tsv

for log in *_megahit_out/log; do
    sample=$(basename "$(dirname "$log")")
    line=$(grep "contigs, total" "$log" | tail -n 1)
    echo -e "${sample}\t${line}" >> megahit_contig_stats.tsv
done

# Assembly statistics (total contigs, N50, Avg length, and total bp) are summarized in Table S2




=======================================================================================================================================
=======================================================================================================================================
############################################## 6. Predict ORFs using Prokka  ##########################################################


===============================================
######## 6.1. Prokka prediction array #########


# Install Prokka in a Conda environment if necessary
# conda create -n prokka -c bioconda prokka
# conda activate prokka

# Create the Prokka script and paste the SLURM batch job below
# touch 6.Prokka_ORF_prediction.sh

# Submit the completed script using
# sbatch 6.Prokka_ORF_prediction.sh


#!/bin/bash

#SBATCH --job-name=Prokka_ORF_Prediction_Array
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --error=prokka_array_error.txt
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
inp_dir="$PROJECT_DIR/5.contigs_megahit_single_array"
out_dir="$PROJECT_DIR/6.prokka_output_single_array"
mkdir -p "$out_dir"

# Activate Prokka environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prokka

# Select one sample per SLURM task
assembly=$(ls "$inp_dir"/*_megahit_out/final.contigs.fa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

if [[ -z "$assembly" ]]; then
    echo "No assembly found for array ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample name
sample_name=$(basename $(dirname "$assembly"))

# Run Prokka on the individual assembly
prokka --outdir "$out_dir/${sample_name}_prokka_out" \
       --prefix "$sample_name" \
       --cpus 12 \
       --metagenome \
       "$assembly"

echo "Prokka ORF prediction completed for ${sample_name}"

###############################################


============================================
#### Summary: ORF prediction statistics ####   

PROJECT_DIR="$HOME/metagenomic_analysis"
cd "$PROJECT_DIR/6.prokka_output_single_array"

echo -e "sample\tcontigs\tbases\tCDS\trRNA\ttRNA" > prokka_annotation_stats.tsv

for d in *_prokka_out; do
    summary_file=$(ls "$d"/*.txt)
    sample=$(basename "$d" _prokka_out)

    # Extract annotation statistics
    contigs=$(grep "^contigs:" "$summary_file" | awk '{print $2}')
    bases=$(grep "^bases:" "$summary_file" | awk '{print $2}')
    cds=$(grep "^CDS:" "$summary_file" | awk '{print $2}')
    rrna=$(grep "^rRNA:" "$summary_file" | awk '{print $2}')
    trna=$(grep "^tRNA:" "$summary_file" | awk '{print $2}')

    # Replace empty fields with NA
    contigs=${contigs:-NA}
    bases=${bases:-NA}
    cds=${cds:-NA}
    rrna=${rrna:-NA}
    trna=${trna:-NA}

    echo -e "${sample}\t${contigs}\t${bases}\t${cds}\t${rrna}\t${trna}" >> prokka_annotation_stats.tsv
done

# Number of CDSs per sample is summarized in Table S2




=======================================================================================================================================
=======================================================================================================================================
############################################## 7. Functional annotation using eggNOG-mapper  ########################################## 


===============================================
########### 7.1 Download eggnog db ############

# Install eggNOG-mapper in a Conda environment if necessary
# conda create -n eggnog -c conda-forge -c bioconda eggnog-mapper -y

# Create the eggNOG database download script and paste the SLURM batch job below
# touch 7.1.eggnog_db_download_wget.sh

# Submit the completed script using
# sbatch 7.1.eggnog_db_download_wget.sh


#!/bin/bash

#SBATCH --job-name=eggnog_db_dl
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=eggnog_db_download.out
#SBATCH --error=eggnog_db_download.err

set -euo pipefail

PROJECT_DIR="$HOME/metagenomic_analysis"
DBDIR="$PROJECT_DIR/databases/eggnog_mapper"
mkdir -p "$DBDIR"
cd "$DBDIR"

BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

echo "Downloading to: $DBDIR"
date

# Resume-capable downloads
wget -c "${BASE_URL}/eggnog.db.gz"
wget -c "${BASE_URL}/eggnog_proteins.dmnd.gz"
wget -c "${BASE_URL}/eggnog.taxa.tar.gz"

echo "Decompressing / extracting..."
gunzip -f eggnog.db.gz
gunzip -f eggnog_proteins.dmnd.gz
tar -xzf eggnog.taxa.tar.gz
rm -f eggnog.taxa.tar.gz

echo "Finished. Key files:"
ls -lh eggnog.db eggnog_proteins.dmnd | sed -n '1,5p'
date


===============================================
######## 7.2. eggNOG annotation of ORFs #######

# Create the eggNOG-mapper script and paste the SLURM batch job below
# touch 7.2.eggnog_orf.sh

# Submit the completed script using
# sbatch 7.2.eggnog_orf.sh


#!/bin/bash

#SBATCH --job-name=eggnog_orf
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)
#SBATCH --output=eggnog.out
#SBATCH --error=eggnog.err

set -eo pipefail

PROJECT_DIR="$HOME/metagenomic_analysis"
DATA_DIR="$PROJECT_DIR/databases/eggnog_mapper"
prokka_dir="$PROJECT_DIR/6.prokka_output_single_array"
out_dir="$PROJECT_DIR/7.2.eggnog_orf"
mkdir -p "$out_dir"

# Activate eggNOG environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate eggnog

# Select one Prokka protein FASTA (.faa) per task
faa_file=$(ls "$prokka_dir"/*_prokka_out/*_megahit_out.faa | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$faa_file" || ! -f "$faa_file" ]]; then
  echo "No .faa found for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

sample_base=$(basename "$faa_file" .faa)

echo "Running eggNOG-mapper on: $faa_file"
echo "Sample base: $sample_base"

# Run eggNOG-mapper
emapper.py \
  -i "$faa_file" \
  --itype proteins \
  --cpu 16 \
  --data_dir "$DATA_DIR" \
  --output "$sample_base" \
  --output_dir "$out_dir" \
  --override

echo "Done eggNOG for: $sample_base"

##############################################


===============================================
########## 7.3. Merge eggNOG annotations #######
# This produces information contained in Table S3

PROJECT_DIR="$HOME/metagenomic_analysis"
out_dir="$PROJECT_DIR/7.2.eggnog_orf"
out_tsv="$out_dir/all_samples.emapper.annotations.tsv"

first=$(ls "$out_dir"/*.emapper.annotations | head -n 1)

# Write correct header (convert "#query" to "query")
{
  printf "sample\t"
  grep -m1 '^#query' "$first" | sed 's/^#//'
} > "$out_tsv"

# Append annotation rows from all samples (skip lines starting with '#')
for f in "$out_dir"/*.emapper.annotations; do
  sample=$(basename "$f" .emapper.annotations)
  grep -v '^#' "$f" | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}'
done >> "$out_tsv"

echo "Wrote: $out_tsv"




======================================================================================================================================
======================================================================================================================================
################################# 8. Build indices and align quality-filtered reads to annotated contigs #############################


===============================================
###### 8.1. Build Bowtie2 index by sample #####

# Create the Bowtie2 script and paste the SLURM batch job below
# touch 8.1.Bowtie2_index_by_sample.sh

# Submit the completed script using
# sbatch 8.1.Bowtie2_index_by_sample.sh


#!/bin/bash

#SBATCH --job-name=Bowtie2_index_by_sample
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --error=bowtie2_index_by_sample_error.txt
#SBATCH --output=bowtie2_index_by_sample.out
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
prokka_dir="$PROJECT_DIR/6.prokka_output_single_array"
bowtie_out="$PROJECT_DIR/8.bowtie2_index_by_sample"
mkdir -p "$bowtie_out"

# Load Bowtie2
module load bowtie2/2.5.1

# Select the corresponding Prokka nucleotide FASTA (.fna) for this sample
FASTA_FILE=$(ls "$prokka_dir"/*_prokka_out/*.fna | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract sample base name (same pattern as phmmer)
sample_base=$(basename "$FASTA_FILE" .fna)

# Define index prefix
INDEX_PREFIX="$bowtie_out/${sample_base}_index"

# Build Bowtie2 index
bowtie2-build "$FASTA_FILE" "$INDEX_PREFIX"

echo "Bowtie2 index created for $sample_base"

###############################################


===============================================
######## 8.2. Align qf reads to contigs #######

# Create the Bowtie2 mapping script and paste the SLURM batch job below
# touch 8.2.bowtie2_mapping_by_sample.sh

# Submit the completed script using
# sbatch 8.2.bowtie2_mapping_by_sample.sh


#!/bin/bash

#SBATCH --job-name=Bowtie2_mapping_by_sample
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --error=bowtie2_mapping_by_sample_error.txt
#SBATCH --output=bowtie2_mapping_by_sample.out
#SBATCH --array=1-119  # One task per metagenomic sample (119 samples selected for processing)

PROJECT_DIR="$HOME/metagenomic_analysis"
READS_DIR="$PROJECT_DIR/3.kneaddata"
INDEX_DIR="$PROJECT_DIR/8.bowtie2_index_by_sample"
OUT_DIR="$PROJECT_DIR/8.bowtie2_mapping_by_sample"
mkdir -p "$OUT_DIR"

module load bowtie2/2.5.1
module load samtools/1.15.1

# Get the sample for this array task
READ1=$(ls "$READS_DIR"/*_paired_1.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")
READ2="${READ1/_paired_1.fastq/_paired_2.fastq}"

# Extract sample name
SAMPLE_NAME=$(basename "$READ1" _paired_1.fastq)

# Define index prefix
INDEX_PREFIX="$INDEX_DIR/${SAMPLE_NAME}_megahit_out_index"

echo "Mapping sample: $SAMPLE_NAME"
echo "Using R1: $READ1"
echo "Using R2: $READ2"
echo "Using index: $INDEX_PREFIX"

# Run Bowtie2 alignment
bowtie2 -x "$INDEX_PREFIX" \
        -1 "$READ1" \
        -2 "$READ2" \
        -S "$OUT_DIR/${SAMPLE_NAME}.sam"

# Convert SAM to BAM, sort and index
samtools view -bS "$OUT_DIR/${SAMPLE_NAME}.sam" | samtools sort -o "$OUT_DIR/${SAMPLE_NAME}_sorted.bam"
samtools index "$OUT_DIR/${SAMPLE_NAME}_sorted.bam"

# Remove large SAM file to save space
rm "$OUT_DIR/${SAMPLE_NAME}.sam"

echo "Bowtie2 mapping completed for $SAMPLE_NAME"

###############################################


=============================================== 
###### 8.2.1. BAM statistics per sample  ######

# Create the script and paste the SLURM batch job below
# touch 8.2.1.mapping_stats.sh

# Submit the script
# sbatch 8.2.1.mapping_stats.sh


#!/bin/bash

#SBATCH --job-name=mapping_bam_stats
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-119
#SBATCH --output=8.2.1.bam_stats.out
#SBATCH --error=8.2.1.bam_stats.txt

PROJECT_DIR="$HOME/metagenomic_analysis"

BAM_DIR="$PROJECT_DIR/8.bowtie2_mapping_by_sample"
STATS_DIR="$BAM_DIR/full_stats"
mkdir -p "$STATS_DIR"

module load samtools/1.15.1

# Get BAM
BAM=$(ls "$BAM_DIR"/*_sorted.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=$(basename "$BAM" _sorted.bam)

# Output files
FLAGSTAT_OUT="$STATS_DIR/${SAMPLE}_flagstat.txt"

# Basic flagstat summary
samtools flagstat "$BAM" > "$FLAGSTAT_OUT"


=============================================== 
########## 8.2.2. Flagstat summary  ###########

# Create the script and paste the SLURM batch job below
# touch 8.2.2.flagstat_summary.sh

# Submit the script after Section 8.2.1 finishes:
# sbatch 8.2.2.flagstat_summary.sh


#!/bin/bash

#SBATCH --job-name=flagstat_summary
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=8.2.2.flagstat_summary.out
#SBATCH --error=8.2.2.flagstat_summary.txt

PROJECT_DIR="$HOME/metagenomic_analysis"
stats_dir="$PROJECT_DIR/8.bowtie2_mapping_by_sample/full_stats"
out_tsv="$PROJECT_DIR/8.2.2.flagstat_summary.tsv"

# Header
echo -e "sample\ttotal_reads\tprimary_reads\tmapped_reads\tmapped_pct\tprimary_mapped\tpaired_reads\tproperly_paired\tboth_mapped\tsingletons\tmates_diff_chr\tmates_diff_chr_mapq_ge_5" > "$out_tsv"

for flagstat in "$stats_dir"/*_flagstat.txt; do
    sample=$(basename "$flagstat" _flagstat.txt)

    total_reads=$(      head -n1 "$flagstat"                               | awk '{print $1}')
    primary_reads=$(    grep ' primary$' "$flagstat"                       | awk '{print $1}')
    mapped_reads=$(     grep ' mapped (' "$flagstat" | head -n1            | awk '{print $1}')
    primary_mapped=$(   grep ' primary mapped (' "$flagstat"               | awk '{print $1}')
    paired_reads=$(     grep ' paired in sequencing' "$flagstat"           | awk '{print $1}')
    properly_paired=$(  grep ' properly paired (' "$flagstat"              | awk '{print $1}')
    both_mapped=$(      grep ' with itself and mate mapped' "$flagstat"    | awk '{print $1}')
    singletons=$(       grep ' singletons (' "$flagstat"                   | awk '{print $1}')
    mates_diff_chr=$(   grep ' with mate mapped to a different chr$' "$flagstat"        | awk '{print $1}')
    mates_diff_chr_q5=$(grep ' with mate mapped to a different chr (mapQ>=5)' "$flagstat" | awk '{print $1}')

    if [[ "$total_reads" -gt 0 ]]; then
        mapped_pct=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{printf "%.4f", (m*100.0)/t}')
    else
        mapped_pct="NA"
    fi

    echo -e "${sample}\t${total_reads}\t${primary_reads}\t${mapped_reads}\t${mapped_pct}\t${primary_mapped}\t${paired_reads}\t${properly_paired}\t${both_mapped}\t${singletons}\t${mates_diff_chr}\t${mates_diff_chr_q5}" >> "$out_tsv"
done

echo "✅ Flagstat summary table written to:"
echo "$out_tsv"

##################################################

# Note: Information in flagstat_summary.tsv (total_reads, mapped_reads, mapped_pct/mapping_rate) is summarized in Table S2.



======================================================================================================================================
======================================================================================================================================
######################################## 9. Count gene abundance per sample ##########################################################


=================================================
########## 9.1. FeatureCounts: All genes ########

# Create the script and paste the SLURM batch job below
# touch 9.featureCounts_all_genes_by_sample.sh

# Submit the script
# sbatch 9.featureCounts_all_genes_by_sample.sh


#!/bin/bash

#SBATCH --job-name=featureCounts_all_array
#SBATCH -A YOUR_ACCOUNT
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --time=04:00:00
#SBATCH --output=featureCounts_all.out
#SBATCH --error=featureCounts_all.err
#SBATCH --array=1-119

set -eo pipefail

module load subread/2.0.3
module load samtools/1.15.1

PROJECT_DIR="$HOME/metagenomic_analysis"
BAM_DIR="$PROJECT_DIR/8.bowtie2_mapping_by_sample"
GFF_DIR="$PROJECT_DIR/6.prokka_output_single_array"
RESULTS_DIR="$PROJECT_DIR/9.gene_counts_by_sample_all"
mkdir -p "$RESULTS_DIR"

# Select BAM for this task
BAM_FILE=$(ls "$BAM_DIR"/*_sorted.bam | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
  echo "ERROR: No BAM found for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

SAMPLE_NAME=$(basename "$BAM_FILE" _sorted.bam)

GFF_FILE="$GFF_DIR/${SAMPLE_NAME}_megahit_out_prokka_out/${SAMPLE_NAME}_megahit_out.gff"

if [[ ! -f "$GFF_FILE" ]]; then
  echo "ERROR: GFF not found: $GFF_FILE"
  exit 1
fi

echo "Sample: $SAMPLE_NAME"
echo "BAM: $BAM_FILE"
echo "GFF: $GFF_FILE"

# Run featureCounts (counts all CDS features)
OUT_COUNTS="$RESULTS_DIR/${SAMPLE_NAME}_gene_counts.txt"

featureCounts \
  -T 12 \
  -p \
  -t CDS \
  -g ID \
  -a "$GFF_FILE" \
  -o "$OUT_COUNTS" \
  "$BAM_FILE"

echo "featureCounts done: $OUT_COUNTS"

# Convert to long format (all genes)
# featureCounts output: header line then column names line then data
# Columns: Geneid Chr Start End Strand Length <sample.bam>
OUT_COUNTS="$RESULTS_DIR/${SAMPLE_NAME}_gene_counts.txt"
OUT_LONG="$RESULTS_DIR/${SAMPLE_NAME}_all_genes_long.tsv"

echo "DEBUG: OUT_COUNTS=$OUT_COUNTS"
ls -lh "$OUT_COUNTS" || { echo "ERROR: counts file missing"; exit 1; }

echo "DEBUG: First 3 lines of counts file:"
head -n 3 "$OUT_COUNTS"

# Write header
echo -e "sample\torf_id\tcontig\tstart\tend\tstrand\tlength\tcount" > "$OUT_LONG"

# Parse counts file (skip metadata + header)
awk -F'\t' -v OFS='\t' -v sample="$SAMPLE_NAME" '
  /^#/ { next }                    # skip Program line(s)
  $1=="Geneid" { next }            # skip header line
  NF>=7 {                          # require at least 7 columns
    print sample, $1, $2, $3, $4, $5, $6, $7
  }
' "$OUT_COUNTS" >> "$OUT_LONG"

echo "Long table done: $OUT_LONG"

###############################################


===============================================
###### Merge gene counts across samples #######

PROJECT_DIR="$HOME/metagenomic_analysis"
RESULTS_DIR="$PROJECT_DIR/9.gene_counts_by_sample_all"
OUT_ALL="$RESULTS_DIR/all_genes_long_all_samples.tsv"

# Write header
echo -e "sample\torf_id\tcontig\tstart\tend\tstrand\tlength\tcount" > "$OUT_ALL"

# Append data from each file (skip header line in each)
for f in "$RESULTS_DIR"/*_all_genes_long.tsv; do
  awk 'NR>1' "$f" >> "$OUT_ALL"
done

echo "Merged file: $OUT_ALL"
echo "Lines (including header): $(wc -l < "$OUT_ALL")"

###############################################

# Note: all_genes_long_all_samples.tsv contains raw read counts assigned to all predicted CDSs (biofilm-associated and non-biofilm-associated) 
# across all samples. This database is used to calculate total ORF-assigned reads per sample and the basis to generate Table S4.























=============================================== 
####### 9.3. Counts (reads) all genes  ########### #Counts number of reads matching to ORFs (including non-bioifilm and biofilm)

BASE_DIR="/dfs7/treseder_lab/echoreno/neon/2017_bacteria_allsites"
RESULTS_DIR="$BASE_DIR/9.gene_counts_by_sample_all"
OUT="$RESULTS_DIR/total_orf_assigned_fragments.tsv"

echo -e "sample\torf_assigned_fragments" > "$OUT"

for fc in "$RESULTS_DIR"/*_gene_counts.txt; do

  # safety: skip summary files if present
  [[ "$fc" == *.summary ]] && continue

  sample=$(basename "$fc" _gene_counts.txt)

  total=$(awk -F'\t' '
    /^#/ { next }            # skip metadata lines
    $1=="Geneid" { next }    # skip header
    NF>=7 { sum += $7 }      # counts column
    END { print sum+0 }
  ' "$fc")

  echo -e "${sample}\t${total}" >> "$OUT"
done

echo "✅ Wrote: $OUT"









################################################ STATISTICS ###################################################################

#### Total raw reads

cd /dfs7/treseder_lab/echoreno/neon/2017_allsites_fungal_genomes/1.raw_reads/

echo -e "SampleID\tRaw_R1_reads" > raw_read_counts.tsv

for f in *_R1*.fastq.gz; do
    sample=$(basename "$f" | sed 's/_R1.*.fastq.gz//')
    reads=$(zcat "$f" | wc -l)
    reads=$((reads/4))
    echo -e "${sample}\t${reads}" >> raw_read_counts.tsv
done



#### Clean reads

touch clean_read_counts.sh

#!/bin/bash

#SBATCH --job-name=clean_read_counts
#SBATCH -A treseder_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=05:00:00
#SBATCH --output=clean_read_countsout
#SBATCH --error=clean_read_countserr
#SBATCH --mail-user=echoreno@uci.edu
#SBATCH --mail-type=ALL


cd /dfs7/treseder_lab/echoreno/neon/2017_allsites_fungal_genomes/3.kneaddata/

echo -e "SampleID\tClean_R1_reads" > clean_read_counts.tsv

for f in *_paired_1.fastq; do
    sample=$(basename "$f" | sed 's/_paired_1.fastq//')
    reads=$(cat "$f" | wc -l)
    reads=$((reads/4))
    echo -e "${sample}\t${reads}" >> clean_read_counts.tsv
done




