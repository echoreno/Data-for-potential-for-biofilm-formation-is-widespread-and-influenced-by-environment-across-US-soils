## SCRIPTS

This folder contains the R and Bash scripts used for metagenomic processing, environmental data preparation, and statistical analyses.

Scripts are numbered according to their general workflow order. Run Section 1 of `1. NEON_metagenomic_environmental_data.R` to generate the sequencing URLs required for `2. Metagenomic_analysis.sh`. Complete the remaining sections of Script 1 to prepare the environmental data and Script 2 to generate the metagenomic data before running `3. Biofilm_gene_analysis.R` for statistical analyses.

---

### 1. NEON_metagenomic_environmental_data.R

**Description:**

This R script retrieves metadata for 146 paired-end NEON metagenomes (292 sequencing URLs) collected from the organic soil horizon in 2017. It then extracts and subsets the sequencing URLs to retain 119 metagenomes (238 URLs) selected for downstream metagenomic processing. It further retrieves and processes discrete and multi-year (2016–2020) environmental data associated with the plots and sites where the metagenomic samples were collected.

**Inputs:**

- `coord_sites.csv`
- NEON data products:
  - `DP1.10107.001`: Soil metagenomics.
  - `DP1.10086.001`: Discrete soil measurements.
  - `DP1.00002.001`: Air temperature.
  - `DP1.00006.001`: Precipitation.
  - `DP1.00041.001`: Soil temperature.
  - `DP1.00094.001`: Soil moisture.

**Outputs:**

- `metagenome_urls_2017.txt`
- `TableS1_metadata2017.csv`

**Execution:**

Run Section 1 before `2. Metagenomic_analysis.sh` to generate the sequencing URLs. Complete the remaining sections before running `3. Biofilm_gene_analysis.R` to prepare the environmental data required for statistical analyses.

---

### 2. Metagenomic_analysis.sh

**Description:**

This Bash pipeline processes 119 paired-end soil metagenomic samples from NEON. It performs quality filtering, contaminant removal, sequence assembly, CDS prediction, functional annotation, read mapping, and quantification of raw read counts assigned to predicted CDSs.

**Inputs:**

- `metagenome_urls_2017.txt`

**Outputs:**

- Summary statistics from the metagenomic processing pipeline presented in `TableS2_summary_pipeline_stats.csv`.
- `all_samples.emapper.annotations.tsv`
- `all_genes_long_all_samples.tsv`

**Execution:**

Run this pipeline after Section 1 of `1. NEON_metagenomic_environmental_data.R` and before `3. Biofilm_gene_analysis.R`.

The processing steps are designed for execution on an HPC cluster using SLURM. Submit each batch script in the indicated order, ensuring that preceding jobs have completed before running dependent steps.

---

### 3. Biofilm_gene_analysis.R

**Description:**

This R script processes raw metagenomic gene abundance data from 112 samples (retained from 119 initially processed) to quantify and summarize the abundance of biofilm-associated genes, both overall and by pathway, across samples. It then evaluates the influence of environmental predictors on biofilm-associated gene abundance at the overall, pathway, and individual KO levels.

**Inputs:**

- `all_samples.emapper.annotations.tsv`
- `all_genes_long_all_samples.tsv`
- `TableS1_metadata2017.csv`
- `TableS2_summary_pipeline_stats.csv`
- `TableS5_KO_list_pathways.csv`

**Outputs:**

- Summary statistics included in `TableS2_summary_pipeline_stats.csv`.
- `TableS6_biofilm_normalized_counts.csv`
- Supplementary Figures S1–S16.
- Figures 1–4.
- Table 1.

**Execution:**

Run this script after completing the metagenomic processing pipeline in `2. Metagenomic_analysis.sh` and retrieving the environmental data using `1. NEON_metagenomic_environmental_data.R`.

All required input files are compiled in the `DATABASES_AND_METADATA` folder and must be available before execution.

Run the sections sequentially, as subsequent analyses depend on objects generated in preceding sections.

