### POTENTIAL FOR BIOFILM FORMATION IS WIDESPREAD AND INFLUENCED BY PRECIPITATION VARIABILITY 
### AND ENVIRONMENTAL EXTREMES ACROSS U.S. SOILS 

# Authors: Eduardo M. Choreño-Parra, Nitish Malhotra, Alejandra Rodríguez-Verdugo, 
# Steven D. Allison, and Kathleen K. Treseder.

# ===================================================================================================================== #

# Description:
# This R script retrieves metadata for 146 paired-end NEON metagenomes (292 sequencing URLs) collected 
# from the organic soil horizon in 2017. It then extracts and subsets the sequencing URLs 
# to retain 119 metagenomes (238 URLs) selected for downstream metagenomic processing. 
# It further retrieves and processes discrete and multi-year environmental data associated with the 
# plots and sites where the metagenomic samples were collected. 

# Inputs:
# - coord_sites.csv: NEON site metadata and coordinates.
# - NEON data products: DP1.10107.001 (soil metagenomics), DP1.10086.001 (discrete soil measurements), 
#   DP1.00002.001 (air temperature), DP1.00006.001 (precipitation), DP1.00041.001 (soil temperature), 
#   and DP1.00094.001 (soil moisture).

# Outputs:
# - metagenome_urls_2017.txt: URLs for downloading the 119 selected paired-end metagenomes (238 URLs) 
# - TableS1.csv: Sample metadata and multi-year summaries of daily environmental conditions 
#   at each site from 2016 to 2020.

# Execution 
# Run Section 1 before “2. Metagenomic_analysis.sh” to generate the sequencing URLs. 
# Complete the remaining sections before running “3. Biofilm_gene_analysis.R” to prepare 
# the environmental data required for statistical analyses.

# Required packages
library(neonUtilities)
library(tidyverse)
library(lubridate)



#### 1. METAGENOMIC DATA ####


##### 1.1. Retrieve NEON metagenomic data #####

metagenomes <- loadByProduct(dpID = 'DP1.10107.001', site="all",
                             startdate = "2017-01", enddate = "2017-12",
                             check.size = FALSE, package = 'expanded',
                             include.provisional = TRUE)


##### 1.2. Extract metadata #####

# Soil metagenomics sequencing metadata
metagenomeSequencingInfo <- metagenomes$mms_metagenomeSequencing
View(metagenomeSequencingInfo)

# Laboratory metadata from soil DNA extraction for metagenomics sequencing
metagenomeDnaExtractionInfo <- metagenomes$mms_metagenomeDnaExtraction
View(metagenomeDnaExtractionInfo)

# Soil metagenomics raw sequence data (includes metadata of the samples)
rawFileInfo <- metagenomes$mms_rawDataFiles
View(rawFileInfo)



##### 1.3. Subset samples from 2017 and organic horizon #####

# Soil metagenomics sequencing metadata (including reads per sample)
metagenomeSequencingInfo2017 <- metagenomeSequencingInfo %>%
  filter(grepl("-O-", dnaSampleID) & grepl("-2017", dnaSampleID))

# Laboratory metadata from soil DNA extraction for metagenomics sequencing
metagenomeDnaExtractionInfo2017 <- metagenomeDnaExtractionInfo %>%
  filter(grepl("-O-", dnaSampleID) & grepl("-2017", dnaSampleID) & grepl("metagenomics",sequenceAnalysisType))

# Raw info for 2017 and organic matter
rawFileInfo2017 <- rawFileInfo %>%
  filter(grepl("-O-", dnaSampleID) & grepl("-2017", dnaSampleID))



##### 1.4. Extract all unique sequencing file URLs #####

urls_2017 <- unique(rawFileInfo2017$rawDataFilePath)


##### 1.5. Subsetting urls #####

# When samples had multiple sequencing runs, we retained the run with the most raw reads

selected_urls <- c(
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellF5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellF5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellF3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellF3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellF2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellF2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellF4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellF4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellF1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellF1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellF6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellF6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellG3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellG3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellF8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellF8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellG1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellG1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellF11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellF11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellG2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellG2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellF12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellF12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellF10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellF10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R1/BMI_HFWM5BGX7_Plate68WellD4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R2/BMI_HFWM5BGX7_Plate68WellD4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellH3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellH3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellF6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellF6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellC6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellC6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellE6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellE6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellH5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellH5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R1/BMI_HT7FTBGX5_Plate35WellD6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HT7FTBGX5_mms_R2/BMI_HT7FTBGX5_Plate35WellD6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellC4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellC4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R2/BMI_HFWM2BGX7_Plate42WellC9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R1/BMI_HFWM2BGX7_Plate42WellC9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R1/BMI_HFWM2BGX7_Plate42WellC10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R2/BMI_HFWM2BGX7_Plate42WellC10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellD1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellD1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R1/BMI_HFWM2BGX7_Plate42WellC11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R2/BMI_HFWM2BGX7_Plate42WellC11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellD8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellD8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellD9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellD9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellD10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellD10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R1/BMI_HVFYCBGX5_Plate42WellE2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVFYCBGX5_mms_R2/BMI_HVFYCBGX5_Plate42WellE2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE7_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE7_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellF2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellF2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellF1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellF1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellF4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellF4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellF3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellF3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellD5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellD5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellD3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellD3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellD1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellD1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellC11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellC11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellD7_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellD7_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellE5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellE5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellD11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellD11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellD9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellD9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R2/BMI_HV3N5BGX5_Plate41WellF6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HV3N5BGX5_mms_R1/BMI_HV3N5BGX5_Plate41WellF6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R2/BMI_HFWM2BGX7_Plate41WellF12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM2BGX7_mms_R1/BMI_HFWM2BGX7_Plate41WellF12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM7C7BGX5_mms_R1/BMI_HM7C7BGX5_Plate25WellH11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM7C7BGX5_mms_R2/BMI_HM7C7BGX5_Plate25WellH11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM7C7BGX5_mms_R2/BMI_HM7C7BGX5_Plate25WellH9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM7C7BGX5_mms_R1/BMI_HM7C7BGX5_Plate25WellH9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R2/BMI_HFWLLBGX7_Plate67WellA3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R1/BMI_HFWLLBGX7_Plate67WellA3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R1/BMI_HFWM5BGX7_Plate67WellC4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R2/BMI_HFWM5BGX7_Plate67WellC4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R1/BMI_HFWLLBGX7_Plate67WellD4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R2/BMI_HFWLLBGX7_Plate67WellD4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate67WellB4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate67WellB4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R2/BMI_HFWLLBGX7_Plate67WellG3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R1/BMI_HFWLLBGX7_Plate67WellG3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R2/BMI_HFWLLBGX7_Plate67WellA4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R1/BMI_HFWLLBGX7_Plate67WellA4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R1/BMI_HFWM5BGX7_Plate67WellE3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R2/BMI_HFWM5BGX7_Plate67WellE3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate67WellB3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate67WellB3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R2/BMI_HFWLLBGX7_Plate67WellD3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLLBGX7_mms_R1/BMI_HFWLLBGX7_Plate67WellD3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R1/BMI_HFWM5BGX7_Plate67WellC3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWM5BGX7_mms_R2/BMI_HFWM5BGX7_Plate67WellC3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellA10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellA10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellH5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellH5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellH2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellH2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellH3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellH3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellH8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellH8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellH6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellH6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellC9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellC9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R1/BMI_HVKLYBGX5_Plate42WellH3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HVKLYBGX5_mms_R2/BMI_HVKLYBGX5_Plate42WellH3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellB11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellB11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R1/BMI_HL7C5BGX7_Plate41WellC5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HL7C5BGX7_mms_R2/BMI_HL7C5BGX7_Plate41WellC5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA7_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA7_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R1/BMI_HFWLVBGX7_Plate61WellA6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HFWLVBGX7_mms_R2/BMI_HFWLVBGX7_Plate61WellA6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellG8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellG8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellG4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellG4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellG6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellG6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellG10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellG10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellG2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellG2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF7_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF7_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellD8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellD8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellD10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellD10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellD12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellD12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellD9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellD9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE4_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE4_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE5_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE5_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellD11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellD11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE9_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE9_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE6_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE6_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF1_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF1_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF3_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF3_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R1/BMI_HM2FCBGX7_Plate61WellF2_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HM2FCBGX7_mms_R2/BMI_HM2FCBGX7_Plate61WellF2_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE12_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE12_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE10_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE10_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE11_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE11_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE8_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE8_R1.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R2/BMI_HLN2JBGX7_Plate61WellE7_R2.fastq.gz",
  "https://storage.neonscience.org/neon-microbial-raw-seq-files/2023/BMI_HLN2JBGX7_mms_R1/BMI_HLN2JBGX7_Plate61WellE7_R1.fastq.gz")

urls_2017_filtered <- selected_urls[selected_urls %in% urls_2017]



##### 1.6. Exporting txt of urls #####

# Export URLs for Bash
writeLines(urls_2017_filtered,"metagenome_urls_2017.txt")



##### 1.7. Merge metagenomic info ######

# rawFileInfo2017 + metagenomeDnaExtractionInfo
rawFileInfo_metagenomeDnaExtractionInfo2017 <- left_join(rawFileInfo2017, metagenomeDnaExtractionInfo,  by = "dnaSampleID")

rawFileInfo_metagenomeDnaExtractionInfo2017 <- rawFileInfo_metagenomeDnaExtractionInfo2017 %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% 
  select_at(vars(-ends_with(".y")))




#### 2. SOIL TEMP AND MOIST (DISCRETE) ####

# This section contains code to download and process discrete measurements of soil temperature and soil moisture
# corresponding to the metagenomic samples.


##### 2.1. Retrieve Soil physical-chemical data #####

phychem2017 <- loadByProduct(dpID = 'DP1.10086.001', site=c("BART","HARV","DSNY","OSBS","STEI","TREE","UNDE","GRSM",
                                                            "ORNL","TALL","WOOD","RMNP","CLBJ","NIWO","ABBY","SJER",
                                                            "BARR","TOOL","DEJU","HEAL"),
                             startdate = "2017-01", enddate = "2017-12",
                             check.size = FALSE, package = 'expanded')

phychem2017_variables <- phychem2017$variables_10086
phychem2017$variables_10086



##### 2.2. Merge soil temperature and soil moisture tables #####

soil_moisture2017 <- phychem2017$sls_soilMoisture                           # Soil moisture multiple samples (O and M) for 2017
sls_soilCoreCollectionInfo <- phychem2017$sls_soilCoreCollection            # Soil temp multiple samples (O and M) and litter depth
sls_metagenomicsPoolingInfo <- phychem2017$sls_metagenomicsPooling          # Info on how to merge moist and temp with metagenomic info
mms_metagenomeDnaExtractionInfo <- metagenomes$mms_metagenomeDnaExtraction  # Metagenomics DNA extraction info


# sls_soilCoreCollectionInfo + soil_moisture2017
merged_soilCoreCollection_moisture2017 <- left_join(sls_soilCoreCollectionInfo, soil_moisture2017,  by = "sampleID")

merged_soilCoreCollection_moisture2017 <- merged_soilCoreCollection_moisture2017 %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% 
  select_at(vars(-ends_with(".y")))

merged_soilCoreCollection_moisture2017_org <- merged_soilCoreCollection_moisture2017 %>%  # Filter only Organic horizon
  filter(grepl("O", horizon))

merged_soilCoreCollection_moisture2017_org
names(merged_soilCoreCollection_moisture2017_org)

# Exclude some variables
var_exclude <- c("subplotID",                     
"biophysicalCriteria",           
"eventID",          
"standingWaterDepth",   
"nTransBoutType",
"samplingImpractical",            
"incubationMethod",               
"incubationCondition",       
"sampleCode",                    
"toxicodendronPossible",          
"horizonDetails",                 
"geneticSampleCode",              
"geneticSampleCondition",         
"geneticSamplePrepMethod",   
"geneticArchiveSample1Code",     
"geneticArchiveSample2Code",   
"geneticArchiveSample3Code",      
"geneticArchiveSample4ID",       
"geneticArchiveSample4Code",     
"geneticArchiveSample5ID",       
"geneticArchiveSample5Code",      
"geneticArchiveSamplePrepMethod", 
"geneticArchiveContainer",       
"biomassID",                     
"biomassCode",                   
"biomassSampleCondition",        
"remarks",   
"dataQF",                        
"sampleCondition",                
"smRemarks",                      
"smDataQF")          

merged_soilCoreCollection_moisture2017_org <- merged_soilCoreCollection_moisture2017_org %>%
  select(-all_of(var_exclude))



##### 2.3. Average across samples per plot ######

merged_soilCoreCollection_moisture2017_org_avg <- merged_soilCoreCollection_moisture2017_org %>%
  group_by(plotID) %>%                                    # Group by plotID
  mutate(
    soilMoisture_avg = mean(soilMoisture, na.rm = TRUE),  # Average soilMoisture from different samples per plot
    soilTemp_avg = mean(soilTemp, na.rm = TRUE),          # Average soilTemperature from different samples per plot
    across( where(~ all(!is.na(.x))),                     # For all non-numeric consistent variables
      ~ first(.x))) %>%                                   # Keep the first value
  ungroup()  

# Mean depth chemistry
mean(merged_soilCoreCollection_moisture2017_org_avg$sampleBottomDepth)
merged_soilCoreCollection_moisture2017_org_avg$sampleBottomDepth




#### 3. LINK METAGENOMIC DATA + SOIL TEMP AND MOIST ####


##### 3.1. Link metagenomic metadata with soil sampling metadata #### 

# Merge metagenomic metadata with pooling information:
# rawFileInfo_metagenomeDnaExtractionInfo2017 + sls_metagenomicsPoolingInfo
key_to_merge_tables2017 <- left_join(rawFileInfo_metagenomeDnaExtractionInfo2017, 
                                     sls_metagenomicsPoolingInfo,  
                                     by = "genomicsSampleID")

key_to_merge_tables2017 <- key_to_merge_tables2017 %>% 
  rename_at(vars(ends_with(".x")), ~str_replace(., "\\..$","")) %>% 
  select_at(vars(-ends_with(".y")))


##### 3.2. Prepare plot-level soil temperature and moisture data #####

# Retain one row per plot with average soil moisture and temperature
aggregated_soilCoreCollection_moisture <- merged_soilCoreCollection_moisture2017_org_avg %>%
  group_by(plotID) %>%
  summarise(
    nlcdClass = first(nlcdClass),                
    soilMoisture_avg = mean(soilMoisture_avg),   # Keep the average value
    soilTemp_avg = mean(soilTemp_avg),           # Keep the average value
    plotType = first(plotType))                 

View(aggregated_soilCoreCollection_moisture)
View(merged_soilCoreCollection_moisture2017_org_avg)




##### 3.3. Merge metagenomic and environmental data ####

# Assign plot-level environmental (soil temp and moist per plot) measurements to metagenomic samples:
# key_to_merge_tables2017 + merged_soilCoreCollection_moisture2017_org_avg$plotID

final_table_2017 <- left_join(key_to_merge_tables2017,aggregated_soilCoreCollection_moisture, by = "plotID")

nrow(final_table_2017)
View(final_table_2017)


##### 3.4. Keeping data from selected samples @ ######

# Keep only samples with URLs included in selected_urls
final_table_2017_dedup <- final_table_2017 %>%
  filter(rawDataFilePath %in% selected_urls)

View(final_table_2017_dedup)

# Confirm number of samples
nrow(final_table_2017_dedup)
# Confirm all retained URLs are in selected_urls
all(final_table_2017_dedup$rawDataFilePath %in% selected_urls)




#### 4. AIR TEMPERATURE (Single) ####

saat <- loadByProduct(dpID="DP1.00002.001", site=c("BART","HARV","DSNY","OSBS","STEI","TREE","UNDE","GRSM",
                                                   "ORNL","TALL","WOOD","RMNP","CLBJ","NIWO","ABBY","SJER",
                                                   "BARR","TOOL","DEJU","HEAL"), 
                      startdate="2016-01", enddate="2020-12", 
                      package="basic", timeIndex="30",
                      check.size = F)

SAAT_30min <- saat$SAAT_30min

# Cleaning NAs and 1 (bad quality data)
SAAT_30min_clean <- SAAT_30min %>%
  filter(finalQF != 1 & !is.na(tempSingleMean))

# Convert to date
SAAT_30min_clean$Date <- as.Date(SAAT_30min_clean$startDateTime)
SAAT_30min_clean$Date



##### 4.1. Defining  valid months-years #####

SAAT_30min_clean <- SAAT_30min_clean %>%
  mutate(Year = year(Date),
         Month = month(Date))

# 1) Count distinct dates per site-year-month
month_counts <- SAAT_30min_clean %>%
  group_by(siteID, Year, Month) %>%
  summarise(n_days = n_distinct(Date), .groups = "drop")

# 2) Valid months: >= 10 days
valid_months <- month_counts %>%
  filter(n_days >= 10)

# 3) Valid years: >= 11 valid months
valid_years <- valid_months %>%
  group_by(siteID, Year) %>%
  summarise(n_valid_months = n(), .groups = "drop") %>%
  filter(n_valid_months >= 11)

# 4) Final filtered 30-min data: only valid months inside valid years
SAAT_30min_full <- SAAT_30min_clean %>%
  semi_join(valid_years,  by = c("siteID", "Year")) %>%
  semi_join(valid_months, by = c("siteID", "Year", "Month"))



##### 4.2. Multi-year summary (min, max, mean, cv) ##### 

temp_day <- SAAT_30min_full %>%
  group_by(domainID, siteID, Date) %>%  # Group by both site and date
  summarise(dayMax = max(tempSingleMaximum, na.rm = TRUE), .groups = "drop", # Calculate daily values in each site 
            dayMean = mean(tempSingleMean, na.rm = TRUE),
            dayMin = min(tempSingleMinimum)) %>%    
  group_by(siteID) %>%                                                       # Multi-year summary
  mutate(avg_dayMax = mean(dayMax, na.rm = TRUE)) %>%  
  group_by(siteID) %>%  
  mutate(avg_dayMean = mean(dayMean, na.rm = TRUE)) %>% 
  group_by(siteID) %>%  
  mutate(avg_dayMin = mean(dayMin, na.rm = TRUE)) %>%
  group_by(siteID) %>%
  mutate(
    sd_dayMax   = sd(dayMax,  na.rm = TRUE),
    sd_dayMean  = sd(dayMean, na.rm = TRUE),
    sd_dayMin   = sd(dayMin,  na.rm = TRUE),
    cv_dayMax   = 100 * sd_dayMax  / avg_dayMax,
    cv_dayMean  = 100 * sd_dayMean / avg_dayMean,
    cv_dayMin   = 100 * sd_dayMin  / avg_dayMin)
  


##### 4.3. Seasonal summaries #####

# Create season subsets from temp_day
# Spring begins March 20 (equinox)
# Summer begins June 21 (solstice)
# Autumn begins September 22 (equinox)
# Winter begins December 21 (solstice)

### Subset 

# Winter (Dec 21 – Mar 19)
SAAT_30min_full_winter <- SAAT_30min_full %>%
  filter(
    (month(Date) == 12 & day(Date) >= 21) |
      (month(Date) %in% c(1, 2)) |
      (month(Date) == 3 & day(Date) < 20))

# Spring (Mar 20 – Jun 20)
SAAT_30min_full_spring <- SAAT_30min_full %>%
  filter(
    (month(Date) == 3 & day(Date) >= 20) |
      (month(Date) %in% c(4, 5)) |
      (month(Date) == 6 & day(Date) < 21))

# Summer (Jun 21 – Sep 21)
SAAT_30min_full_summer <- SAAT_30min_full %>%
  filter(
    (month(Date) == 6 & day(Date) >= 21) |
      (month(Date) %in% c(7, 8)) |
      (month(Date) == 9 & day(Date) < 22))

# Fall (Sep 22 – Dec 20)
SAAT_30min_full_fall <- SAAT_30min_full %>%
  filter(
    (month(Date) == 9 & day(Date) >= 22) |
      (month(Date) %in% c(10, 11)) |
      (month(Date) == 12 & day(Date) < 21))


### Summary per season 

###### Winter ######
temp_day_winter <- SAAT_30min_full_winter %>%
  group_by(domainID, siteID, Date) %>%  # Group by both site and date
  summarise(dayMax = max(tempSingleMaximum, na.rm = TRUE), .groups = "drop", # Calculate max temp at each day (daily max temp)
            dayMean = mean(tempSingleMean, na.rm = TRUE),
            dayMin = min(tempSingleMinimum)) %>%    # Calculate min temp at each day
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMax_winter = mean(dayMax, na.rm = TRUE)) %>%  
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMean_winter = mean(dayMean, na.rm = TRUE)) %>% 
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMin_winter = mean(dayMin, na.rm = TRUE))

###### Spring ######
temp_day_spring <- SAAT_30min_full_spring %>%
  group_by(domainID, siteID, Date) %>%  # Group by both site and date
  summarise(dayMax = max(tempSingleMaximum, na.rm = TRUE), .groups = "drop", # Calculate max temp at each day (daily max temp)
            dayMean = mean(tempSingleMean, na.rm = TRUE),
            dayMin = min(tempSingleMinimum)) %>%    # Calculate min temp at each day
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMax_spring = mean(dayMax, na.rm = TRUE)) %>%  
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMean_spring = mean(dayMean, na.rm = TRUE)) %>% 
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMin_spring = mean(dayMin, na.rm = TRUE))

###### Summer ######
temp_day_summer <- SAAT_30min_full_summer %>%
  group_by(domainID, siteID, Date) %>%  # Group by both site and date
  summarise(dayMax = max(tempSingleMaximum, na.rm = TRUE), .groups = "drop", # Calculate max temp at each day (daily max temp)
            dayMean = mean(tempSingleMean, na.rm = TRUE),
            dayMin = min(tempSingleMinimum)) %>%    # Calculate min temp at each day
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMax_summer = mean(dayMax, na.rm = TRUE)) %>%  
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMean_summer = mean(dayMean, na.rm = TRUE)) %>% 
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMin_summer = mean(dayMin, na.rm = TRUE))

###### Fall ######
temp_day_fall <- SAAT_30min_full_fall %>%
  group_by(domainID, siteID, Date) %>%  # Group by both site and date
  summarise(dayMax = max(tempSingleMaximum, na.rm = TRUE), .groups = "drop", # Calculate max temp at each day (daily max temp)
            dayMean = mean(tempSingleMean, na.rm = TRUE),
            dayMin = min(tempSingleMinimum)) %>%    # Calculate min temp at each day
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMax_fall = mean(dayMax, na.rm = TRUE)) %>%  
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMean_fall = mean(dayMean, na.rm = TRUE)) %>% 
  group_by(siteID) %>%  # Group by siteID for the average
  mutate(avg_dayMin_fall = mean(dayMin, na.rm = TRUE))




#### 5. PRECIPITATION ####


precipitation <- loadByProduct(dpID="DP1.00006.001", site=c("BART","HARV","DSNY","OSBS","STEI","TREE","UNDE","GRSM",
                                                   "ORNL","TALL","WOOD","RMNP","CLBJ","NIWO","ABBY","SJER",
                                                   "BARR","TOOL","DEJU","HEAL"), 
                      startdate="2016-01", enddate="2020-12", 
                      package="basic", timeIndex="30",
                      check.size = F)

# Working with secondary precipitation
SECPRE_30min <- precipitation$SECPRE_30min

# Cleaning NAs and 1
SECPRE_30min_clean <- SECPRE_30min %>%
  filter(secPrecipRangeQF != 1 & secPrecipRangeQF != -1 & !is.na(secPrecipBulk))

# Convert to date
SECPRE_30min_clean$Date <- as.Date(SECPRE_30min_clean$startDateTime)



##### 5.1.Defining valid months-years #####

SECPRE_30min_clean <- SECPRE_30min_clean %>%
  mutate(Year = year(Date),
         Month = month(Date))

# 1) Count distinct dates per site-year-month
sec_month_counts <- SECPRE_30min_clean %>%
  group_by(siteID, Year, Month) %>%
  summarise(n_days = n_distinct(Date), .groups = "drop")

# 2) Valid months: >= 10 days
sec_valid_months <- sec_month_counts %>%
  filter(n_days >= 10)

# 3) Valid years: >= 11 valid months
sec_valid_years <- sec_valid_months %>%
  group_by(siteID, Year) %>%
  summarise(n_valid_months = n(), .groups = "drop") %>%
  filter(n_valid_months >= 11)

# 4) Final filtered 30-min data: only valid months inside valid years
SECPRE_30min_full <- SECPRE_30min_clean %>%
  semi_join(sec_valid_years,  by = c("siteID", "Year")) %>%
  semi_join(sec_valid_months, by = c("siteID", "Year", "Month"))


##### 5.2. Multi-year summary (mean, cv) #####

tot_secprecip_day <- SECPRE_30min_full %>%
  group_by(siteID, Date) %>%
  summarise(dayTot = sum(secPrecipBulk, na.rm = TRUE), .groups = "drop") %>%  # Calculate daily values in each site
  group_by(siteID) %>%                                                        # Multi-year summary
  mutate(
    avg_dayTot = mean(dayTot, na.rm = TRUE),
    sd_dayTot  = sd(dayTot,   na.rm = TRUE),
    cv_dayTot  = 100 * sd_dayTot / avg_dayTot)



##### 5.3 Seasonal summaries #####

# Create season subsets
# Spring begins March 20 (equinox)
# Summer begins June 21 (solstice)
# Autumn begins September 22 (equinox)
# Winter begins December 21 (solstice)

### Subset 

# Winter (Dec 21 – Mar 19)
SECPRE_30min_full_winter <- SECPRE_30min_full %>%
  filter(
    (month(Date) == 12 & day(Date) >= 21) |
      (month(Date) %in% c(1, 2)) |
      (month(Date) == 3 & day(Date) < 20))

# Spring (Mar 20 – Jun 20)
SECPRE_30min_full_spring <- SECPRE_30min_full %>%
  filter(
    (month(Date) == 3 & day(Date) >= 20) |
      (month(Date) %in% c(4, 5)) |
      (month(Date) == 6 & day(Date) < 21))

# Summer (Jun 21 – Sep 21)
SECPRE_30min_full_summer <- SECPRE_30min_full %>%
  filter(
    (month(Date) == 6 & day(Date) >= 21) |
      (month(Date) %in% c(7, 8)) |
      (month(Date) == 9 & day(Date) < 22))

# Fall (Sep 22 – Dec 20)
SECPRE_30min_full_fall <- SECPRE_30min_full %>%
  filter(
    (month(Date) == 9 & day(Date) >= 22) |
      (month(Date) %in% c(10, 11)) |
      (month(Date) == 12 & day(Date) < 21))


### Summary per season

###### Winter ######
tot_secprecip_day_winter <- SECPRE_30min_full_winter %>%
  group_by(siteID, Date) %>%
  summarise(dayTot = sum(secPrecipBulk, na.rm = TRUE), .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(avg_dayTot_winter = mean(dayTot, na.rm = TRUE))

###### Spring ######
tot_secprecip_day_spring <- SECPRE_30min_full_spring %>%
  group_by(siteID, Date) %>%
  summarise(dayTot = sum(secPrecipBulk, na.rm = TRUE), .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(avg_dayTot_spring = mean(dayTot, na.rm = TRUE))

###### Summer ######
tot_secprecip_day_summer <- SECPRE_30min_full_summer %>%
  group_by(siteID, Date) %>%
  summarise(dayTot = sum(secPrecipBulk, na.rm = TRUE), .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(avg_dayTot_summer = mean(dayTot, na.rm = TRUE))

###### Fall ######
tot_secprecip_day_fall <- SECPRE_30min_full_fall %>%
  group_by(siteID, Date) %>%
  summarise(dayTot = sum(secPrecipBulk, na.rm = TRUE), .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(avg_dayTot_fall = mean(dayTot, na.rm = TRUE))




#### 6. SOIL TEMPERATURE PROFILE ####

soil_temp <- loadByProduct(dpID="DP1.00041.001", site=c("BART","HARV","DSNY","OSBS","STEI","TREE","UNDE","GRSM",
                                                            "ORNL","TALL","WOOD","RMNP","CLBJ","NIWO","ABBY","SJER",
                                                            "BARR","TOOL","DEJU","HEAL"), 
                               startdate="2016-01", enddate="2020-12", 
                               package="basic", timeIndex="30",
                               check.size = F, nCores=4)

ST_30_minute <- soil_temp$ST_30_minute

# Cleaning NAs and 1
ST_30_minute_clean <- ST_30_minute %>%
  filter(finalQF != 1 & !is.na(soilTempMaximum))

# Only shallow soil depth (level 501)
ST_30_minute_clean_501 <- ST_30_minute_clean %>%
  filter(verticalPosition == 501)

# Convert to date
ST_30_minute_clean_501$Date <- as.Date(ST_30_minute_clean_501$startDateTime)


###### 6.1. Defining valid months-years ######

ST_30_minute_clean_501 <- ST_30_minute_clean_501 %>%
  mutate(Year = year(Date),
         Month = month(Date))

# 1) Count distinct dates per site-year-month
st_month_counts <- ST_30_minute_clean_501 %>%
  group_by(siteID, Year, Month) %>%
  summarise(n_days = n_distinct(Date), .groups = "drop")

# 2) Valid months: >= 10 days
st_valid_months <- st_month_counts %>%
  filter(n_days >= 10)

# 3) Valid years: >= 11 valid months
st_valid_years <- st_valid_months %>%
  group_by(siteID, Year) %>%
  summarise(n_valid_months = n(), .groups = "drop") %>%
  filter(n_valid_months >= 11)

# 4) Final filtered 30-min data: only valid months inside valid years
ST_30min_full_501 <- ST_30_minute_clean_501 %>%
  semi_join(st_valid_years,  by = c("siteID", "Year")) %>%
  semi_join(st_valid_months, by = c("siteID", "Year", "Month"))



##### 6.2. Multi-year summary (min, max, mean, cv) #####

ST_day_by_sensor_501 <- ST_30min_full_501 %>%
  group_by(siteID, Date, horizontalPosition) %>%     # Calculate daily values for each sensor within each site
  summarise(
    dayMax  = max(soilTempMaximum, na.rm = TRUE),
    dayMean = mean(soilTempMean,   na.rm = TRUE),
    dayMin  = min(soilTempMinimum, na.rm = TRUE),
    .groups = "drop")

max_ST_day_501 <- ST_day_by_sensor_501 %>%
  group_by(siteID, Date) %>%                         # Average daily values across sensor positions within each site and date            
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  group_by(siteID) %>%                               # Multi-year summary
  mutate(
    avg_dayMax  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean = mean(dayMean, na.rm = TRUE),
    avg_dayMin  = mean(dayMin,  na.rm = TRUE),
    sd_dayMax   = sd(dayMax,  na.rm = TRUE),
    sd_dayMean  = sd(dayMean, na.rm = TRUE),
    sd_dayMin   = sd(dayMin,  na.rm = TRUE),
    cv_dayMax   = 100 * sd_dayMax  / avg_dayMax,
    cv_dayMean  = 100 * sd_dayMean / avg_dayMean,
    cv_dayMin   = 100 * sd_dayMin  / avg_dayMin)



##### 6.3. Seasonal summaries #####

# Create season subsets
# Spring begins March 20 (equinox)
# Summer begins June 21 (solstice)
# Autumn begins September 22 (equinox)
# Winter begins December 21 (solstice)

### Subset 

# Winter (Dec 21 – Mar 19)
ST_30min_full_501_winter <- ST_30min_full_501 %>%
  filter(
    (month(Date) == 12 & day(Date) >= 21) |
      (month(Date) %in% c(1, 2)) |
      (month(Date) == 3 & day(Date) < 20))

# Spring (Mar 20 – Jun 20)
ST_30min_full_501_spring <- ST_30min_full_501 %>%
  filter(
    (month(Date) == 3 & day(Date) >= 20) |
      (month(Date) %in% c(4, 5)) |
      (month(Date) == 6 & day(Date) < 21))

# Summer (Jun 21 – Sep 21)
ST_30min_full_501_summer <- ST_30min_full_501 %>%
  filter(
    (month(Date) == 6 & day(Date) >= 21) |
      (month(Date) %in% c(7, 8)) |
      (month(Date) == 9 & day(Date) < 22))

# Fall (Sep 22 – Dec 20)
ST_30min_full_501_fall <- ST_30min_full_501 %>%
  filter(
    (month(Date) == 9 & day(Date) >= 22) |
      (month(Date) %in% c(10, 11)) |
      (month(Date) == 12 & day(Date) < 21))


### Summary per season (daily metrics + site averages)

###### Winter ######
ST_day_by_sensor_501_winter <- ST_30min_full_501_winter %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(soilTempMaximum, na.rm = TRUE),
    dayMean = mean(soilTempMean,   na.rm = TRUE),
    dayMin  = min(soilTempMinimum, na.rm = TRUE),
    .groups = "drop")

max_ST_day_501_winter <- ST_day_by_sensor_501_winter %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_winter  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_winter = mean(dayMean, na.rm = TRUE),
    avg_dayMin_winter  = mean(dayMin,  na.rm = TRUE))


####### Spring ######
ST_day_by_sensor_501_spring <- ST_30min_full_501_spring %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(soilTempMaximum, na.rm = TRUE),
    dayMean = mean(soilTempMean,   na.rm = TRUE),
    dayMin  = min(soilTempMinimum, na.rm = TRUE),
    .groups = "drop")

max_ST_day_501_spring <- ST_day_by_sensor_501_spring %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_spring  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_spring = mean(dayMean, na.rm = TRUE),
    avg_dayMin_spring  = mean(dayMin,  na.rm = TRUE))


###### Summer ###### 
ST_day_by_sensor_501_summer <- ST_30min_full_501_summer %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(soilTempMaximum, na.rm = TRUE),
    dayMean = mean(soilTempMean,   na.rm = TRUE),
    dayMin  = min(soilTempMinimum, na.rm = TRUE),
    .groups = "drop")

max_ST_day_501_summer <- ST_day_by_sensor_501_summer %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_summer  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_summer = mean(dayMean, na.rm = TRUE),
    avg_dayMin_summer  = mean(dayMin,  na.rm = TRUE))


###### Fall ######
ST_day_by_sensor_501_fall <- ST_30min_full_501_fall %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(soilTempMaximum, na.rm = TRUE),
    dayMean = mean(soilTempMean,   na.rm = TRUE),
    dayMin  = min(soilTempMinimum, na.rm = TRUE),
    .groups = "drop")

max_ST_day_501_fall <- ST_day_by_sensor_501_fall %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_fall  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_fall = mean(dayMean, na.rm = TRUE),
    avg_dayMin_fall  = mean(dayMin,  na.rm = TRUE))



#### 7. SOIL MOISTURE PROFILE ####

# Download in chunks
soil_moist_1 <- loadByProduct(
  dpID="DP1.00094.001",
  site=c("BART","HARV","DSNY","OSBS","STEI"),
  startdate="2016-01", enddate="2020-12",
  package="basic", timeIndex="30",
  check.size=FALSE,
  nCores=1)

list2env(soil_moist_1, .GlobalEnv)
SWS_30_minute_1 <- SWS_30_minute
rm(SWS_30_minute)
gc()


soil_moist_2 <- loadByProduct(
  dpID="DP1.00094.001",
  site=c("TREE","UNDE","GRSM","ORNL","TALL"),
  startdate="2016-01", enddate="2020-12",
  package="basic", timeIndex="30",
  check.size=FALSE,
  nCores=1)

list2env(soil_moist_2, .GlobalEnv)
SWS_30_minute_2 <- SWS_30_minute
rm(SWS_30_minute)
gc()

soil_moist_3 <- loadByProduct(
  dpID="DP1.00094.001",
  site=c("WOOD","RMNP","CLBJ","NIWO","ABBY"),
  startdate="2016-01", enddate="2020-12",
  package="basic", timeIndex="30",
  check.size=FALSE,
  nCores=1)

list2env(soil_moist_3, .GlobalEnv)
SWS_30_minute_3 <- SWS_30_minute
rm(SWS_30_minute)
gc()


soil_moist_4 <- loadByProduct(
  dpID="DP1.00094.001",
  site=c("SJER","BARR","TOOL","DEJU","HEAL"),
  startdate="2016-01", enddate="2020-12",
  package="basic", timeIndex="30",
  check.size=FALSE,
  nCores=1)

list2env(soil_moist_4, .GlobalEnv)
SWS_30_minute_4 <- SWS_30_minute
rm(SWS_30_minute)
gc()

# Merge
SWS_30_minute_all <- bind_rows(SWS_30_minute_1, SWS_30_minute_2,
                               SWS_30_minute_3, SWS_30_minute_4)


# Cleaning NAs and 1
SWS_30_minute_clean <- SWS_30_minute_all %>%
  filter(VSWCFinalQF != 1 & !is.na(VSWCMaximum))

# Only shallow soil depth (level 501)
SWS_30_minute_clean_501 <- SWS_30_minute_clean %>%
  filter(verticalPosition == 501)

# convert to date
SWS_30_minute_clean_501$Date <- as.Date(SWS_30_minute_clean_501$startDateTime)


###### 7.1. Defining valid months-years ######

SWS_30_minute_clean_501 <- SWS_30_minute_clean_501 %>%
  mutate(
    Year  = year(Date),
    Month = month(Date))

# 1) Count distinct dates per site-year-month
swc_month_counts <- SWS_30_minute_clean_501 %>%
  group_by(siteID, Year, Month) %>%
  summarise(n_days = n_distinct(Date), .groups = "drop")

# 2) Valid months: >= 10 days
swc_valid_months <- swc_month_counts %>%
  filter(n_days >= 10)

# 3) Valid years: >= 11 valid months
swc_valid_years <- swc_valid_months %>%
  group_by(siteID, Year) %>%
  summarise(n_valid_months = n(), .groups = "drop") %>%
  filter(n_valid_months >= 11)

# 4) Final filtered 30-min data: only valid months inside valid years
SWS_30min_full_501 <- SWS_30_minute_clean_501 %>%
  semi_join(swc_valid_years,  by = c("siteID", "Year")) %>%
  semi_join(swc_valid_months, by = c("siteID", "Year", "Month"))


##### 7.2. Multi-year summary (min, max, mean, cv) #####

swc_day_by_sensor_501 <- SWS_30min_full_501 %>%
  group_by(siteID, Date, horizontalPosition) %>%    # Calculate daily values for each sensor within each site
  summarise(
    dayMax  = max(VSWCMaximum, na.rm = TRUE),
    dayMean = mean(VSWCMean,   na.rm = TRUE),
    dayMin  = min(VSWCMinimum, na.rm = TRUE),
    .groups = "drop")

swc_day_501 <- swc_day_by_sensor_501 %>%
  group_by(siteID, Date) %>%                        # Average daily values across sensor positions within each site and date
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    dayMax  = dayMax  * 100,
    dayMean = dayMean * 100,
    dayMin  = dayMin  * 100) %>%
  group_by(siteID) %>%                              # Multi-year summary  
  mutate(
    avg_dayMax  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean = mean(dayMean, na.rm = TRUE),
    avg_dayMin  = mean(dayMin,  na.rm = TRUE),
    sd_dayMax   = sd(dayMax,  na.rm = TRUE),
    sd_dayMean  = sd(dayMean, na.rm = TRUE),
    sd_dayMin   = sd(dayMin,  na.rm = TRUE),
    cv_dayMax   = 100 * sd_dayMax  / avg_dayMax,
    cv_dayMean  = 100 * sd_dayMean / avg_dayMean,
    cv_dayMin   = 100 * sd_dayMin  / avg_dayMin)



##### 7.3. Seasonal summaries #####

# Create season subsets
# Spring begins March 20 (equinox)
# Summer begins June 21 (solstice)
# Autumn begins September 22 (equinox)
# Winter begins December 21 (solstice)

### Subset 

# Winter (Dec 21 – Mar 19)
SWS_30min_full_501_winter <- SWS_30min_full_501 %>%
  filter(
    (month(Date) == 12 & day(Date) >= 21) |
      (month(Date) %in% c(1, 2)) |
      (month(Date) == 3 & day(Date) < 20))

# Spring (Mar 20 – Jun 20)
SWS_30min_full_501_spring <- SWS_30min_full_501 %>%
  filter(
    (month(Date) == 3 & day(Date) >= 20) |
      (month(Date) %in% c(4, 5)) |
      (month(Date) == 6 & day(Date) < 21))

# Summer (Jun 21 – Sep 21)
SWS_30min_full_501_summer <- SWS_30min_full_501 %>%
  filter(
    (month(Date) == 6 & day(Date) >= 21) |
      (month(Date) %in% c(7, 8)) |
      (month(Date) == 9 & day(Date) < 22))

# Fall (Sep 22 – Dec 20)
SWS_30min_full_501_fall <- SWS_30min_full_501 %>%
  filter(
    (month(Date) == 9 & day(Date) >= 22) |
      (month(Date) %in% c(10, 11)) |
      (month(Date) == 12 & day(Date) < 21))


### Summary per season (daily metrics + site averages)

###### Winter ###### 
swc_day_by_sensor_501_winter <- SWS_30min_full_501_winter %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(VSWCMaximum, na.rm = TRUE),
    dayMean = mean(VSWCMean,   na.rm = TRUE),
    dayMin  = min(VSWCMinimum, na.rm = TRUE),
    .groups = "drop")

swc_day_501_winter <- swc_day_by_sensor_501_winter %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    dayMax  = dayMax  * 100,
    dayMean = dayMean * 100,
    dayMin  = dayMin  * 100) %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_winter  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_winter = mean(dayMean, na.rm = TRUE),
    avg_dayMin_winter  = mean(dayMin,  na.rm = TRUE))


###### Spring ###### 
swc_day_by_sensor_501_spring <- SWS_30min_full_501_spring %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(VSWCMaximum, na.rm = TRUE),
    dayMean = mean(VSWCMean,   na.rm = TRUE),
    dayMin  = min(VSWCMinimum, na.rm = TRUE),
    .groups = "drop")

swc_day_501_spring <- swc_day_by_sensor_501_spring %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    dayMax  = dayMax  * 100,
    dayMean = dayMean * 100,
    dayMin  = dayMin  * 100) %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_spring  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_spring = mean(dayMean, na.rm = TRUE),
    avg_dayMin_spring  = mean(dayMin,  na.rm = TRUE))


###### Summer ###### 
swc_day_by_sensor_501_summer <- SWS_30min_full_501_summer %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(VSWCMaximum, na.rm = TRUE),
    dayMean = mean(VSWCMean,   na.rm = TRUE),
    dayMin  = min(VSWCMinimum, na.rm = TRUE),
    .groups = "drop")

swc_day_501_summer <- swc_day_by_sensor_501_summer %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    dayMax  = dayMax  * 100,
    dayMean = dayMean * 100,
    dayMin  = dayMin  * 100) %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_summer  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_summer = mean(dayMean, na.rm = TRUE),
    avg_dayMin_summer  = mean(dayMin,  na.rm = TRUE))


###### Fall ######
swc_day_by_sensor_501_fall <- SWS_30min_full_501_fall %>%
  group_by(siteID, Date, horizontalPosition) %>%
  summarise(
    dayMax  = max(VSWCMaximum, na.rm = TRUE),
    dayMean = mean(VSWCMean,   na.rm = TRUE),
    dayMin  = min(VSWCMinimum, na.rm = TRUE),
    .groups = "drop")

swc_day_501_fall <- swc_day_by_sensor_501_fall %>%
  group_by(siteID, Date) %>%
  summarise(
    dayMax  = mean(dayMax,  na.rm = TRUE),
    dayMean = mean(dayMean, na.rm = TRUE),
    dayMin  = mean(dayMin,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    dayMax  = dayMax  * 100,
    dayMean = dayMean * 100,
    dayMin  = dayMin  * 100) %>%
  group_by(siteID) %>%
  mutate(
    avg_dayMax_fall  = mean(dayMax,  na.rm = TRUE),
    avg_dayMean_fall = mean(dayMean, na.rm = TRUE),
    avg_dayMin_fall  = mean(dayMin,  na.rm = TRUE))







#### 8. PREPARING SUMMARY TABLES ####


# Soil Temperature Profile 
max_ST_day_501
max_ST_day_501_winter
max_ST_day_501_spring
max_ST_day_501_summer
max_ST_day_501_fall

# Secondary precipitation 
tot_secprecip_day
tot_secprecip_day_winter
tot_secprecip_day_spring
tot_secprecip_day_summer
tot_secprecip_day_fall

# Air Temperature (single aspired)
temp_day
temp_day_winter
temp_day_spring
temp_day_summer
temp_day_fall

# Soil Moisture profile
swc_day_501
swc_day_501_winter
swc_day_501_spring
swc_day_501_summer
swc_day_501_fall



##### 1. Soil Temperature #####

soil_temp_summary <- max_ST_day_501 %>%
  select(siteID, starts_with("avg_"), starts_with("cv_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilTemp_", .),
              c(starts_with("avg_"), starts_with("cv_")))

soil_temp_summary_winter <- max_ST_day_501_winter %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilTemp_", .), starts_with("avg_"))

soil_temp_summary_spring <- max_ST_day_501_spring %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilTemp_", .), starts_with("avg_"))

soil_temp_summary_summer <- max_ST_day_501_summer %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilTemp_", .), starts_with("avg_"))

soil_temp_summary_fall <- max_ST_day_501_fall %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilTemp_", .), starts_with("avg_"))



##### 2. Secondary Precipitation  #####

secPrecip_summary <- tot_secprecip_day %>%
  select(siteID, avg_dayTot, cv_dayTot) %>%
  distinct() %>%
  rename(secPrecip_avg_dayTot = avg_dayTot,
         secPrecip_cv_dayTot = cv_dayTot)
secPrecip_summary

secPrecip_summary_winter <- tot_secprecip_day_winter %>%
  select(siteID, avg_dayTot_winter) %>%
  distinct() %>%
  rename(secPrecip_avg_dayTot_winter = avg_dayTot_winter)

secPrecip_summary_spring <- tot_secprecip_day_spring %>%
  select(siteID, avg_dayTot_spring) %>%
  distinct() %>%
  rename(secPrecip_avg_dayTot_spring = avg_dayTot_spring)

secPrecip_summary_summer <- tot_secprecip_day_summer %>%
  select(siteID, avg_dayTot_summer) %>%
  distinct() %>%
  rename(secPrecip_avg_dayTot_summer = avg_dayTot_summer)

secPrecip_summary_fall <- tot_secprecip_day_fall %>%
  select(siteID, avg_dayTot_fall) %>%
  distinct() %>%
  rename(secPrecip_avg_dayTot_fall = avg_dayTot_fall)



##### 3. Air Temperature ####

air_temp_summary <- temp_day %>%
  select(siteID, starts_with("avg_"), starts_with("cv_") ) %>%
  distinct() %>%
  rename_with(~ paste0("airTemp_", .), 
              c(starts_with("avg_"), starts_with("cv_")))

air_temp_summary_winter <- temp_day_winter %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("airTemp_", .), starts_with("avg_"))

air_temp_summary_spring <- temp_day_spring %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("airTemp_", .), starts_with("avg_"))

air_temp_summary_summer <- temp_day_summer %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("airTemp_", .), starts_with("avg_"))

air_temp_summary_fall <- temp_day_fall %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("airTemp_", .), starts_with("avg_"))


##### 4. Soil Moisture profile #####
soil_moist_summary <- swc_day_501 %>%
  select(siteID, starts_with("avg_"), starts_with("cv_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilMoist_", .), 
              c(starts_with("avg_"), starts_with("cv_")))
soil_moist_summary

soil_moist_summary_winter <- swc_day_501_winter %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilMoist_", .), starts_with("avg_"))
soil_moist_summary_winter

soil_moist_summary_spring <- swc_day_501_spring %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilMoist_", .), starts_with("avg_"))

soil_moist_summary_summer <- swc_day_501_summer %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilMoist_", .), starts_with("avg_"))

soil_moist_summary_fall <- swc_day_501_fall %>%
  select(siteID, starts_with("avg_")) %>%
  distinct() %>%
  rename_with(~ paste0("soilMoist_", .), starts_with("avg_"))



##### Summary environmental variables @ ####

envir_variables_summary <- soil_temp_summary %>%
  full_join(soil_temp_summary_winter, by = "siteID") %>%
  full_join(soil_temp_summary_spring, by = "siteID") %>%
  full_join(soil_temp_summary_summer, by = "siteID") %>%
  full_join(soil_temp_summary_fall, by = "siteID") %>%
  
  full_join(soil_moist_summary, by = "siteID") %>%
  full_join(soil_moist_summary_winter, by = "siteID") %>%
  full_join(soil_moist_summary_spring, by = "siteID") %>%
  full_join(soil_moist_summary_summer, by = "siteID") %>%
  full_join(soil_moist_summary_fall, by = "siteID") %>%
  
  full_join(air_temp_summary, by = "siteID") %>%
  full_join(air_temp_summary_winter, by = "siteID") %>%
  full_join(air_temp_summary_spring, by = "siteID") %>%
  full_join(air_temp_summary_summer, by = "siteID") %>%
  full_join(air_temp_summary_fall, by = "siteID") %>%
  
  full_join(secPrecip_summary, by = "siteID") %>%
  full_join(secPrecip_summary_winter, by = "siteID") %>%
  full_join(secPrecip_summary_spring, by = "siteID") %>%
  full_join(secPrecip_summary_summer, by = "siteID") %>%
  full_join(secPrecip_summary_fall, by = "siteID")

colnames(envir_variables_summary)
View(envir_variables_summary)




#### MERGING METAGENOMIC + SOIL + ENVIR. SUMMARY ####


# Import coordinates
coord_sites <- read.csv("coord_sites.csv", header = T)    
coord_sites <- coord_sites %>% select(-domainID, -Tot_samples, -Final_samples)


##### Table S1 @@@ #####

metadata2017 <- final_table_2017_dedup %>%
  mutate(sample_clean = str_remove(rawDataFileName, "_R[12]\\.fastq\\.gz")) %>%
  left_join(envir_variables_summary, by = "siteID") %>%
  left_join(coord_sites, by = "siteID")

metadata2017 <- metadata2017 %>%
  select(-uid, -namedLocation, -laboratoryName, -sequencingFacilityID, -dnaSampleCode,
         -remarks, -dataQF, -release, -collectDate, -processedDate, -deprecatedVialID, -sequenceAnalysisType,
         -nucleicAcidQuantMethod, -nucleicAcidPurity, -nucleicAcidRange, -dnaPooledStatus, -qaqcStatus, -dnaProcessedBy,
         -eventID, -genomicsSampleCode, -sampleCondition, -sampleTotalReadNumber, -sampleFilteredReadNumber, -TotalReadsAllSamples)

metadata2017 <- metadata2017 %>%
  relocate(siteName, Area, Latitude, Longitude, Type, .after = siteID) %>%
  relocate(soilMoisture_avg, soilTemp_avg, .before = soilTemp_avg_dayMax)

write.csv(metadata2017,"TableS1_metadata2017.csv", row.names = FALSE)

