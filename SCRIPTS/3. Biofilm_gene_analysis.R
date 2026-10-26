### POTENTIAL FOR BIOFILM FORMATION IS WIDESPREAD AND INFLUENCED BY PRECIPITATION VARIABILITY 
### AND ENVIRONMENTAL EXTREMES ACROSS U.S. SOILS 

# Authors: Eduardo M. Choreño-Parra, Nitish Malhotra, Alejandra Rodríguez-Verdugo, 
# Steven D. Allison, and Kathleen K. Treseder.

# ===================================================================================================================== #

# Description: This R script processes raw metagenomic gene abundance data from 112 samples 
# (retained from 119 initially processed) to quantify and summarize the abundance of biofilm-associated genes,
# both overall and by pathway, across samples. 
# It then evaluates the influence of environmental predictors on biofilm-associated gene abundance at 
# the overall, pathway, and individual KO levels.

# Inputs:
# - all_samples.emapper.annotations.tsv.
# - all_genes_long_all_samples.tsv.
# - TableS1_metadata2017.csv.
# - TableS2_summary_pipeline_stats.csv.
# - TableS5_KO_list_pathways.csv.

# Outputs:
# - Summary statistics included in TableS2_summary_pipeline_stats.csv.
# - TableS6_biofilm_normalized_counts.csv 
# - Supplementary Figures S1-S16.
# - Figures 1-4.
# - Table 1.

# Execution:
# Run this script after completing the metagenomic processing pipeline in “2. Metagenomic_analysis.sh” 
# and retrieving the environmental data using “1. NEON_metagenomic_environmental_data.R”. 
# All required input files are compiled in the “DATABASES_AND_METADATA” folder and must be available before execution. 
# Run the sections sequentially, as subsequent analyses depend on objects generated in preceding sections.

# Required packages
# Required packages
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(readr)
library(purrr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(tidytext)
library(broom)
library(glue)
library(ggpubr)
library(lme4)
library(lmerTest)
library(performance)
library(broom.mixed)
library(influence.ME)
library(scales)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(ggnewscale)
library(ggtext)
library(cowplot)
library(grid)




#### 1. IMPORT AND PROCESS DATA ####

# 1) Raw read counts for all CDSs (biofilm-associated and non-biofilm-associated) (TSV version of Table S4)
orf_counts <- read_tsv("all_genes_long_all_samples.tsv", show_col_types = FALSE, progress = FALSE)

# 2) Functional annotations of CDSs including KOs (TSV version of Table S3)
kos_samples <- read_tsv("all_samples.emapper.annotations.tsv", show_col_types = FALSE, progress = FALSE)

# 3) Biofilm KO database (CSV) from Thorpe et al 2025 (CSV version of Table S5)
biofilm_db <- read_csv("TableS5_KO_list_pathways.csv", show_col_types = FALSE)

# 4) Sample metadata and environmental variables (CSV version of Table S1)
metadata2017 <- read.csv("TableS1_metadata2017.csv", header = TRUE)


##### 1.1. Table S4 @@@ #####

# Make raw counts wide and export to create Table S4

orf_counts_wide <- orf_counts %>%
  pivot_wider(
    id_cols = c(orf_id, contig, start, end, strand, length),
    names_from = sample,
    values_from = count,
    values_fill = 0)

write.csv(orf_counts_wide, "TableS4_orf_raw_counts_wide.csv", row.names = FALSE)



##### 1.2. Merge ORF counts with KO annotations #####

# Ensure sample identifiers are stored as character
orf_counts <- orf_counts %>%
  mutate(sample = as.character(sample))

# Clean sample name for functional annotations
kos_samples <- kos_samples %>%
  mutate(sample = str_remove(sample, "_megahit_out$"), KEGG_ko = as.character(KEGG_ko))

# Expand ORFs with multiple KO annotations into separate rows, one per KO
kos_samples_long <- kos_samples %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "", KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = "\\s*,\\s*") %>%
  mutate(
    KO = str_remove(KEGG_ko, "^ko:"),
    KO = str_trim(KO)) %>%
  filter(!is.na(KO), KO != "", KO != "-") %>%
  transmute(
    sample,
    orf_id = query,
    KO,
    ann_preferred = Preferred_name,
    ann_desc = Description) %>%
  distinct()
# Each KO in the functional annotation database can now be linked to its corresponding ORF count and subsequently
# to biofilm metadata, including biofilm pathways.

# Merge KO annotations with raw ORF counts using sample and ORF identifiers
orf_counts_kos <- orf_counts %>%
  left_join(kos_samples_long, by = c("sample", "orf_id"))

# Check the number of distinct ORFs and ORFs with at least one KO annotation
orf_counts_kos %>%
  summarise(
    n_orfs_total = n_distinct(orf_id),
    n_orfs_with_KO = n_distinct(orf_id[!is.na(KO)]))



##### 1.3. Merge ORF counts and annotations with the biofilm database #####

# Match KO annotations with the biofilm database to assign biofilm pathways
# This allows to identify biofilm-associated ORFs while retaining all other ORFs
orf_counts_kos_biofilm <- orf_counts_kos %>%
  left_join(biofilm_db %>% select(KO, pathway), by = "KO") %>%
  mutate(is_biofilm_KO = !is.na(pathway))



##### 1.4. Subset biofilm-associated ORFs #####

# Retain biofilm-associated ORFs, with one row per sample and ORF.
# ORFs with multiple KO annotations are retained in a single row, preserving
# all associated KOs and biofilm pathways while preventing duplication of their read counts.

biofilm_counts_orf <- orf_counts_kos_biofilm %>%
  filter(is_biofilm_KO) %>%
  group_by(sample, orf_id) %>%
  summarise(
    Contigid = dplyr::first(contig),
    start  = dplyr::first(start),
    end    = dplyr::first(end),
    strand = dplyr::first(strand),
    length = dplyr::first(length),
    Count  = dplyr::first(count),
    biofilm_KOs = paste(sort(unique(KO)), collapse = "; "),
    biofilm_pathways = paste(sort(unique(na.omit(pathway))), collapse = "; "),
    ann_preferred = paste(sort(unique(na.omit(ann_preferred))), collapse = "; "),
    ann_desc = paste(sort(unique(na.omit(ann_desc))), collapse = "; "),
    .groups = "drop")

# biofilm_counts_orf: Raw read counts and functional annotations of biofilm-associated
# ORFs, with one row per sample and ORF, preserving multiple KO and pathway assignments.

# Note: This database contains 112 of the original 119 metagenomic samples.
# Seven samples were excluded because their assemblies did not meet the
# minimum quality criteria of 40 contigs and 50,000 bp of assembly length (Table S2)
# Therefore, only the 112 retained samples are used for biofilm gene analyses.



##### 1.5. Calculate summaries for Table S2 @ #####

# Pipeline statistics for TableS2

###### ORF assigned reads ######

orf_assigned_reads_db <- orf_counts %>%
  group_by(sample) %>%
  summarise(
    orf_assigned_reads = sum(count, na.rm = TRUE),
    .groups = "drop")

# Exported these data to be manually added to Table S2
write.csv(orf_assigned_reads_db, "orf_assigned_reads.csv", row.names = F)


###### Total biofilm-related raw reads per sample #######

total_biofilm_reads <- biofilm_counts_orf %>%
  group_by(sample) %>%
  summarise(
    total_biofilm_reads = sum(Count, na.rm = TRUE),
    .groups = "drop")

# Exported these data to be manually added to Table S2
write.csv(total_biofilm_reads, "biofilm_reads.csv", row.names = F)


###### Number of distinct biofilm-associated ORFs per sample ######

n_unique_orfs_per_sample <- biofilm_counts_orf %>%
  group_by(sample) %>%
  summarise(
    n_unique_orfs_biofilm = n_distinct(orf_id)) 

# Exported these data to be manually added to Table S2
write.csv(n_unique_orfs_per_sample, "n_unique_biofilm_orfs_sample.csv", row.names = F)


###### Total number of distinct ORFs per sample ######

total_orfs_per_sample <- orf_counts %>%
  group_by(sample) %>%
  summarise(
    n_total_orfs = n_distinct(orf_id),
    .groups = "drop")

# Exported these data to be manually added to Table S2
write.csv(total_orfs_per_sample, "total_orfs_per_sample.csv", row.names = F)


###### Proportion of biofilm-associated ORFs ######

# Combine and calculate proportion
biofilm_orf_proportion <- total_orfs_per_sample %>%
  left_join(n_unique_orfs_per_sample, by = "sample") %>%
  mutate(
    n_biofilm_orfs = coalesce(n_unique_orfs_biofilm, 0L),
    n_other_orfs = n_total_orfs - n_biofilm_orfs,
    biofilm_orf_percent = n_biofilm_orfs / n_total_orfs * 100)

# Exported these data to be manually added to Table S2
write.csv(biofilm_orf_proportion, "biofilm_orf_proportion.csv", row.names = FALSE)



##### 1.6. Merge sample metadata with pipeline stats  ######

# This step updates metadata2017

# Import completed metagenomic processing pipeline statistics (CSV version of Table S2)
pipeline_summary <- read.csv("TableS2_summary_pipeline_stats.csv", header = TRUE)

# Merge sample metadata and environmental variables with pipeline statistics
metadata2017 <- metadata2017 %>%
  left_join(pipeline_summary %>% select(-genomicsSampleID, -siteID),
    by = c("sample_clean" = "sample"))

colnames(metadata2017)



##### 1.7. Assign sample identifiers and prepare biofilm-associated counts #####

# Create a mapping between raw sample names and genomics sample identifiers
metadata_mapping <- metadata2017 %>%
  mutate(rawDataFileName = str_remove(rawDataFileName, "_R[12]\\.fastq\\.gz")) %>%
  distinct(rawDataFileName, genomicsSampleID,siteID)

# Assign genomics sample identifiers to biofilm-associated counts
biofilm_counts_named <- biofilm_counts_orf %>%
  left_join(
    metadata_mapping,
    by = c("sample" = "rawDataFileName"))

biofilm_counts_named <- biofilm_counts_named %>%
  mutate(sample_site = genomicsSampleID)

# Convert raw counts to wide format, with one column per sample
biofilm_wide <- biofilm_counts_named %>%
  mutate(
    Count  = as.numeric(Count),
    length = as.numeric(length)) %>%
  select(orf_id, Contigid, biofilm_KOs, biofilm_pathways, ann_preferred, ann_desc, sample_site, sample, Count, length) %>%
  pivot_wider(
    id_cols     = c(orf_id, Contigid, length, biofilm_KOs, biofilm_pathways, ann_preferred, ann_desc),
    names_from  = sample_site,
    values_from = Count,
    values_fill = 0)


#### 2. NORMALIZING ####

##### 2.1. BIOFILM-ASSOCIATED ORFs COUNTS #####

# Here, ORF-assigned reads are normalized by gene length and library size


###### 2.1.1. Prepare raw counts and normalization denominators ######

# Store raw counts of biofilm-associated ORFs
counts_raw_full <- biofilm_wide

# Obtain total ORF-assigned reads per sample from the pipeline metadata
mapped_read_counts_named <- metadata2017 %>%
  mutate(rawDataFileName = sample_clean)

# Obtain total ORF-assigned reads per sample for normalization (denominator for normalization by library size)
mapped_read_counts_named$orf_assigned_reads <- as.numeric(mapped_read_counts_named$orf_assigned_reads)

# Create a named vector names = genomicsSampleID, values = orf_assigned_reads
mapped_reads_map <- setNames(
  mapped_read_counts_named$orf_assigned_reads,
  mapped_read_counts_named$genomicsSampleID)

# Extract sample count columns, excluding ORF identifiers and annotations
biofilm_counts_numeric <- counts_raw_full %>%
  select(-orf_id, -Contigid, -length ,-biofilm_KOs, -biofilm_pathways, -ann_preferred, -ann_desc)

# Match normalization denominators to the sample column order
mapped_reads_map <- mapped_reads_map[colnames(biofilm_counts_numeric)]

# Verify that normalization denominators are available for all samples (no NAs in the denominator mapping)
stopifnot(!any(is.na(mapped_reads_map)))
stopifnot(nrow(biofilm_counts_numeric) == length(counts_raw_full$length))



###### 2.1.2. Normalize biofilm-associated ORF raw counts #####

# 1) BY GENE LENGTH (kb)
gene_length_kb <- counts_raw_full$length / 1000
biofilm_counts_numeric <- sweep(biofilm_counts_numeric, 1, gene_length_kb, FUN = "/")

# 2) BY LIBRARY SIZE (total ORF-assigned reads): (gene count / total reads per sample) * 1e6
normalized_counts_numeric <- sweep(biofilm_counts_numeric, 2, mapped_reads_map, FUN = "/") * 1e6

# Reattach ORF identifiers and functional annotations
normalized_counts_df <- cbind(
  biofilm_wide[, c("orf_id", "Contigid", "biofilm_KOs", "biofilm_pathways", 
                   "ann_preferred", "ann_desc")],
  normalized_counts_numeric)

# Remove rows without ORF identifiers
normalized_counts_df <- normalized_counts_df %>%
  filter(!is.na(orf_id))


###### 2.1.3. Verify normalized ORF-assigned counts #####

# annotation columns
anno_cols <- c("orf_id","Contigid","biofilm_KOs",
               "biofilm_pathways","ann_preferred","ann_desc")

sample_cols <- setdiff(names(normalized_counts_df), anno_cols)

# Check for missing normalized abundance values
colSums(is.na(normalized_counts_df[sample_cols]))

# total normalized counts per sample
total_per_sample <- colSums(normalized_counts_df[sample_cols], na.rm = TRUE)
total_per_sample

# range
range(total_per_sample)



##### 2.2. ALL PREDICTED CDSs ####

# Normalize raw read counts of all CDSs by gene length (kb) and
# total ORF-assigned reads per sample, using the same method as for biofilm-associated ORFs
all_cds_normalized <- orf_counts %>%
  filter(sample %in% unique(biofilm_counts_orf$sample)) %>%
  mutate(length_kb = length / 1000) %>%
  left_join(orf_assigned_reads_db, by = "sample") %>%
  mutate(normalized_count = (count / length_kb) / orf_assigned_reads * 1e6)




#### FINAL TABLES ####

##### 1. Raw counts of biofilm-associated ORFs #####

biofilm_wide



##### 2. Normalized abundance of biofilm-associated ORFs (Table S6) @ #####  

normalized_counts_df
# Each ORF is retained as a single row, therefore this table contains ORF-level abundance
# An ORF can contain multiple KO annotations that are treated as a single unit to avoid
# duplication of normalized abundance per ORF

# Export to create Table S6
write.csv(normalized_counts_df, "TableS6_biofilm_normalized_counts.csv", row.names = FALSE)



##### 3. Normalized abundance summarized by KO ##### 

# Sum normalized abundances of ORFs sharing the same KO annotation or
# combination of KO annotations, treating multiple KOs as a single unit.

normalized_counts_summary_ko <- normalized_counts_df %>%
  mutate(
    biofilm_KOs = ifelse(is.na(biofilm_KOs) | biofilm_KOs == "", "unknown", biofilm_KOs),
    biofilm_pathways = ifelse(is.na(biofilm_pathways) | biofilm_pathways == "", "unknown", biofilm_pathways),
    ann_preferred = ifelse(is.na(ann_preferred) | ann_preferred == "", "unknown", ann_preferred),
    ann_desc = ifelse(is.na(ann_desc) | ann_desc == "", "unknown", ann_desc)) %>%
  group_by(biofilm_KOs) %>%
  summarise(
    biofilm_pathways = paste(sort(unique(biofilm_pathways)), collapse = "; "),
    ann_preferred = paste(sort(unique(ann_preferred)), collapse = "; "),
    ann_desc = paste(sort(unique(ann_desc)), collapse = " | "),
    n_orfs = n_distinct(orf_id),
    n_contigs = n_distinct(Contigid),
    across(any_of(sample_cols), ~sum(.x, na.rm = TRUE)),  # use sample_cols from 2.1.3
    .groups = "drop")

# Since multiple KOs within an ORF are treated as a single unit, as in the ORF-level table, 
# summarizing by KO does not duplicate normalized counts. Therefore, both tables should have 
# identical total normalized abundance per sample, despite having different numbers of rows.

# Verify that KO-level aggregation preserves total abundance per sample
all.equal(
  colSums(normalized_counts_df[sample_cols]),
  colSums(normalized_counts_summary_ko[sample_cols]))


##### 4. Normalized abundance of ORFs by pathway #####

# Assign each ORF to all its associated biofilm pathways.
# ORFs associated with multiple pathways contribute their full normalized
# abundance to each pathway, but are counted only once within each pathway.


normalized_counts_pathway_df <- normalized_counts_df %>%
  mutate(
    biofilm_pathways = ifelse(is.na(biofilm_pathways), "", biofilm_pathways),
    biofilm_pathways = str_replace_all(biofilm_pathways, "\\|", ";"),
    biofilm_pathways = str_squish(biofilm_pathways)) %>%
  # Create list column of unique pathways for each ORF
  rowwise() %>%
  mutate(
    biofilm_pathways = list(
      unique(str_trim(unlist(str_split(biofilm_pathways, ";\\s*")))))) %>%
  ungroup() %>%
  # Expand rows (create one row per ORF and associated pathway)
  unnest(biofilm_pathways) %>%
  filter(biofilm_pathways != "")

# Since each ORF can contribute its full normalized abundance to multiple pathways,
# the sum across pathways may exceed the total biofilm-associated gene abundance
# per sample in the original ORF-level table.

# Verify whether pathway-level totals exceed ORF-level totals
colSums(normalized_counts_pathway_df[sample_cols]) >
  colSums(normalized_counts_df[sample_cols])



##### 5. Normalized abundance summarized by pathway ######

# Sum normalized abundance of all ORFs associated with each pathway
pathway_summary <- normalized_counts_pathway_df %>%
  group_by(biofilm_pathways) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

pathway_summary

##### 6. Normalized abundance summarized by pathway and KO #####

# Sum normalized abundances of ORFs sharing the same pathway and KO annotation
# or combination of KOs. Multiple KOs within an ORF are treated as a single unit,
# while ORFs associated with multiple pathways contribute to each pathway.

pathway_ko_summary <- normalized_counts_pathway_df %>%
  mutate(
    biofilm_KOs = ifelse(is.na(biofilm_KOs) | biofilm_KOs == "", "unknown", biofilm_KOs),
    biofilm_pathways = ifelse(is.na(biofilm_pathways) | biofilm_pathways == "", "unknown", biofilm_pathways),
    ann_preferred = ifelse(is.na(ann_preferred) | ann_preferred == "", "unknown", ann_preferred),
    ann_desc = ifelse(is.na(ann_desc) | ann_desc == "", "unknown", ann_desc)) %>%
  group_by(biofilm_pathways, biofilm_KOs) %>%
  summarise(
    ann_preferred = paste(sort(unique(ann_preferred)), collapse = "; "),
    ann_desc      = paste(sort(unique(ann_desc)), collapse = " | "),
    n_orfs        = n_distinct(orf_id),
    n_contigs     = n_distinct(Contigid),
    across(any_of(sample_cols), ~sum(.x, na.rm = TRUE)),
    .groups = "drop")

# Since multiple KOs are treated as a single unit, summarizing by pathway and KO preserves the 
# total normalized abundance per sample from the pathway-level ORF table.

# Verify that total normalized abundance per sample is preserved
all.equal(
  colSums(normalized_counts_pathway_df[sample_cols]),
  colSums(pathway_ko_summary[sample_cols]))



##### Summary of final abundance tables #####

# 1. Raw counts of biofilm-associated ORFs
biofilm_wide
# 2. Normalized abundance of biofilm-associated ORFs
normalized_counts_df
# 3. Normalized abundance summarized by KO
normalized_counts_summary_ko
# 4. Normalized abundance of ORFs by pathway
normalized_counts_pathway_df
# 5. Normalized abundance summarized by pathway
pathway_summary
# 6. Normalized abundance summarized by pathway and KO
pathway_ko_summary





#### ! #####
#### DESCRIPTIVE ANALYSES ####

# This section summarizes the abundance of biofilm-associated genes across metagenomic samples and NEON sites.
# Analyses are conducted at two levels: overall biofilm-associated gene abundance and the relative abundance of 
# individual biofilm-associated pathways.


##### Geographic distribution of NEON sampling sites (Fig1) #####


### Prepare site metadata

# Retain metadata for the 112 metagenomes used in the analyses
metadata2017_final <- metadata2017 %>%
  filter(sample_clean %in% unique(biofilm_counts_orf$sample))

# Verify that all retained samples have matching metadata
stopifnot(all(unique(biofilm_counts_orf$sample) %in% metadata2017$sample_clean))

# Verify the number of retained samples
stopifnot(n_distinct(metadata2017_final$sample_clean) == 112)

# Calculate the number of samples per site
site_table_coord <- metadata2017_final %>%
  group_by(siteID) %>%
  summarise(
    siteName = dplyr::first(siteName),
    Latitude = dplyr::first(Latitude),
    Longitude = dplyr::first(Longitude),
    Type = dplyr::first(Type),
    secPrecip_cv_dayTot = dplyr::first(secPrecip_cv_dayTot),
    Final_samples = n_distinct(sample_clean),
    .groups = "drop")

# Verify the total number of samples and sites
stopifnot(sum(site_table_coord$Final_samples) == 112)
stopifnot(n_distinct(site_table_coord$siteID) == 20)

# Create site labels
site_table_coord <- site_table_coord %>%
  mutate(
    site_label = paste0(siteID, " - ", siteName, " (N = ", Final_samples, ")"),
    site_label_md = paste0("<b>", siteID, "</b>", " - ", siteName, " (N = ", Final_samples, ")"))

### Base map

# Load USA map
usa_map <- ne_states(country = "united states of america", returnclass = "sf")

# Separate Alaska
alaska <- usa_map[usa_map$name == "Alaska", ]
mainland <- usa_map

### Separate sampling sites

# Alaska sites
alaska_sites <- site_table_coord %>%
  filter(Type == "Alaska")

# All sites for the main map
non_alaska_sites <- site_table_coord

# Site colors for the site-name legend
site_colors <- setNames(
  rep("black", length(unique(site_table_coord$site_label))),
  unique(site_table_coord$site_label))

site_colors

### Mainland map

# Extract colored legend
precip_legend_plot <- ggplot(
  non_alaska_sites, aes(x = Longitude, y = Latitude, fill = secPrecip_cv_dayTot)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_viridis_c(name = "Precipitation Day Total CV (%)",
    limits = range( non_alaska_sites$secPrecip_cv_dayTot, na.rm = TRUE),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3, "cm"),
      barheight = unit(0.25, "cm"))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7))

# Convert legend plot to grob
g <- ggplotGrob(precip_legend_plot)

# Find all guide boxes
guide_ids <- grep("guide-box", g$layout$name)

# Extract precipitation legend
precip_legend <- g$grobs[[
  which(g$layout$name == "guide-box-bottom")]]

### Make plot

# Map limits
lon_limits_main <- c(-145, -66)
lat_limits_main <- c(19, 52)

# Create USA outline
usa_outline <- st_union(mainland)

p_mainland <- ggplot() +
  # State boundaries
  geom_sf(data = mainland, fill = "grey95", color = "grey20", linewidth = 0.08) +
  # Country outline
  geom_sf(data = usa_outline, fill = NA, color = "black", linewidth = 0.23) +
  coord_sf(xlim = lon_limits_main, ylim = lat_limits_main, expand = FALSE) +
  # Points colored by precipitation variability
  geom_point(data = non_alaska_sites, aes(x = Longitude, y = Latitude, fill = secPrecip_cv_dayTot),
    alpha = 1, size = 3, shape = 21, color = "black", stroke = 0.6) +
  # Continuous scale for precipitation variability
  scale_fill_viridis_c(name = "Precipitation Day Total CV (%)",
    limits = range(non_alaska_sites$secPrecip_cv_dayTot, na.rm = TRUE), guide = "none", na.value = "white") +
  # Start a new fill scale for the site-name legend
  ggnewscale::new_scale_fill() +
  # Invisible points used to create site-name legend
  geom_point(data = non_alaska_sites, aes(x = Longitude, y = Latitude, fill = site_label),
    alpha = 0, size = 3, shape = 21, color = NA, stroke = 0) +
  # North arrow
  annotation_north_arrow(location = "tl", which_north = "true",
    pad_x = unit(0.35, "cm"),
    pad_y = unit(0.35, "cm"),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering) +
  # Site labels
  geom_label_repel(data = non_alaska_sites %>% filter(!siteID %in% c("HEAL", "DEJU", "TOOL", "BARR")),
    aes(x = Longitude, y = Latitude, label = siteID),
    size = 3.2, fill = "white", fontface = "bold", color = "black",
    seed = 234,
    direction = "both",
    label.padding = unit(0.15, "lines"),
    min.segment.length = 0,
    segment.size = 0.5,
    max.overlaps = 41,
    segment.linetype = 1,
    box.padding = 0.35,
    point.padding = 0.5,
    label.size = 0.3) +
  # Site-name legend
  scale_fill_manual(values = site_colors, name = "", 
                    breaks = site_table_coord$site_label, labels = site_table_coord$site_label_md) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(
    shape = NA, size = 0, alpha = 0, color = NA, fill = NA))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = ggtext::element_markdown(size = 8),
    legend.key = element_blank(),
    legend.key.height = unit(0, "cm"),
    legend.key.width = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.margin = margin(10, 1, 10, 1),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.justification = "center",
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA))

### Add precipitation legend to mainland plot

p_mainland_final <- cowplot::ggdraw() +
  cowplot::draw_plot(p_mainland) +
  cowplot::draw_grob(precip_legend, x = 0.41, y = 0.43, width = 0.25, height = 0.10)

### Alaska inset

lon_limits_ak <- c(-170, -130)
lat_limits_ak <- c(51, 72)

p_alaska <- ggplot() +
  geom_sf(data = alaska, fill = "grey95", color = "grey20", linewidth = 0.23) +
  coord_sf(xlim = lon_limits_ak, ylim = lat_limits_ak, expand = FALSE) +
  # Alaska sites colored by precipitation variability
  geom_point(data = alaska_sites, aes(x = Longitude, y = Latitude, fill = secPrecip_cv_dayTot),
    alpha = 1, size = 3, shape = 21, color = "black", stroke = 0.6) +
  # Use the same precipitation scale as the main map
  scale_fill_viridis_c(limits = range(non_alaska_sites$secPrecip_cv_dayTot, na.rm = TRUE),
    guide = "none", na.value = "white") +
  # Alaska site labels
  geom_label_repel(data = alaska_sites,aes(x = Longitude, y = Latitude, label = siteID),
    seed = 456,
    size = 3.0,
    fill = "white",
    fontface = "bold",
    color = "black",
    direction = "both",
    label.padding = unit(0.1, "lines"),
    min.segment.length = 0,
    segment.size = 0.5,
    max.overlaps = 20,
    box.padding = 0.25,
    point.padding = 0.3,
    label.size = 0.3) +
  scale_x_continuous(breaks = c(-165, -150, -135)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA))

### Combine mainland map and Alaska inset

final_map_w_color <- ggdraw() +
  draw_plot(p_mainland_final, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(p_alaska, x = 0.08, y = 0.43, width = 0.23, height = 0.30)

##### Fig. 1 @ #####

ggsave("Figure_1.pdf", plot = final_map_w_color, device = cairo_pdf, units = "in", width = 7, height = 6, dpi = 1000)
ggsave("Figure_1.png", plot = final_map_w_color, units = "in", width = 7, height = 6, dpi = 1000)




#### ! #####
#### 3. OVERALL-LEVEL ####


##### 3.1. Prepare normalized biofilm abundance data #####

# Convert normalized KO-level abundance to long format

anno_cols <- c("biofilm_KOs","biofilm_pathways",
               "ann_preferred","ann_desc",
               "n_orfs","n_contigs")

sample_cols <- setdiff(names(normalized_counts_summary_ko), anno_cols)

biofilm_long <- normalized_counts_summary_ko %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "genomicsSampleID",
    values_to = "count")

# Add site information to each sample
metadata_dedup <- metadata2017 %>%
  distinct(genomicsSampleID, siteID)

biofilm_long <- biofilm_long %>%
  left_join(metadata_dedup, by = "genomicsSampleID")



##### 3.2. Total abundance per sample ##### 

# Summarize normalized abundance across all KO groups within each sample (total normalized abundance per sample)
biofilm_signal_per_sample <- biofilm_long %>%
  group_by(genomicsSampleID, siteID) %>%
  summarise(
    total_count_sample = sum(count, na.rm = TRUE),   # total biofilm signal per sample
    .groups = "drop")

sort(biofilm_signal_per_sample$total_count_sample)


###### Fig.S1A @ ######

# Create shortened sample labels
biofilm_signal_per_sample_plot <- biofilm_signal_per_sample %>%
  mutate(sample_short = sub("-O-.*", "", genomicsSampleID))

# Plot
p_total_sample <- ggplot(biofilm_signal_per_sample_plot,
                         aes(x = sample_short, y = total_count_sample)) +
  geom_col(width = 0.95, fill = "black", color = "white", linewidth = 0.2) +
  facet_grid(. ~ siteID, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, 250000),
    breaks = seq(0, 250000, by = 50000)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        axis.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(
          size = 7, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 9, face = "bold", angle = 90),
        strip.background = element_blank(),
        axis.line = element_line(linewidth = 0.3)) +
  labs(x = "Sample", y = "Total Gene Abundance\n(counts / kb / million reads)")

plot(p_total_sample)



##### 3.3. Relative abundance of biofilm-associated genes @ #####

# Total normalized abundance of all predicted CDSs per sample
total_all_cds_normalized <- all_cds_normalized %>%  # use from section 2.2
  group_by(sample) %>%
  summarise(total_normalized_all_cds = sum(normalized_count, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(metadata_mapping, by = c("sample" = "rawDataFileName"))

# Calculate the proportion of normalized abundance associated with biofilm genes
biofilm_proportion <- biofilm_signal_per_sample %>%
  left_join(
    total_all_cds_normalized %>% select(-siteID),
    by = "genomicsSampleID") %>%
  mutate(
    biofilm_percent =
      total_count_sample / total_normalized_all_cds * 100,
    nonbiofilm_percent =
      100 - biofilm_percent)


###### Fig.S1B @ ######

# Prepare relative abundance data for plotting
biofilm_proportion_long <- biofilm_proportion %>%
  select(
    genomicsSampleID,
    siteID,
    biofilm_percent,
    nonbiofilm_percent) %>%
  pivot_longer(
    cols = c(biofilm_percent, nonbiofilm_percent),
    names_to = "gene_category",
    values_to = "percent") %>%
  mutate(
    gene_category = recode(
      gene_category,
      biofilm_percent    = "Biofilm-associated",
      nonbiofilm_percent = "Other"),
    gene_category = factor(
      gene_category,
      levels = c("Other", "Biofilm-associated")),
    sample_short = sub("-O-.*", "", genomicsSampleID))

# Define colors
biofilm_gene_colors <- c(
  "Biofilm-associated" = "black",
  "Other"              = "gray90")

# Plot
p_relative_biofilm_sample <- ggplot(
  biofilm_proportion_long,
  aes(x = sample_short, y = percent / 100, fill = gene_category)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.2) +
  facet_grid(. ~ siteID, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = biofilm_gene_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold", angle = 90),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right") +
  labs(x = "Sample", y = "Relative Gene Abundance", fill = "Gene Category")



#### ! #####
#### 4. PATHWAY-LEVEL ####


##### 4.1 Rel. abundance of biofilm-associated pathways ####

# Convert pathway-level normalized abundance per sample to long format
pathway_long <- pathway_summary %>%
  pivot_longer(
    cols = -biofilm_pathways,
    names_to = "genomicsSampleID",
    values_to = "total_count_sample")

# Add site information
pathway_long <- pathway_long %>%
  left_join(metadata_dedup, by = "genomicsSampleID")


# Define clean pathway names and set order
pathway_labels <- c(
  "abc_transporters"     = "ABC Transporters",
  "BBSdb"                = "BBSdb",
  "biofilm_formation"    = "Biofilm Formation",
  "chemotaxis"           = "Chemotaxis",
  "flagellar_assembly"   = "Flagellar Assembly",
  "polysac_eps"          = "EPS Biosynthesis",
  "quorum_sensing"       = "Quorum Sensing",
  "two_component_system" = "Two-Component System")

pathway_long <- pathway_long %>%
  mutate(
    biofilm_pathways_clean = dplyr::recode(
      biofilm_pathways,
      !!!pathway_labels,
      .default = biofilm_pathways),
    biofilm_pathways_clean = factor(
      biofilm_pathways_clean,
      levels = c(
        "ABC Transporters",
        "Two-Component System",
        "Biofilm Formation",
        "EPS Biosynthesis",
        "Quorum Sensing",
        "Chemotaxis",
        "Flagellar Assembly",
        "BBSdb")))


###### FigS1C @ ######

# Calculate relative pathway abundance within each sample
pathway_sample_relative <- pathway_long %>%
  group_by(genomicsSampleID) %>%
  mutate(
    sample_total = sum(total_count_sample, na.rm = TRUE),
    rel_abund = ifelse(
      sample_total > 0,
      total_count_sample / sample_total,
      NA_real_)) %>%
  ungroup() %>%
  mutate(
    sample_short = sub("-O-.*", "", genomicsSampleID))

# Define pathway colors
pathway_colors <- c(
  "ABC Transporters"     = "#6889C0",
  "Two-Component System" = "#A0CBE8",
  "Biofilm Formation"    = "#FDF2C7",
  "EPS Biosynthesis"     = "#FBCE1D",
  "Quorum Sensing"       = "#59A14F",
  "Chemotaxis"           = "#8CD17D",
  "Flagellar Assembly"   = "#B54745",
  "BBSdb"                = "#D68A8A")

# Plot
p_relative_pathway_sample <- ggplot(pathway_sample_relative, aes(x = sample_short,
                                                                 y = rel_abund, fill = biofilm_pathways_clean)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pathway_colors) +
  facet_grid(. ~ siteID, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold", angle = 90),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right") +
  labs(x = "Sample", y = "Relative Pathway Abundance", fill = "Pathway")



#### 3 + 4 = FigS1 @ ####

sample_plots_all3 <- ggarrange(
  p_total_sample + labs(x = NULL , title = "Site")+
    theme(plot.title = element_text(
      size = 10, face = "bold", hjust = 0.5)),
  p_relative_biofilm_sample + labs(x = NULL) + theme(strip.text.x = element_blank()),
  p_relative_pathway_sample  + theme(strip.text.x = element_blank()),
  ncol = 1,
  nrow = 3,
  labels = LETTERS[1:3],
  align = "v")

plot(sample_plots_all3)

ggsave("Figure_S1.png", plot = sample_plots_all3,
       width = 15, height = 10, units = "in", dpi = 1000)

# Summary of abundance and pathways
summary(biofilm_proportion$biofilm_percent)
summary_pathways <- pathway_sample_relative %>%
  group_by(biofilm_pathways_clean) %>%
  summarise(
    mean_rel_abund = mean(rel_abund, na.rm = TRUE),
    median_rel_abund = median(rel_abund, na.rm = TRUE),
    min_rel_abund = min(rel_abund, na.rm = TRUE),
    max_rel_abund = max(rel_abund, na.rm = TRUE),
    n = sum(!is.na(rel_abund)),
    sd = sd(rel_abund, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop") %>%
  arrange(desc(mean_rel_abund))




#### ! #####
#### ABUNDANCE VS ENVIRONMENTAL PARAMETERS ####

# This section evaluates relationships between environmental variables and the normalized abundance of
# biofilm-associated genes at the overall, pathway, and individual KO levels.


#### ! #####
#### 5. OVERALL - LEVEL ####

# This section evaluates the effects of environmental variables on the total abundance 
# of biofilm-associated genes per sample.

# Sections 5.1-5.5. fit and compare linear and quadratic mixed-effects models for each environmental variable.
# For each variable, the model with the lowest AIC is selected and plotted. Then, the predictor with the lowest
# AIC within each environmental category (soil temperature, air temperature, soil mositure, precipitation) 
# is identified and plotted.

# Section 5.6. evaluates the combined effects of selected environmental predictors using multiple mixed-effects models.
# Additive and two-way interaction models are fitted and compared using AIC.


##### 5.1. Prepare data for modeling #####

# Merge biofilm-associated gene abundance with environmental variables
biofilm_env <- biofilm_long %>%
  left_join(
    metadata2017 %>%
      distinct(genomicsSampleID, .keep_all = TRUE) %>%
      select(-siteID),
    by = "genomicsSampleID")

# Verify that the join did not duplicate observations
stopifnot(nrow(biofilm_env) == nrow(biofilm_long))

# Verify the number of metagenomic samples
n_distinct(biofilm_env$genomicsSampleID)

# Calculate total biofilm-associated gene abundance per sample
# and summarize environmental variables for subsequent modeling

biofilm_env_samples <- biofilm_env %>%
  group_by(genomicsSampleID) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    siteID = dplyr::first(siteID),
    
    # Average environmental variables within each sample (they should be constant)
    across(
      c(
        soilTemp_avg,
        starts_with("soilTemp_avg_"),
        starts_with("soilTemp_cv_"),
        soilMoisture_avg,
        starts_with("soilMoist_avg_"),
        starts_with("soilMoist_cv_"),
        starts_with("airTemp_avg_"),
        starts_with("airTemp_cv_"),
        starts_with("secPrecip_avg_"),
        secPrecip_cv_dayTot),
      ~ mean(.x, na.rm = TRUE)),
    .groups = "drop") %>%
  mutate(
    # coefficients of variation as nonnegative magnitudes
    across(
      c(
        starts_with("soilTemp_cv_"),
        starts_with("soilMoist_cv_"),
        starts_with("airTemp_cv_"),
        secPrecip_cv_dayTot),
      abs))

# Inspect the resulting dataset
summary(biofilm_env_samples)

# Verify the number of samples
nrow(biofilm_env_samples)

# Center environmental parameters for quadratic models ###
biofilm_env_samples <- biofilm_env_samples %>%
  mutate(
    across(
      .cols = -c(genomicsSampleID, siteID, total_count),
      .fns = ~ as.numeric(scale(.x, center = TRUE, scale = FALSE)),
      .names = "{.col}_c"))

# List environmental predictors, excluding centered copies
predictors <- biofilm_env_samples %>%
  select(-genomicsSampleID, -siteID, -total_count) %>%
  select(-ends_with("_c")) %>%
  names()

length(predictors)
predictors

# Define descriptive x-axis labels for environmental predictors
x_labels <- list(
  
  # Soil Temperature
  soilTemp_avg = "Soil Temperature (°C)",
  soilTemp_avg_dayMax = "Soil Temperature Day Max (°C)",
  soilTemp_avg_dayMean = "Soil Temperature Day Mean (°C)",
  soilTemp_avg_dayMin = "Soil Temperature Day Min (°C)",
  soilTemp_cv_dayMax = "Soil Temperature Day Max CV (%)",
  soilTemp_cv_dayMean = "Soil Temperature Day Mean CV (%)",
  soilTemp_cv_dayMin = "Soil Temperature Day Min CV (%)",
  soilTemp_avg_dayMax_winter = "Soil Temperature Day Max (Winter, °C)",
  soilTemp_avg_dayMean_winter = "Soil Temperature Day Mean (Winter, °C)",
  soilTemp_avg_dayMin_winter = "Soil Temperature Day Min (Winter, °C)",
  soilTemp_avg_dayMax_spring = "Soil Temperature Day Max (Spring, °C)",
  soilTemp_avg_dayMean_spring = "Soil Temperature Day Mean (Spring, °C)",
  soilTemp_avg_dayMin_spring = "Soil Temperature Day Min (Spring, °C)",
  soilTemp_avg_dayMax_summer = "Soil Temperature Day Max (Summer, °C)",
  soilTemp_avg_dayMean_summer = "Soil Temperature Day Mean (Summer, °C)",
  soilTemp_avg_dayMin_summer = "Soil Temperature Day Min (Summer, °C)",
  soilTemp_avg_dayMax_fall = "Soil Temperature Day Max (Fall, °C)",
  soilTemp_avg_dayMean_fall = "Soil Temperature Day Mean (Fall, °C)",
  soilTemp_avg_dayMin_fall = "Soil Temperature Day Min (Fall, °C)",
  
  # Soil Moisture
  soilMoisture_avg = "Soil Moisture (%)",
  soilMoist_avg_dayMax = "Soil Moisture Day Max (%)",
  soilMoist_avg_dayMean = "Soil Moisture Day Mean (%)",
  soilMoist_avg_dayMin = "Soil Moisture Day Min (%)",
  soilMoist_cv_dayMax = "Soil Moisture Day Max CV (%)",
  soilMoist_cv_dayMean = "Soil Moisture Day Mean CV (%)",
  soilMoist_cv_dayMin = "Soil Moisture Day Min CV (%)",
  soilMoist_avg_dayMax_winter = "Soil Moisture Day Max (Winter, %)",
  soilMoist_avg_dayMean_winter = "Soil Moisture Day Mean (Winter, %)",
  soilMoist_avg_dayMin_winter = "Soil Moisture Day Min (Winter, %)",
  soilMoist_avg_dayMax_spring = "Soil Moisture Day Max (Spring, %)",
  soilMoist_avg_dayMean_spring = "Soil Moisture Day Mean (Spring, %)",
  soilMoist_avg_dayMin_spring = "Soil Moisture Day Min (Spring, %)",
  soilMoist_avg_dayMax_summer = "Soil Moisture Day Max (Summer, %)",
  soilMoist_avg_dayMean_summer = "Soil Moisture Day Mean (Summer, %)",
  soilMoist_avg_dayMin_summer = "Soil Moisture Day Min (Summer, %)",
  soilMoist_avg_dayMax_fall = "Soil Moisture Day Max (Fall, %)",
  soilMoist_avg_dayMean_fall = "Soil Moisture Day Mean (Fall, %)",
  soilMoist_avg_dayMin_fall = "Soil Moisture Day Min (Fall, %)",
  
  # Air Temperature
  airTemp_avg_dayMax = "Air Temperature Day Max (°C)",
  airTemp_avg_dayMean = "Air Temperature Day Mean (°C)",
  airTemp_avg_dayMin = "Air Temperature Day Min (°C)",
  airTemp_cv_dayMax = "Air Temperature Day Max CV (%)",
  airTemp_cv_dayMean = "Air Temperature Day Mean CV (%)",
  airTemp_cv_dayMin = "Air Temperature Day Min CV (%)",
  airTemp_avg_dayMax_winter = "Air Temperature Day Max (Winter, °C)",
  airTemp_avg_dayMean_winter = "Air Temperature Day Mean (Winter, °C)",
  airTemp_avg_dayMin_winter = "Air Temperature Day Min (Winter, °C)",
  airTemp_avg_dayMax_spring = "Air Temperature Day Max (Spring, °C)",
  airTemp_avg_dayMean_spring = "Air Temperature Day Mean (Spring, °C)",
  airTemp_avg_dayMin_spring = "Air Temperature Day Min (Spring, °C)",
  airTemp_avg_dayMax_summer = "Air Temperature Day Max (Summer, °C)",
  airTemp_avg_dayMean_summer = "Air Temperature Day Mean (Summer, °C)",
  airTemp_avg_dayMin_summer = "Air Temperature Day Min (Summer, °C)",
  airTemp_avg_dayMax_fall = "Air Temperature Day Max (Fall, °C)",
  airTemp_avg_dayMean_fall = "Air Temperature Day Mean (Fall, °C)",
  airTemp_avg_dayMin_fall = "Air Temperature Day Min (Fall, °C)",
  
  # Precipitation
  secPrecip_avg_dayTot = "Precipitation Day Total (mm)",
  secPrecip_cv_dayTot = "Precipitation Day Total CV (%)",
  secPrecip_avg_dayTot_winter = "Precipitation Day Total (Winter, mm)",
  secPrecip_avg_dayTot_spring = "Precipitation Day Total (Spring, mm)",
  secPrecip_avg_dayTot_summer = "Precipitation Day Total (Summer, mm)",
  secPrecip_avg_dayTot_fall = "Precipitation Day Total (Fall, mm)")

# Verify that all predictors have corresponding axis labels
stopifnot(setequal(predictors, names(x_labels)))


##### 5.2. Screening linear and quadratic mixed-effects models #####

###### 5.2.1. Helper functions ###### 

# Fit linear mixed-effects models and handle fitting errors
safe_lmer <- function(formula, data) {
  tryCatch(lmer(formula, data = data, REML = FALSE), error = function(e) NULL)}

# Extract fixed-effect statistics
safe_tidy_fixed <- function(model) {
  if (is.null(model)) return(tibble(term = character(), estimate = numeric(), std.error = numeric(), p.value = numeric()))
  tryCatch(broom.mixed::tidy(model, effects = "fixed"),
           error = function(e) tibble(term = character(), estimate = numeric(), std.error = numeric(), p.value = numeric()))}

# Extract marginal R-squared
safe_r2_marginal <- function(model) {
  if (is.null(model)) return(NA_real_)
  out <- tryCatch(performance::r2(model)$R2_marginal, error = function(e) NA_real_)
  as.numeric(out)}

# Extract AIC
safe_AIC <- function(model) {
  if (is.null(model)) return(NA_real_)
  tryCatch(AIC(model), error = function(e) NA_real_)}

# Check for singular model fits
safe_singular <- function(model) {
  if (is.null(model)) return(NA)
  tryCatch(isSingular(model), error = function(e) NA)}

# Extract coefficient, standard error, and p-value for a model term
get_term <- function(tt, term) {
  if (!term %in% tt$term) return(list(est = NA_real_, se = NA_real_, p = NA_real_))
  list(
    est = tt$estimate[tt$term == term][1],
    se  = tt$std.error[tt$term == term][1],
    p   = tt$p.value[tt$term == term][1])}


###### 5.2.2. Function for screening ######

# Screen one environmental predictor using linear and quadratic mixed-effects models,
# with siteID as a random intercept to account for non-independence of samples within sites

# screen one predictor (returns ONE row with BOTH models)
screen_one <- function(dat, response, v) {
  # Identify the centered predictor created in Section 5.1
  v_c <- paste0(v, "_c") 
  # Center
  # dat2 <- dat %>%
  #   mutate(!!v_c := as.numeric(scale(.data[[v]], center = TRUE, scale = FALSE)))
  # Define linear and quadratic models
  f_lin  <- as.formula(paste0(response, " ~ ", v,   " + (1 | siteID)"))
  f_quad <- as.formula(paste0(response, " ~ ", v_c, " + I(", v_c, "^2) + (1 | siteID)"))
  # Fit models
  m_lin  <- safe_lmer(f_lin,  dat)
  m_quad <- safe_lmer(f_quad, dat)
  # Extract fixed-effect statistics
  tt_lin  <- safe_tidy_fixed(m_lin)
  tt_quad <- safe_tidy_fixed(m_quad)
  # Linear model term
  lin_term_lin <- get_term(tt_lin, v)
  # Quadratic model terms
  lin_term_quad <- get_term(tt_quad, v_c)
  quad_term     <- get_term(tt_quad, paste0("I(", v_c, "^2)"))
  # Summarize results from both models
  tibble(
    predictor = v,
    response = response,
    # AIC
    AIC_linear = safe_AIC(m_lin),
    AIC_quadratic = safe_AIC(m_quad),
    deltaAIC_quad_minus_lin = AIC_quadratic - AIC_linear,
    # R2 marginal
    R2m_linear = safe_r2_marginal(m_lin),
    R2m_quadratic = safe_r2_marginal(m_quad),
    # Singularity flags
    singular_linear = safe_singular(m_lin),
    singular_quadratic = safe_singular(m_quad),
    
    # Linear model coefficients
    coef_linear_model = round(lin_term_lin$est, 4),
    se_linear_model   = round(lin_term_lin$se,  4),
    p_linear_model_num = lin_term_lin$p,
    p_linear_model    = ifelse(is.na(lin_term_lin$p), NA_character_, scales::pvalue(lin_term_lin$p)),
    
    # Quadratic model linear term (centered)
    coef_linear_quad_model = round(lin_term_quad$est, 4),
    se_linear_quad_model   = round(lin_term_quad$se,  4),
    p_linear_quad_model    = ifelse(is.na(lin_term_quad$p), NA_character_, scales::pvalue(lin_term_quad$p)),
    
    # Quadratic term
    coef_quadratic = round(quad_term$est, 6),
    se_quadratic   = round(quad_term$se,  6),
    p_quadratic_num = quad_term$p,
    p_quadratic    = ifelse(is.na(quad_term$p), NA_character_, scales::pvalue(quad_term$p)),
    
    # Model fitting notes
    note = case_when(
      is.na(AIC_linear) & is.na(AIC_quadratic) ~ "both models failed",
      is.na(AIC_quadratic) ~ "quadratic model failed",
      is.na(AIC_linear) ~ "linear model failed",
      TRUE ~ ""))}


###### 5.2.3. Screen environmental predictors ######

# Define response variable
response <- "total_count"   

# Fit linear and quadratic mixed-effects models for all 62 predictors
screening_models_total <- bind_rows(lapply(predictors, function(v) {
  screen_one(biofilm_env_samples, response = response, v = v)}))

# Classify predictors by environmental category
screening_models_total <- screening_models_total %>%
  mutate(
    variable_type = case_when(
      grepl("^soilTemp_", predictor) ~ "Soil Temperature",
      grepl("^airTemp_", predictor) ~ "Air Temperature",
      grepl("^soilMoist", predictor) ~ "Soil Moisture",
      grepl("^secPrecip_", predictor) ~ "Secondary Precipitation",
      TRUE ~ "Other"))

View(screening_models_total)



##### 5.3. Model selection based on AIC ######

# Select the linear or quadratic model with the lowest AIC for each predictor
screening_models_total <- screening_models_total %>%
  mutate(
    best_model = case_when(
      is.na(AIC_linear) & is.na(AIC_quadratic) ~ NA_character_,
      is.na(AIC_quadratic) ~ "Linear",
      is.na(AIC_linear) ~ "Polynomial",
      AIC_quadratic < AIC_linear ~ "Polynomial",
      TRUE ~ "Linear"),
    # AIC of the selected model
    AIC_selected = case_when(
      best_model == "Polynomial" ~ AIC_quadratic,
      best_model == "Linear" ~ AIC_linear,
      TRUE ~ NA_real_),
    # Calculate delta AIC relative to the best model across all predictors
    deltaAIC_selected_vs_best = 
      AIC_selected - min(AIC_selected, na.rm = TRUE)) 

View(screening_models_total)


###### 5.3.1. Plots for all predictors ######

# Function: Generate plots using the AIC-selected linear or quadratic model

format_p <- function(p) {
  p <- trimws(as.character(p))
  ifelse(startsWith(p, "<"), paste0("P ", p), paste0("P = ", p))}

make_selected_plot <- function(v, show_AIC = TRUE) {
  # Get model selection information
  model_info <- screening_models_total %>%
    filter(predictor == v)
  selected_model <- model_info$best_model
  
  # Choose color according to environmental category
  smooth_color <- if (
    model_info$variable_type %in% c(
      "Soil Temperature",
      "Air Temperature")) {
    "#f1a226"} else {"#298c8c"}
  
  # Extract statistics from the selected model
  if (selected_model == "Linear") {
    # Formula
    smooth_formula <- y ~ x
    # Is selected relationship significant?
    significant_relationship <- 
      !is.na(model_info$p_linear_model_num) &&
      model_info$p_linear_model_num < 0.05
    # Model information
    coef_lin <- sprintf("%.2f", model_info$coef_linear_model)
    p_lin    <- model_info$p_linear_model
    r2_value  <- model_info$R2m_linear
    aic_value <- model_info$AIC_selected
    # Format R2
    r2_label <- ifelse(
      r2_value < 0.001,
      "<0.001",
      sprintf("%.3f", r2_value))
    # Format AIC
    aic_label <- sprintf("%.2f", aic_value)
    # Annotation
    model_label <- glue(
      "Linear: {coef_lin} ({format_p(p_lin)})\n",
      "R²: {r2_label}")
    
    if (show_AIC) {
      model_label <- paste0(model_label, "\nAIC: ", aic_label)}
    
    } else if (selected_model == "Polynomial") {
    # Formula
    smooth_formula <- y ~ x + I(x^2)
    # Is quadratic relationship significant?
    significant_relationship <- 
      !is.na(model_info$p_quadratic_num) &&
      model_info$p_quadratic_num < 0.05
    # Model information
    coef_lin <- sprintf(
      "%.2f",
      model_info$coef_linear_quad_model)
    p_lin <- model_info$p_linear_quad_model
    coef_quad <- sprintf(
      "%.2f",
      model_info$coef_quadratic)
    p_quad <- model_info$p_quadratic
    r2_value  <- model_info$R2m_quadratic
    aic_value <- model_info$AIC_selected
    # Format R2
    r2_label <- ifelse(
      r2_value < 0.001,
      "<0.001",
      sprintf("%.3f", r2_value))
    # Format AIC
    aic_label <- sprintf("%.2f", aic_value)
    # Annotation
    model_label <- glue(
      "Linear: {coef_lin} ({format_p(p_lin)})\n",
      "Quadratic: {coef_quad} ({format_p(p_quad)})\n",
      "R²: {r2_label}")
    
    if (show_AIC) {
      model_label <- paste0(model_label, "\nAIC: ", aic_label)}
    }
  
  # Get environmental variable labels
  xlab <- if (!is.null(x_labels[[v]])) {
    x_labels[[v]]
  } else {
    v
  }
  
  # Remove "Secondary" from precipitation axis labels
  if (grepl("^secPrecip_", v)) {
    xlab <- gsub("Secondary ", "", xlab)}
  
  # Base plot
  p <- ggplot(
    biofilm_env_samples,
    aes(x = .data[[v]], y = total_count))
  
  # Add fitted curve only for significant relationships
  if (significant_relationship) {
    p <- p +
      geom_smooth(
        method = "lm",
        formula = smooth_formula,
        se = TRUE,
        color = "black",
        fill = smooth_color,
        linewidth = 0.7,
        alpha = 0.4)}
  
  # Add points, model statistics, and formatting
  p <- p + geom_point(size = 1.5, alpha = 1, shape = 21, color = "grey10", fill = "gray50", stroke = 0.35) +
    annotate("label", x = Inf, y = -Inf, label = model_label,
      hjust = 1, vjust = -0.2, size = 2.3, linewidth = 0,
      fill = scales::alpha("white", 0),
      color = "black") +
    labs(x = xlab, y = "Total Gene Abundance\n(counts / kb / million reads)") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold"),
      legend.margin = margin(0, 0, 0, 0))
  return(p)}


# Create individual plots

for (v in predictors) {plot_name <- paste0("total_counts_", v, "_plot")
  assign(plot_name, make_selected_plot(v), envir = .GlobalEnv)}



###### Soil Temperature (FigS2) @ #####

soilTemp_all_total <- ggarrange(
  
  # Overall (non-seasonal)
  total_counts_soilTemp_avg_plot,
  total_counts_soilTemp_avg_dayMean_plot,
  total_counts_soilTemp_avg_dayMax_plot,
  total_counts_soilTemp_avg_dayMin_plot,
  
  # DayMean by season
  total_counts_soilTemp_avg_dayMean_winter_plot,
  total_counts_soilTemp_avg_dayMean_spring_plot,
  total_counts_soilTemp_avg_dayMean_summer_plot,
  total_counts_soilTemp_avg_dayMean_fall_plot,
  
  # DayMax by season
  total_counts_soilTemp_avg_dayMax_winter_plot,
  total_counts_soilTemp_avg_dayMax_spring_plot,
  total_counts_soilTemp_avg_dayMax_summer_plot,
  total_counts_soilTemp_avg_dayMax_fall_plot,
  
  # DayMin by season
  total_counts_soilTemp_avg_dayMin_winter_plot,
  total_counts_soilTemp_avg_dayMin_spring_plot,
  total_counts_soilTemp_avg_dayMin_summer_plot,
  total_counts_soilTemp_avg_dayMin_fall_plot,
  
  # CV (non-seasonal)
  total_counts_soilTemp_cv_dayMax_plot,
  total_counts_soilTemp_cv_dayMean_plot,
  total_counts_soilTemp_cv_dayMin_plot,
  
  ncol = 4,
  nrow = 5,
  labels = LETTERS[1:19])

plot(soilTemp_all_total)

ggsave("Figure_S2.png", plot = soilTemp_all_total,
       width = 14, height = 12, units = "in",dpi = 1000)


###### Air Temperature (FigS3) @  #####

airTemp_all_total <- ggarrange(
  # Overall (non-seasonal)
  total_counts_airTemp_avg_dayMean_plot,
  total_counts_airTemp_avg_dayMax_plot,
  total_counts_airTemp_avg_dayMin_plot,
  # Daytotal by season
  total_counts_airTemp_avg_dayMean_winter_plot,
  total_counts_airTemp_avg_dayMean_spring_plot,
  total_counts_airTemp_avg_dayMean_summer_plot,
  total_counts_airTemp_avg_dayMean_fall_plot,
  # DayMax by season
  total_counts_airTemp_avg_dayMax_winter_plot,
  total_counts_airTemp_avg_dayMax_spring_plot,
  total_counts_airTemp_avg_dayMax_summer_plot,
  total_counts_airTemp_avg_dayMax_fall_plot,
  # DayMin by season
  total_counts_airTemp_avg_dayMin_winter_plot,
  total_counts_airTemp_avg_dayMin_spring_plot,
  total_counts_airTemp_avg_dayMin_summer_plot,
  total_counts_airTemp_avg_dayMin_fall_plot,
  # CV (non-seasonal)
  total_counts_airTemp_cv_dayMax_plot,
  total_counts_airTemp_cv_dayMean_plot,
  total_counts_airTemp_cv_dayMin_plot,
  ncol = 4, nrow = 5,
  labels = LETTERS[1:18])

airTemp_all_total

ggsave("Figure_S3.png", plot = airTemp_all_total,
       width = 14, height = 12, units = "in", dpi = 1000)


###### Soil Moisture (FigS4) @  #####

soilMoist_all_total <- ggarrange(
  # Overall (non-seasonal)
  total_counts_soilMoisture_avg_plot,
  total_counts_soilMoist_avg_dayMean_plot,
  total_counts_soilMoist_avg_dayMax_plot,
  total_counts_soilMoist_avg_dayMin_plot,
  # DayMean by season
  total_counts_soilMoist_avg_dayMean_winter_plot,
  total_counts_soilMoist_avg_dayMean_spring_plot,
  total_counts_soilMoist_avg_dayMean_summer_plot,
  total_counts_soilMoist_avg_dayMean_fall_plot,
  # DayMax by season
  total_counts_soilMoist_avg_dayMax_winter_plot,
  total_counts_soilMoist_avg_dayMax_spring_plot,
  total_counts_soilMoist_avg_dayMax_summer_plot,
  total_counts_soilMoist_avg_dayMax_fall_plot,
  # DayMin by season
  total_counts_soilMoist_avg_dayMin_winter_plot,
  total_counts_soilMoist_avg_dayMin_spring_plot,
  total_counts_soilMoist_avg_dayMin_summer_plot,
  total_counts_soilMoist_avg_dayMin_fall_plot,
  # CV (non-seasonal)
  total_counts_soilMoist_cv_dayMax_plot,
  total_counts_soilMoist_cv_dayMean_plot,
  total_counts_soilMoist_cv_dayMin_plot,
  ncol = 4, nrow = 5,
  labels = LETTERS[1:19])

plot(soilMoist_all_total)

ggsave("Figure_S4.png", plot = soilMoist_all_total,
       width = 14, height = 12, units = "in", dpi = 1000)


###### Precipitation (FigS5) @  ######

precip_all_total <- ggarrange(
  # Secondary precip overall + CV
  total_counts_secPrecip_avg_dayTot_plot,
  total_counts_secPrecip_cv_dayTot_plot,
  # Secondary precip by season
  total_counts_secPrecip_avg_dayTot_winter_plot,
  total_counts_secPrecip_avg_dayTot_spring_plot,
  total_counts_secPrecip_avg_dayTot_summer_plot,
  total_counts_secPrecip_avg_dayTot_fall_plot,
  ncol = 3, nrow = 2,
  labels = LETTERS[1:6])

precip_all_total

ggsave("Figure_S5.png", plot = precip_all_total,
       width = 9, height = 4.8, units = "in", dpi = 1000)



##### 5.4. Best predictor within environmental category #####

# Select the predictor with the lowest AIC within each environmental category 
# (from all the selected variables based on AIC)
best_per_type_total <- screening_models_total %>%
  group_by(variable_type) %>%
  slice_min(AIC_selected, n = 1, with_ties = FALSE) %>%
  arrange(AIC_selected)

view(best_per_type_total)
best_per_type_total %>% select(predictor, best_model)








##### 5.5. Plot the best environmental predictors #####

# Extract the best predictor from each environmental category
best_predictors_overall <- best_per_type_total$predictor
best_predictors_overall

# Generate individual plots without AIC annotations
best_soilMoist_plot <- make_selected_plot(
  "soilMoist_avg_dayMin", show_AIC = FALSE)
best_precip_plot <- make_selected_plot(
  "secPrecip_cv_dayTot", show_AIC = FALSE)
best_soilTemp_plot <- make_selected_plot(
  "soilTemp_avg", show_AIC = FALSE)

###### Fig2 @ ######

# Arrange plots for the three selected environmental predictors
combined_OVERALL_envir_importantvar <- ggpubr::ggarrange(
  best_soilMoist_plot,
  best_precip_plot,
  best_soilTemp_plot,
  nrow = 2, ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

# Display the combined figure
combined_OVERALL_envir_importantvar

# Export the figure
ggsave(filename = "Figure_2.png", plot = combined_OVERALL_envir_importantvar,
  width = 7, height = 5, dpi = 1000)



##### 5.6. Multiple regression using selected predictors ######

# This section evaluates the combined effects of selected environmental predictors on total 
# biofilm-associated gene abundance. It compares additive and two-way interaction mixed-effects models using AIC.

# Check correlation among selected predictors
vars_total_full_model <- c(
  "soilMoist_avg_dayMin",
  "secPrecip_cv_dayTot",
  "soilTemp_avg",
  "airTemp_avg_dayMax_fall")

cor_vars_total_full_model <- cor(
  biofilm_env_samples[, vars_total_full_model],
  use = "complete.obs")

round(cor_vars_total_full_model, 2)

# Compare AIC of selected predictors
best_per_type_total %>% select(predictor, best_model, AIC_selected)

# Exclude highly correlated predictors (>= 0.7), retaining the one with the lowest individual AIC:
# soilTemp_avg and airTemp_avg_dayMax_fall are highly correlated.
# soilTemp_avg has the lower AIC and is therefore retained.


###### 5.6.1. Additive model ####

# Centered predictors were generated in Section 5.1.
# Missing soilMoist_avg_dayMin values substantially reduced sample size and coverage 
# across environmental gradients. Therefore, multiple regression are restricted to
# soilTemp_avg and secPrecip_cv_dayTot.


multiple_model_overall <- lmer(
  total_count ~
    secPrecip_cv_dayTot_c + I(secPrecip_cv_dayTot_c^2) +
    soilTemp_avg_c + I(soilTemp_avg_c^2) +
    (1 | siteID),
  data = biofilm_env_samples,
  REML = FALSE)

summary(multiple_model_overall)
anova(multiple_model_overall)
performance::r2(multiple_model_overall)


###### 5.6.2. Two-way interaction model ######

multiple_model_overall_twoway <- lmer(
  total_count ~
    secPrecip_cv_dayTot_c + I(secPrecip_cv_dayTot_c^2) +
    soilTemp_avg_c + I(soilTemp_avg_c^2) +
    secPrecip_cv_dayTot_c:soilTemp_avg_c +
    (1 | siteID),
  data = biofilm_env_samples,
  REML = FALSE)

summary(multiple_model_overall_twoway)
anova(multiple_model_overall_twoway)
performance::r2(multiple_model_overall_twoway)


###### 5.6.3. Comparing models (AIC) ######

AIC(multiple_model_overall, multiple_model_overall_twoway)
# The two-way interaction model had the lowest AIC compared with the aditive model




#### ! ####
#### 6. PATHWAY-LEVEL ####

# This section evaluates the relationships between environmental conditions and biofilm-associated gene abundance at 
# the pathway level.

# Sections 6.2–6.4 fit individual quadratic mixed-effects models to assess environmental responses within each pathway.
# Section 6.5 fits unified mixed-effects models to compare responses among pathways.

# Section 6.6 evaluates the combined effects of selected environmental predictors on gene abundance within each pathway 
# using multiple regression models. Additive and two-way interaction models are fitted and compared using AIC.


##### 6.1. Prepare pathway-level data #####

# Merge environmental variables with pathway-level abundances
env_vars <- biofilm_env_samples %>%
  select(-total_count)

biofilm_env_pathway_samples <- pathway_long %>%
  left_join(env_vars, by = c("genomicsSampleID", "siteID"))
# biofilm_env_pathway_samples is 8x the size of biofilm_env_samples.
# No need to summarize as in biofilm_env_samples, since pathway_long and biofilm_env_pathway_samples already 
# contain the total counts per pathway. 

# Master dataset for pathway-level analyses
biofilm_env_pathway_samples_all <- biofilm_env_pathway_samples

# Three environmental predictors retained from the overall analysis
predictors_pathway <- c(
  "soilMoist_avg_dayMin",
  "secPrecip_cv_dayTot",
  "soilTemp_avg")

# Identify biofilm-associated pathways
pathways_vec <- levels(biofilm_env_pathway_samples_all$biofilm_pathways_clean)

# Check pathway-level observations
biofilm_env_pathway_samples_all %>%
  dplyr::count(biofilm_pathways_clean)



##### 6.2. Quadratic mixed-effects models by pathway #####

# Fit individual quadratic mixed-effects models for each pathway and environmental
# predictor, with siteID as a random intercept to account for non-independence
# of samples within sites.

# Store fitted models and statistical results
pathway_models <- list()
pathway_results_list <- list()

# Fit one quadratic model for each pathway and predictor
for (pw in pathways_vec) {
  # Subset data for the current pathway
  dat_pw <- biofilm_env_pathway_samples_all %>%
    dplyr::filter(biofilm_pathways_clean == pw)
  pathway_models[[pw]] <- list()
  for (v in predictors_pathway) {
    # Center predictor within the current pathway
    v_c <- paste0(v, "_c")
    dat_current <- dat_pw %>%
      dplyr::filter(
        !is.na(.data[[v]]),
        !is.na(total_count_sample),
        !is.na(siteID))
    dat_current[[v_c]] <- as.numeric(scale(
      dat_current[[v]],
      center = TRUE,
      scale = FALSE))
    # Quadratic mixed-effects model
    formula_quad <- as.formula(paste0(
      "total_count_sample ~ ",
      v_c, " + I(", v_c, "^2) + (1 | siteID)"))
    model_quad <- lmerTest::lmer(
      formula_quad,
      data = dat_current,
      REML = FALSE)
    # Store fitted model
    pathway_models[[pw]][[v]] <- model_quad
    # Extract fixed-effect statistics
    model_summary <- broom.mixed::tidy(
      model_quad,
      effects = "fixed")
    # Extract marginal R2
    r2_marginal <- as.numeric(
      performance::r2(model_quad)$R2_marginal)
    # Identify linear and quadratic terms
    linear_term <- v_c
    quadratic_term <- paste0("I(", v_c, "^2)")
    linear_results <- model_summary %>%
      dplyr::filter(term == linear_term)
    quadratic_results <- model_summary %>%
      dplyr::filter(term == quadratic_term)
    # Store results for this pathway and predictor
    pathway_results_list[[paste(pw, v, sep = "_")]] <-
      tibble::tibble(
        pathway = pw,
        predictor = v,
        n_samples = nrow(dat_current),
        coef_linear = linear_results$estimate,
        se_linear = linear_results$std.error,
        p_linear = linear_results$p.value,
        coef_quadratic = quadratic_results$estimate,
        se_quadratic = quadratic_results$std.error,
        p_quadratic = quadratic_results$p.value,
        R2_marginal = r2_marginal)
  }
}

# Combine statistical results
pathway_results <- dplyr::bind_rows(pathway_results_list)
# Display results
View(pathway_results)



##### 6.3. Prepare pathway-level figure annotations #####

# Format model statistics for figures
pathway_plot_results <- pathway_results %>%
  dplyr::mutate(
    coef_linear_label = round(coef_linear, 2),
    coef_quadratic_label = round(coef_quadratic, 2),
    p_linear_label = scales::pvalue(p_linear),
    p_quadratic_label = scales::pvalue(p_quadratic),
    R2_label = dplyr::if_else(
      R2_marginal < 0.001,
      "<0.001",
      as.character(round(R2_marginal, 3))),
    # Flag significant relationships
    significant = !is.na(p_quadratic) & p_quadratic < 0.05,
    # Figure annotations
    model_label = paste0(
      "Linear: ", coef_linear_label,
      " (P ", ifelse(p_linear < 0.001, "", "= "), p_linear_label, ")\n",
      "Quadratic: ", coef_quadratic_label,
      " (P ", ifelse(p_quadratic < 0.001, "", "= "), p_quadratic_label, ")\n",
      "R²: ", R2_label))

# Inspect results
View(pathway_plot_results)
pathway_plot_results %>%
  dplyr::select(pathway, predictor, p_linear, p_quadratic, significant, model_label) %>%
  print(n = Inf)





##### 6.4. Generate individual pathway-level plots #####

# Original colors for environmental predictors
pathway_fill_colors <- c(
  soilMoist_avg_dayMin = "#298c8c",
  secPrecip_cv_dayTot = "#298c8c",
  soilTemp_avg = "#f1a226")

# Original environmental axis labels
pathway_x_labels <- c(
  soilMoist_avg_dayMin = "Soil Moisture Day Min (%)",
  secPrecip_cv_dayTot = "Precipitation Day Total CV (%)",
  soilTemp_avg = "Soil Temperature (°C)")

# Function to generate individual pathway-level plots
plot_pathway <- function(pw, v) {
  # Subset observations for the current pathway and predictor
  dat_plot <- biofilm_env_pathway_samples_all %>%
    dplyr::filter(
      biofilm_pathways_clean == pw,
      !is.na(.data[[v]]),
      !is.na(total_count_sample),
      !is.na(siteID))
  # Retrieve model statistics and figure annotation
  plot_info <- pathway_plot_results %>%
    dplyr::filter(
      pathway == pw,
      predictor == v)
  # Create scatterplot
  p <- ggplot2::ggplot(
    dat_plot,
    ggplot2::aes(
      x = .data[[v]],
      y = total_count_sample))
  # Add fitted quadratic curve only for significant quadratic relationships
  if (plot_info$significant) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x + I(x^2),
        se = TRUE,
        color = "black",
        fill = pathway_fill_colors[[v]],
        linewidth = 0.7,
        alpha = 0.4)}
  # Add original points, model annotation, and formatting
  p <- p +
    ggplot2::geom_point(size = 1.5, alpha = 1, shape = 21,
      color = "grey10", fill = "gray50", stroke = 0.35) +
    ggplot2::annotate("label", x = Inf, y = -Inf, label = plot_info$model_label,
      hjust = 1, vjust = -4.2, size = 2.3, linewidth = 0,
      fill = scales::alpha("white", 0), color = "black") +
    ggplot2::labs(x = pathway_x_labels[[v]],
      y = paste0(
        pw,
        "\nTotal Gene Abundance",
        "\n(counts / kb / million reads)")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(
        size = 9,
        face = "bold"),
      axis.text = ggplot2::element_text(size = 7))
  return(p)}

# Generate 24 individual plots using their original object names
for (pw in pathways_vec) {
  
  # Convert pathway names to original naming convention
  pw_name <- gsub("[- ]", "_", pw)
  
  for (v in predictors_pathway) {
    
    # Construct original plot object name
    plot_name <- paste0(
      "pathway_counts_",
      pw_name,
      "_",
      v,
      "_plot"
    )
    
    # Generate and store individual plot
    assign(
      plot_name,
      plot_pathway(pw, v)
    )
  }
}

# List generated plot objects
ls(pattern = "^pathway_counts_.*_plot$")


###### ABC transporters FigS6@ ######

combined_ABC_transporters_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_ABC_Transporters_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 50000, 10000)) +
    coord_cartesian(ylim = c(0, 50500)),
  
  pathway_counts_ABC_Transporters_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 50000, 10000)) +
    coord_cartesian(ylim = c(0, 50500)),
  
  pathway_counts_ABC_Transporters_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 50000, 10000)) +
    coord_cartesian(ylim = c(0, 50500)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S6.png",
       plot = combined_ABC_transporters_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### Two-Component System FigS7@ ######

combined_Two_Component_System_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_Two_Component_System_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 80000, 20000)) +
    coord_cartesian(ylim = c(0, 82500)),
  
  pathway_counts_Two_Component_System_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 80000, 20000)) +
    coord_cartesian(ylim = c(0, 82500)),
  
  pathway_counts_Two_Component_System_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 80000, 20000)) +
    coord_cartesian(ylim = c(0, 82500)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S7.png",
       plot = combined_Two_Component_System_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### Biofilm formation FigS8@ ######

combined_Biofilm_Formation_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_Biofilm_Formation_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 15000, 3000)) +
    coord_cartesian(ylim = c(0, 15000)),
  
  pathway_counts_Biofilm_Formation_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 15000, 3000)) +
    coord_cartesian(ylim = c(0, 15000)),
  
  pathway_counts_Biofilm_Formation_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 15000, 3000)) +
    coord_cartesian(ylim = c(0, 15000)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S8.png",
       plot = combined_Biofilm_Formation_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### EPS biosynthesis FigS9@ ######

combined_EPS_Biosynthesis_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_EPS_Biosynthesis_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 2000, 500)) +
    coord_cartesian(ylim = c(0, 2200)),
  
  pathway_counts_EPS_Biosynthesis_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 2000, 500)) +
    coord_cartesian(ylim = c(0, 2200)),
  
  pathway_counts_EPS_Biosynthesis_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 2000, 500)) +
    coord_cartesian(ylim = c(0, 2200)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S9.png",
       plot = combined_EPS_Biosynthesis_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### Quorum sensing FigS10@ ######


combined_Quorum_Sensing_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_Quorum_Sensing_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 40000, 10000)) +
    coord_cartesian(ylim = c(0, 42500)),
  
  pathway_counts_Quorum_Sensing_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 40000, 10000)) +
    coord_cartesian(ylim = c(0, 42500)),
  
  pathway_counts_Quorum_Sensing_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 40000, 10000)) +
    coord_cartesian(ylim = c(0, 42500)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S10.png",
       plot = combined_Quorum_Sensing_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### Chemotaxis FigS11@ ######

combined_Chemotaxis_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_Chemotaxis_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 27000, 5000)) +
    coord_cartesian(ylim = c(0, 26500)),
  
  pathway_counts_Chemotaxis_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 27000, 5000)) +
    coord_cartesian(ylim = c(0, 26500)),
  
  pathway_counts_Chemotaxis_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 27000, 5000)) +
    coord_cartesian(ylim = c(0, 26500)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S11.png",
       plot = combined_Chemotaxis_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### Flagellar assembly FigS12@ ######

combined_Flagellar_Assembly_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_Flagellar_Assembly_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 9000, 2000)) +
    coord_cartesian(ylim = c(0, 8750)),
  
  pathway_counts_Flagellar_Assembly_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 9000, 2000)) +
    coord_cartesian(ylim = c(0, 8750)),
  
  pathway_counts_Flagellar_Assembly_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 9000, 2000)) +
    coord_cartesian(ylim = c(0, 8750)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S12.png",
       plot = combined_Flagellar_Assembly_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)


###### BBSdb FigS13@ ######

combined_BBSdb_envir_importantvar_SHARED <- ggarrange(
  pathway_counts_BBSdb_soilMoist_avg_dayMin_plot +
    scale_y_continuous(breaks = seq(0, 180000, 40000)) +
    coord_cartesian(ylim = c(0, 173500)),
  
  pathway_counts_BBSdb_secPrecip_cv_dayTot_plot +
    scale_y_continuous(breaks = seq(0, 180000, 40000)) +
    coord_cartesian(ylim = c(0, 173500)),
  
  pathway_counts_BBSdb_soilTemp_avg_plot +
    scale_y_continuous(breaks = seq(0, 180000, 40000)) +
    coord_cartesian(ylim = c(0, 173500)),
  
  nrow = 2,
  ncol = 2,
  labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

ggsave(filename = "Figure_S13.png",
       plot = combined_BBSdb_envir_importantvar_SHARED,
       width = 7, height = 5, dpi = 1000)



##### 6.5. Comparing slopes ######

# This section compares environmental responses among biofilm-associated pathways using unified mixed-effects models
# with pathway identity and its interactions with each environmental predictor. 
# Site is included as a random effect to account for non-independence of samples within sites, 
# while sample identity is included as another random effect to account for non-independence of counts contributing to
# multiple pathways. 
# BBSdb is excluded because it represents a curated database rather than a specific biofilm-associated pathway.


###### 6.5.1. Prepare data for pathway comparisons ######

# Exclude BBSdb from pathway comparisons
biofilm_env_pathway_samples_nobbsdb <- biofilm_env_pathway_samples_all %>%
  dplyr::filter(biofilm_pathways_clean != "BBSdb") %>%
  droplevels()

# Define colors without BBSdb
pathway_colors_nobbsdb <- pathway_colors[names(pathway_colors) != "BBSdb"]


###### ...soilMoist_avg_dayMin  Fig3A@ ######

# Fit unified quadratic interaction model
model_soilMoist_avg_dayMin_pathway_poly_nobbsdb <- lmer(
  total_count_sample ~
    soilMoist_avg_dayMin_c +
    I(soilMoist_avg_dayMin_c^2) +
    biofilm_pathways_clean +
    soilMoist_avg_dayMin_c:biofilm_pathways_clean +
    I(soilMoist_avg_dayMin_c^2):biofilm_pathways_clean +
    (1 | siteID) +
    (1 | genomicsSampleID),
  data = biofilm_env_pathway_samples_nobbsdb,
  REML = FALSE)

# Results
summary(model_soilMoist_avg_dayMin_pathway_poly_nobbsdb)
anova(model_soilMoist_avg_dayMin_pathway_poly_nobbsdb)
performance::r2(model_soilMoist_avg_dayMin_pathway_poly_nobbsdb)

# Extract interaction P-values for plot annotation
anova_soilMoist_avg_dayMin_nobbsdb <- anova(model_soilMoist_avg_dayMin_pathway_poly_nobbsdb)
p_linear_soilMoist_avg_dayMin_nobbsdb <- anova_soilMoist_avg_dayMin_nobbsdb["soilMoist_avg_dayMin_c:biofilm_pathways_clean", "Pr(>F)"]
p_quad_soilMoist_avg_dayMin_nobbsdb <- anova_soilMoist_avg_dayMin_nobbsdb["I(soilMoist_avg_dayMin_c^2):biofilm_pathways_clean", "Pr(>F)"]
p_linear_lab_soilMoist_avg_dayMin_nobbsdb <- scales::pvalue(p_linear_soilMoist_avg_dayMin_nobbsdb)
p_quad_lab_soilMoist_avg_dayMin_nobbsdb   <- scales::pvalue(p_quad_soilMoist_avg_dayMin_nobbsdb)

interaction_label_soilMoist_avg_dayMin_nobbsdb <- paste0(
  "Pathway × soil moisture interaction:\n",
  "Linear: P = ", p_linear_lab_soilMoist_avg_dayMin_nobbsdb, "\n",
  "Quadratic: P ", p_quad_lab_soilMoist_avg_dayMin_nobbsdb)

# Define significant pathways for soilMoist_avg_dayMin
# Check in pathway_plot_results from section 6.3
significant_pathways_soilMoist_avg_dayMin_nobbsdb <- c(
  "ABC Transporters",
  "Quorum Sensing")

# Split data: all vs significant subset
biofilm_env_pathway_samples_soilMoist_avg_dayMin_sig_nobbsdb <- biofilm_env_pathway_samples_nobbsdb %>%
  dplyr::filter(biofilm_pathways_clean %in% significant_pathways_soilMoist_avg_dayMin_nobbsdb)

# Plot
p_soilMoist_avg_dayMin_curves_nobbsdb <- ggplot(
  biofilm_env_pathway_samples_nobbsdb,
  aes(x = soilMoist_avg_dayMin, y = total_count_sample,
      color = biofilm_pathways_clean, fill = biofilm_pathways_clean)) +
  geom_point(aes(fill = biofilm_pathways_clean),
             color = "black", size = 1.2, alpha = 1, shape = 21, stroke = 0.25) +
  geom_smooth(
    data = biofilm_env_pathway_samples_soilMoist_avg_dayMin_sig_nobbsdb,
    method = "lm", formula = y ~ x + I(x^2),
    se = TRUE, linewidth = 0.7, alpha = 0.4) +
  annotate("label", x = Inf, y = -Inf,
           label = interaction_label_soilMoist_avg_dayMin_nobbsdb,
           hjust = 1.05, vjust = -4.7, size = 2.3, linewidth = 0,
           fill = scales::alpha("white", 0), color = "black") +
  scale_color_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  scale_fill_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  labs(x = "Soil Moisture Day Min (%)",
       y = "Total Gene Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"))

plot(p_soilMoist_avg_dayMin_curves_nobbsdb)


###### ...secPrecip_cv_dayTot Fig3B@ ######

# Fit unified quadratic interaction model
model_secPrecip_cv_dayTot_pathway_poly_nobbsdb <- lmer(
  total_count_sample ~
    secPrecip_cv_dayTot_c +
    I(secPrecip_cv_dayTot_c^2) +
    biofilm_pathways_clean +
    secPrecip_cv_dayTot_c:biofilm_pathways_clean +
    I(secPrecip_cv_dayTot_c^2):biofilm_pathways_clean +
    (1 | siteID) +
    (1 | genomicsSampleID),
  data = biofilm_env_pathway_samples_nobbsdb,
  REML = FALSE)

# Results
summary(model_secPrecip_cv_dayTot_pathway_poly_nobbsdb)
anova(model_secPrecip_cv_dayTot_pathway_poly_nobbsdb)
performance::r2(model_secPrecip_cv_dayTot_pathway_poly_nobbsdb)

# Extract interaction P-values for plot annotation
anova_secPrecip_cv_dayTot_nobbsdb <- anova(model_secPrecip_cv_dayTot_pathway_poly_nobbsdb)
p_linear_secPrecip_cv_dayTot_nobbsdb <- anova_secPrecip_cv_dayTot_nobbsdb["secPrecip_cv_dayTot_c:biofilm_pathways_clean", "Pr(>F)"]
p_quad_secPrecip_cv_dayTot_nobbsdb   <- anova_secPrecip_cv_dayTot_nobbsdb["I(secPrecip_cv_dayTot_c^2):biofilm_pathways_clean", "Pr(>F)"]
p_linear_lab_secPrecip_cv_dayTot_nobbsdb <- scales::pvalue(p_linear_secPrecip_cv_dayTot_nobbsdb)
p_quad_lab_secPrecip_cv_dayTot_nobbsdb   <- scales::pvalue(p_quad_secPrecip_cv_dayTot_nobbsdb)

interaction_label_secPrecip_cv_dayTot_nobbsdb <- paste0(
  "Pathway × precipitation interaction:\n",
  "Linear: P ", p_linear_lab_secPrecip_cv_dayTot_nobbsdb, "\n",
  "Quadratic: P ", p_quad_lab_secPrecip_cv_dayTot_nobbsdb)

# Define significant pathways for secPrecip_cv_dayTot
# Check in pathway_plot_results from section 6.3
significant_pathways_secPrecip_cv_dayTot_nobbsdb <- c(
  "ABC Transporters",
  "Flagellar Assembly",
  "Quorum Sensing",
  "Two-Component System")

# Subset significant pathways for plotting
biofilm_env_pathway_samples_secPrecip_cv_dayTot_sig_nobbsdb <- biofilm_env_pathway_samples_nobbsdb %>%
  dplyr::filter(biofilm_pathways_clean %in% significant_pathways_secPrecip_cv_dayTot_nobbsdb)

# Plot
p_secPrecip_cv_dayTot_curves_nobbsdb <- ggplot(
  biofilm_env_pathway_samples_nobbsdb,
  aes(x = secPrecip_cv_dayTot, y = total_count_sample,
      color = biofilm_pathways_clean, fill = biofilm_pathways_clean)) +
  geom_point(aes(fill = biofilm_pathways_clean),
             color = "black", size = 1.2, alpha = 1, shape = 21, stroke = 0.25) +
  geom_smooth(
    data = biofilm_env_pathway_samples_secPrecip_cv_dayTot_sig_nobbsdb,
    method = "lm", formula = y ~ x + I(x^2),
    se = TRUE, linewidth = 0.7, alpha = 0.4) +
  annotate("label", x = Inf, y = -Inf,
           label = interaction_label_secPrecip_cv_dayTot_nobbsdb,
           hjust = 1.05, vjust = -4.7, size = 2.3, linewidth = 0,
           fill = scales::alpha("white", 0), color = "black") +
  scale_color_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  scale_fill_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  labs(x = "Precipitation Day Tot CV (%)",
       y = "Total Gene Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"))

plot(p_secPrecip_cv_dayTot_curves_nobbsdb)


###### ...soilTemp_avg Fig3C@ ######

# Fit unified quadratic interaction model
model_soilTemp_avg_pathway_poly_nobbsdb <- lmer(
  total_count_sample ~
    soilTemp_avg_c +
    I(soilTemp_avg_c^2) +
    biofilm_pathways_clean +
    soilTemp_avg_c:biofilm_pathways_clean +
    I(soilTemp_avg_c^2):biofilm_pathways_clean +
    (1 | siteID) +
    (1 | genomicsSampleID),
  data = biofilm_env_pathway_samples_nobbsdb,
  REML = FALSE)

# Results
summary(model_soilTemp_avg_pathway_poly_nobbsdb)
anova(model_soilTemp_avg_pathway_poly_nobbsdb)
performance::r2(model_soilTemp_avg_pathway_poly_nobbsdb)

# Extract interaction P-values for plot annotation
anova_soilTemp_avg_nobbsdb <- anova(model_soilTemp_avg_pathway_poly_nobbsdb)
p_linear_soilTemp_avg_nobbsdb <- anova_soilTemp_avg_nobbsdb["soilTemp_avg_c:biofilm_pathways_clean", "Pr(>F)"]
p_quad_soilTemp_avg_nobbsdb   <- anova_soilTemp_avg_nobbsdb["I(soilTemp_avg_c^2):biofilm_pathways_clean", "Pr(>F)"]
p_linear_lab_soilTemp_avg_nobbsdb <- scales::pvalue(p_linear_soilTemp_avg_nobbsdb)
p_quad_lab_soilTemp_avg_nobbsdb   <- scales::pvalue(p_quad_soilTemp_avg_nobbsdb)

interaction_label_soilTemp_avg_nobbsdb <- paste0(
  "Pathway × temperature interaction:\n",
  "Linear: P ", p_linear_lab_soilTemp_avg_nobbsdb, "\n",
  "Quadratic: P ", p_quad_lab_soilTemp_avg_nobbsdb)

# Define significant pathways for soilTemp_avg
# Check in pathway_plot_results from section 6.3
significant_pathways_soilTemp_avg_nobbsdb <- c(
  "Quorum Sensing",
  "Two-Component System")

# Subset significant pathways for plotting
biofilm_env_pathway_samples_soilTemp_avg_sig_nobbsdb <- biofilm_env_pathway_samples_nobbsdb %>%
  dplyr::filter(biofilm_pathways_clean %in% significant_pathways_soilTemp_avg_nobbsdb)

# Plot
p_soilTemp_avg_curves_nobbsdb <- ggplot(
  biofilm_env_pathway_samples_nobbsdb,
  aes(x = soilTemp_avg, y = total_count_sample,
      color = biofilm_pathways_clean, fill = biofilm_pathways_clean)) +
  geom_point(aes(fill = biofilm_pathways_clean),
             color = "black", size = 1.2, alpha = 1, shape = 21, stroke = 0.25) +
  geom_smooth(
    data = biofilm_env_pathway_samples_soilTemp_avg_sig_nobbsdb,
    method = "lm", formula = y ~ x + I(x^2),
    se = TRUE, linewidth = 0.7, alpha = 0.4) +
  annotate("label", x = Inf, y = -Inf,
           label = interaction_label_soilTemp_avg_nobbsdb,
           hjust = 1.05, vjust = -4.7, size = 2.3, linewidth = 0,
           fill = scales::alpha("white", 0), color = "black") +
  scale_color_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  scale_fill_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  labs(x = "Soil Temperature (°C)",
       y = "Total Gene Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"))

plot(p_soilTemp_avg_curves_nobbsdb)



###### Fig3 @ ######

# Create a dummy plot to extract the shared legend (it excludes BBSdb)
legend_plot_nobbsdb <- ggplot(
  biofilm_env_pathway_samples_nobbsdb,
  aes(x = soilTemp_avg, y = total_count_sample,
      color = biofilm_pathways_clean, fill = biofilm_pathways_clean)) +
  geom_point(shape = 21, size = 2, stroke = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2),
              se = TRUE, linewidth = 0.7, alpha = 0.4) +
  scale_color_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  scale_fill_manual(values = pathway_colors_nobbsdb, name = "Pathway") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"))

# Extract the shared legend
leg_full_nobbsdb <- get_legend(legend_plot_nobbsdb)

# Remove legends from individual plots
p1_nobbsdb <- p_soilMoist_avg_dayMin_curves_nobbsdb + theme(legend.position = "none")
p2_nobbsdb <- p_secPrecip_cv_dayTot_curves_nobbsdb  + theme(legend.position = "none")
p3_nobbsdb <- p_soilTemp_avg_curves_nobbsdb         + theme(legend.position = "none")


# Arrange panels
combined_ALL_PATHWAYS_envir_importantvar_SHARED_MULTIPLE_panels_nobbsdb <- ggarrange(
  p1_nobbsdb, p2_nobbsdb, p3_nobbsdb,
  nrow = 2, ncol = 2, labels = LETTERS,
  font.label = list(size = 11, face = "bold"))

# Add shared legend
Figure3_nobbsdb <- ggarrange(
  combined_ALL_PATHWAYS_envir_importantvar_SHARED_MULTIPLE_panels_nobbsdb,
  leg_full_nobbsdb,
  ncol = 1,
  heights = c(1, 0.12))

ggsave(filename = "Figure_3.png", plot = Figure3_nobbsdb,
       width = 7, height = 6, dpi = 1000)


##### 6.6 Multiple regression per pathway @@@ ######

# This section evaluates the combined effects of precipitation variability and soil temperature on
# biofilm-associated gene abundance per pathway.
# Additive and interaction mixed-effects models are compared using AIC.
# Site is included as a random intercept to account for non-independence of samples within sites. 
# The model with the lowest AIC is selected for each pathway.

# Split dataset by pathway
pathway_list <- biofilm_env_pathway_samples_all %>%
  split(.$biofilm_pathways_clean)

# Function to fit and compare models within each pathway
fit_pathway_summary <- function(df, pathway_name) {
  # Remove observations with missing values
  df <- df %>%
    filter(
      !is.na(total_count_sample),
      !is.na(secPrecip_cv_dayTot),
      !is.na(soilTemp_avg),
      !is.na(siteID))
  if (nrow(df) < 10) return(NULL)
  # Center predictors within pathway
  df <- df %>%
    mutate(
      secPrecip_cv_dayTot_c = as.numeric(scale(secPrecip_cv_dayTot, center = TRUE, scale = FALSE)),
      soilTemp_avg_c        = as.numeric(scale(soilTemp_avg, center = TRUE, scale = FALSE)))
  # Fit additive model
  mod_add <- lmer(
    total_count_sample ~
      secPrecip_cv_dayTot_c + I(secPrecip_cv_dayTot_c^2) +
      soilTemp_avg_c + I(soilTemp_avg_c^2) +
      (1 | siteID),
    data = df,
    REML = FALSE)
  # Fit two-way interaction model
  mod_int <- lmer(
    total_count_sample ~
      secPrecip_cv_dayTot_c + I(secPrecip_cv_dayTot_c^2) +
      soilTemp_avg_c + I(soilTemp_avg_c^2) +
      secPrecip_cv_dayTot_c:soilTemp_avg_c +
      (1 | siteID),
    data = df,
    REML = FALSE)
  # Compare models using AIC and select the model with the lowest AIC
  aic_add <- AIC(mod_add)
  aic_int <- AIC(mod_int)
  if (aic_add <= aic_int) {
    best_model <- mod_add
    best_model_name <- "Additive"
    best_aic <- aic_add} else {
    best_model <- mod_int
    best_model_name <- "Interaction"
    best_aic <- aic_int}
  # Extract marginal R2 from the selected model
  r2_best <- performance::r2(best_model)$R2_marginal
  
  # Extract results for Table1
  
  # Extract fixed effects from selected model
  tidy_best <- broom.mixed::tidy(best_model, effects = "fixed") %>%
    filter(term != "(Intercept)")
  # Keep only significant terms and relabel them
  sig_terms <- tidy_best %>%
    filter(p.value < 0.05) %>%
    mutate(
      predictor = case_when(
        term == "secPrecip_cv_dayTot_c" ~ "Precipitation CV",
        term == "I(secPrecip_cv_dayTot_c^2)" ~ "Precipitation CV",
        term == "soilTemp_avg_c" ~ "Soil temperature",
        term == "I(soilTemp_avg_c^2)" ~ "Soil temperature",
        term == "secPrecip_cv_dayTot_c:soilTemp_avg_c" ~ "Precipitation CV × Soil temperature"),
      effect_type = case_when(
        term == "secPrecip_cv_dayTot_c" ~ "Linear",
        term == "I(secPrecip_cv_dayTot_c^2)" ~ "Quadratic",
        term == "soilTemp_avg_c" ~ "Linear",
        term == "I(soilTemp_avg_c^2)" ~ "Quadratic",
        term == "secPrecip_cv_dayTot_c:soilTemp_avg_c" ~ "Interaction"),
      term_text = paste0(
        predictor,
        " (", effect_type, " = ", round(estimate, 3),
        ", P = ", format.pval(p.value, digits = 2, eps = 0.001), ")"))
  # Combine significant terms into a single summary
  sig_terms_text <- if (nrow(sig_terms) == 0) {
    "None"
  } else {
    paste(sig_terms$term_text, collapse = "; ")}
  # Return model selection and statistical results for each pathway
  tibble(
    `Pathway` = pathway_name,
    `Best model` = best_model_name,
    `AIC` = round(best_aic, 3),
    `Significant terms` = sig_terms_text,
    `Marginal R2` = round(r2_best, 3))}



##### Table 1 @ ####

# Apply function to all pathways
pathway_summary_table <- imap_dfr(pathway_list, fit_pathway_summary)

# View final table
view(pathway_summary_table)

# Save table and edit manually
write.csv(pathway_summary_table, "Table1.csv", row.names = FALSE)




####!####
#### 7. KO-LEVEL ####

# This section evaluates the relationships between environmental conditions and the abundance of individual 
# KOs associated with biofilm-related pathways.
# Quadratic mixed-effects models are fitted for the three selected environmental predictors, 
# and representative KOs are selected for visualization.


##### 7.1. Prepare KO-level data #####

# Convert KO abundance data from wide to long format
anno_cols <- c("biofilm_pathways", "biofilm_KOs",
               "ann_preferred", "ann_desc",
               "n_orfs", "n_contigs")

sample_cols <- setdiff(names(pathway_ko_summary), anno_cols)

biofilm_long_pathwayKO <- pathway_ko_summary %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "genomicsSampleID",
    values_to = "count")

# Merge sample metadata
biofilm_long_pathwayKO <- biofilm_long_pathwayKO %>%
  left_join(metadata_dedup, by = "genomicsSampleID")

# Merge environmental data
biofilm_long_pathwayKO_env <- biofilm_long_pathwayKO %>%
  left_join(
    biofilm_env_samples %>% 
      select(-siteID),  
    by = "genomicsSampleID")

# Standardize pathway names using labels defined in Section 4.1
biofilm_long_pathwayKO_env <- biofilm_long_pathwayKO_env %>%
  mutate(
    biofilm_pathways_clean = dplyr::recode(
      biofilm_pathways,
      !!!pathway_labels,
      .default = biofilm_pathways))

# Set pathway order for figures
biofilm_long_pathwayKO_env <- biofilm_long_pathwayKO_env %>%
  mutate(
    biofilm_pathways_clean = factor(
      biofilm_pathways_clean,
      levels = c(
        "ABC Transporters",
        "Two-Component System",
        "Biofilm Formation",
        "EPS Biosynthesis",
        "Quorum Sensing",
        "Chemotaxis",
        "Flagellar Assembly",
        "BBSdb")))



##### 7.2. Select the 20 most abundant KOs per pathway #####

# Calculate mean KO abundance across samples within each pathway
top_kos_by_pathway <- biofilm_long_pathwayKO_env %>%
  group_by(biofilm_pathways, biofilm_KOs) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop") %>%
  group_by(biofilm_pathways) %>%
  slice_max(mean_count, n = 20, with_ties = FALSE) %>%
  ungroup()

# Retain the selected KOs and their corresponding sample-level data
biofilm_long_pathwayKO_env_top <- biofilm_long_pathwayKO_env %>%
  semi_join(top_kos_by_pathway, by = c("biofilm_pathways", "biofilm_KOs"))



##### 7.3. Quadratic mixed-effects models by KO #####


###### 7.3.1. soilMoist_avg_dayMin #####


# Fit a quadratic mixed-effects model for each KO
fit_KO_soilMoist <- function(df){
  # Center soil moisture within each KO subset
  df <- df %>% mutate(soilMoist_avg_dayMin_c = as.numeric(scale(soilMoist_avg_dayMin, center = TRUE, scale = FALSE)))
  # Fit quadratic model with site as a random intercept
  model_quad <- lmer(
    count ~ soilMoist_avg_dayMin_c + I(soilMoist_avg_dayMin_c^2) + (1 | siteID),
    data = df, REML = FALSE)
  # Extract model statistics
  a_quad <- anova(model_quad)
  coefs <- summary(model_quad)$coefficients
  r2_vals <- tryCatch(performance::r2(model_quad), error = function(e) NULL)
  r2_marg <- if (!is.null(r2_vals)) r2_vals$R2_marginal else NA_real_
  r2_cond <- if (!is.null(r2_vals)) r2_vals$R2_conditional else NA_real_
  # Return model results
  tibble(
    biofilm_pathways = df$biofilm_pathways[1],
    biofilm_KOs      = df$biofilm_KOs[1],
    p_linear         = a_quad["soilMoist_avg_dayMin_c", "Pr(>F)"],
    p_quad           = a_quad["I(soilMoist_avg_dayMin_c^2)", "Pr(>F)"],
    AIC_quad         = AIC(model_quad),
    beta_linear      = coefs["soilMoist_avg_dayMin_c", "Estimate"],
    se_linear        = coefs["soilMoist_avg_dayMin_c", "Std. Error"],
    beta_quad        = coefs["I(soilMoist_avg_dayMin_c^2)", "Estimate"],
    se_quad          = coefs["I(soilMoist_avg_dayMin_c^2)", "Std. Error"],
    R2_marginal      = r2_marg,
    R2_conditional   = r2_cond)}

# Apply the model to each KO within each pathway
ko_soilMoist_results <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  group_modify(~fit_KO_soilMoist(.x)) %>%
  ungroup()
view(ko_soilMoist_results)

# Identify KOs with significant quadratic relationships
ko_soilMoist_sig <- ko_soilMoist_results %>%
  filter(p_quad < 0.05)
view(ko_soilMoist_sig)

# Plot

# Calculate mean KO abundance across samples within each pathway
ko_soilMoist_avg_dayMin_mean_abundance_pathway <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop")

# Select the most abundant KO with a significant quadratic response per pathway
ko_soilMoist_avg_dayMin_sig_abundant_pathway <- ko_soilMoist_results %>%
  filter(p_quad < 0.05) %>%
  left_join(
    ko_soilMoist_avg_dayMin_mean_abundance_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs")) %>%
  group_by(biofilm_pathways_clean) %>%
  slice_max(mean_count, n = 1, with_ties = FALSE) %>%
  ungroup()

# Prepare regression statistics for figure annotations
ko_soilMoist_avg_dayMin_sig_abundant_pathway <- ko_soilMoist_avg_dayMin_sig_abundant_pathway %>%
  mutate(
    coef_linear = round(beta_linear, 2),
    coef_quad   = round(beta_quad, 2),
    p_linear_label = ifelse(
      p_linear < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_linear))),
    p_quad_label = ifelse(
      p_quad < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_quad))),
    r2_label = ifelse(R2_marginal < 0.001, "<0.001", round(R2_marginal, 3)),
    lm_label = glue::glue(
      "Linear: {coef_linear} ({p_linear_label})\n",
      "Quadratic: {coef_quad} ({p_quad_label})\n",
      "R²: {r2_label}"))

# Keep sample-level data for the selected KOs
biofilm_long_pathwayKO_env_top_soilMoist_selected <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_soilMoist_avg_dayMin_sig_abundant_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs"))

# Arrange pathways and create facet labels
ko_soilMoist_avg_dayMin_sig_abundant_pathway <- ko_soilMoist_avg_dayMin_sig_abundant_pathway %>%
  arrange(biofilm_pathways_clean) %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_soilMoist_selected <- biofilm_long_pathwayKO_env_top_soilMoist_selected %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_soilMoist_avg_dayMin_sig_abundant_pathway$pathway_KO_label)))


###### ...FigS14 @ #####
KO_soilMoist_avg_dayMin_selected_plot <- ggplot(biofilm_long_pathwayKO_env_top_soilMoist_selected,
                                                aes(x = soilMoist_avg_dayMin, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean), method = "lm", formula = y ~ x + I(x^2),
              se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE,
              color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50", size = 1.5, stroke = 0.35) +
  geom_label(data = ko_soilMoist_avg_dayMin_sig_abundant_pathway,
             aes(x = Inf, y = -Inf, label = lm_label),
             inherit.aes = FALSE, hjust = 1, vjust = -4, size = 2.3,
             label.size = 0, fill = scales::alpha("white",0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(x = "Soil Moisture Day Min (%)", y = "Total KO Abundance (counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text  = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

plot(KO_soilMoist_avg_dayMin_selected_plot)

ggsave("Figure_S14.png", plot = KO_soilMoist_avg_dayMin_selected_plot, units = "in", height = 10, width = 7, dpi = 1000)



###### 7.3.2. secPrecip_cv_dayTot #####

# Fit a quadratic mixed-effects model for each KO
fit_KO_secPrecip <- function(df, keys){
  # Center precipitation variability within each KO subset
  df <- df %>% mutate(secPrecip_cv_dayTot_c = as.numeric(scale(secPrecip_cv_dayTot, center = TRUE, scale = FALSE)))
  # Fit quadratic model with site as a random intercept
  model_quad <- tryCatch(
    lmer(
      count ~ secPrecip_cv_dayTot_c + I(secPrecip_cv_dayTot_c^2) + (1 | siteID),
      data = df, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  # Return missing values if model fitting fails
  if (is.null(model_quad)) {
    return(tibble(
      p_linear       = NA_real_,
      p_quad         = NA_real_,
      AIC_quad       = NA_real_,
      beta_linear    = NA_real_,
      se_linear      = NA_real_,
      beta_quad      = NA_real_,
      se_quad        = NA_real_,
      R2_marginal    = NA_real_,
      R2_conditional = NA_real_))}
  # Extract model statistics
  a_quad <- tryCatch(anova(model_quad), error = function(e) NULL)
  coefs <- tryCatch(summary(model_quad)$coefficients, error = function(e) NULL)
  r2_vals <- tryCatch(performance::r2(model_quad), error = function(e) NULL)
  r2_marg <- if (!is.null(r2_vals)) r2_vals$R2_marginal else NA_real_
  r2_cond <- if (!is.null(r2_vals)) r2_vals$R2_conditional else NA_real_
  # Return model results
  tibble(
    p_linear       = if (!is.null(a_quad)) a_quad["secPrecip_cv_dayTot_c", "Pr(>F)"] else NA_real_,
    p_quad         = if (!is.null(a_quad)) a_quad["I(secPrecip_cv_dayTot_c^2)", "Pr(>F)"] else NA_real_,
    AIC_quad       = AIC(model_quad),
    beta_linear    = if (!is.null(coefs)) coefs["secPrecip_cv_dayTot_c", "Estimate"] else NA_real_,
    se_linear      = if (!is.null(coefs)) coefs["secPrecip_cv_dayTot_c", "Std. Error"] else NA_real_,
    beta_quad      = if (!is.null(coefs)) coefs["I(secPrecip_cv_dayTot_c^2)", "Estimate"] else NA_real_,
    se_quad        = if (!is.null(coefs)) coefs["I(secPrecip_cv_dayTot_c^2)", "Std. Error"] else NA_real_,
    R2_marginal    = r2_marg,
    R2_conditional = r2_cond)}

# Apply the model to each KO within each pathway
ko_secPrecip_results <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  group_modify(~fit_KO_secPrecip(.x, .y)) %>%
  ungroup()
view(ko_secPrecip_results)

# Identify KOs with significant quadratic relationships
ko_secPrecip_sig <- ko_secPrecip_results %>%
  filter(p_quad < 0.05)
view(ko_secPrecip_sig)

# Plot

# Calculate mean KO abundance across samples within each pathway
ko_secPrecip_cv_dayTot_mean_abundance_pathway <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop")

# Select the most abundant KO with a significant quadratic response per pathway
ko_secPrecip_cv_dayTot_sig_abundant_pathway <- ko_secPrecip_results %>%
  filter(p_quad < 0.05) %>%
  left_join(
    ko_secPrecip_cv_dayTot_mean_abundance_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs")) %>%
  group_by(biofilm_pathways_clean) %>%
  slice_max(mean_count, n = 1, with_ties = FALSE) %>%
  ungroup()

# Prepare regression statistics for figure annotations
ko_secPrecip_cv_dayTot_sig_abundant_pathway <- ko_secPrecip_cv_dayTot_sig_abundant_pathway %>%
  mutate(
    coef_linear = round(beta_linear, 2),
    coef_quad   = round(beta_quad, 2),
    p_linear_label = ifelse(
      p_linear < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_linear))),
    p_quad_label = ifelse(
      p_quad < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_quad))),
    r2_label = ifelse(R2_marginal < 0.001, "<0.001", round(R2_marginal, 3)),
    lm_label = glue::glue(
      "Linear: {coef_linear} ({p_linear_label})\n",
      "Quadratic: {coef_quad} ({p_quad_label})\n",
      "R²: {r2_label}"))

# Keep sample-level data for the selected KOs
biofilm_long_pathwayKO_env_top_secPrecip_selected <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_secPrecip_cv_dayTot_sig_abundant_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs"))

# Arrange pathways and create facet labels
ko_secPrecip_cv_dayTot_sig_abundant_pathway <- ko_secPrecip_cv_dayTot_sig_abundant_pathway %>%
  arrange(biofilm_pathways_clean) %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_secPrecip_selected <- biofilm_long_pathwayKO_env_top_secPrecip_selected %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_secPrecip_cv_dayTot_sig_abundant_pathway$pathway_KO_label)))


###### ...FigS15 @ #####
KO_secPrecip_cv_dayTot_selected_plot <- ggplot(biofilm_long_pathwayKO_env_top_secPrecip_selected,
                                               aes(x = secPrecip_cv_dayTot, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean), method = "lm", formula = y ~ x + I(x^2),
              se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE,
              color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50", size = 1.5, stroke = 0.35) +
  geom_label(data = ko_secPrecip_cv_dayTot_sig_abundant_pathway,
             aes(x = Inf, y = -Inf, label = lm_label),
             inherit.aes = FALSE, hjust = 1, vjust = -4, size = 2.3,
             label.size = 0, fill = scales::alpha("white",0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(x = "Precipitation Day Total CV (%)", y = "Total KO Abundance (counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text  = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

plot(KO_secPrecip_cv_dayTot_selected_plot)

ggsave("Figure_S15.png", plot = KO_secPrecip_cv_dayTot_selected_plot, units = "in", height = 10, width = 7, dpi = 1000)



###### 7.3.3. soilTemp_avg #####

# Fit a quadratic mixed-effects model for each KO
fit_KO_soilTemp <- function(df, keys){
  # Center soil temperature within each KO subset
  df <- df %>% mutate(soilTemp_avg_c = as.numeric(scale(soilTemp_avg, center = TRUE, scale = FALSE)))
  # Fit quadratic model with site as a random intercept
  model_quad <- tryCatch(
    lmer(
      count ~ soilTemp_avg_c + I(soilTemp_avg_c^2) + (1 | siteID),
      data = df, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  # Return missing values if model fitting fails
  if (is.null(model_quad)) {
    return(tibble(
      p_linear       = NA_real_,
      p_quad         = NA_real_,
      AIC_quad       = NA_real_,
      beta_linear    = NA_real_,
      se_linear      = NA_real_,
      beta_quad      = NA_real_,
      se_quad        = NA_real_,
      R2_marginal    = NA_real_,
      R2_conditional = NA_real_))}
  # Extract model statistics
  a_quad <- tryCatch(anova(model_quad), error = function(e) NULL)
  coefs <- tryCatch(summary(model_quad)$coefficients, error = function(e) NULL)
  r2_vals <- tryCatch(performance::r2(model_quad), error = function(e) NULL)
  r2_marg <- if (!is.null(r2_vals)) r2_vals$R2_marginal else NA_real_
  r2_cond <- if (!is.null(r2_vals)) r2_vals$R2_conditional else NA_real_
  # Return model results
  tibble(
    p_linear       = if (!is.null(a_quad)) a_quad["soilTemp_avg_c", "Pr(>F)"] else NA_real_,
    p_quad         = if (!is.null(a_quad)) a_quad["I(soilTemp_avg_c^2)", "Pr(>F)"] else NA_real_,
    AIC_quad       = AIC(model_quad),
    beta_linear    = if (!is.null(coefs)) coefs["soilTemp_avg_c", "Estimate"] else NA_real_,
    se_linear      = if (!is.null(coefs)) coefs["soilTemp_avg_c", "Std. Error"] else NA_real_,
    beta_quad      = if (!is.null(coefs)) coefs["I(soilTemp_avg_c^2)", "Estimate"] else NA_real_,
    se_quad        = if (!is.null(coefs)) coefs["I(soilTemp_avg_c^2)", "Std. Error"] else NA_real_,
    R2_marginal    = r2_marg,
    R2_conditional = r2_cond)}

# Apply the model to each KO within each pathway
ko_soilTemp_results <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  group_modify(~fit_KO_soilTemp(.x, .y)) %>%
  ungroup()
view(ko_soilTemp_results)

# Identify KOs with significant quadratic relationships
ko_soilTemp_sig <- ko_soilTemp_results %>%
  filter(p_quad < 0.05)
view(ko_soilTemp_sig)

# Plot

# Calculate mean KO abundance across samples within each pathway
ko_soilTemp_avg_mean_abundance_pathway <- biofilm_long_pathwayKO_env_top %>%
  group_by(biofilm_pathways_clean, biofilm_KOs) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop")

# Select the most abundant KO with a significant quadratic response per pathway
ko_soilTemp_avg_sig_abundant_pathway <- ko_soilTemp_results %>%
  filter(p_quad < 0.05) %>%
  left_join(
    ko_soilTemp_avg_mean_abundance_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs")) %>%
  group_by(biofilm_pathways_clean) %>%
  slice_max(mean_count, n = 1, with_ties = FALSE) %>%
  ungroup()

# Prepare regression statistics for figure annotations
ko_soilTemp_avg_sig_abundant_pathway <- ko_soilTemp_avg_sig_abundant_pathway %>%
  mutate(
    coef_linear = round(beta_linear, 2),
    coef_quad   = round(beta_quad, 2),
    p_linear_label = ifelse(
      p_linear < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_linear))),
    p_quad_label = ifelse(
      p_quad < 0.001,
      "P <0.001",
      paste0("P = ", scales::pvalue(p_quad))),
    r2_label = ifelse(R2_marginal < 0.001, "<0.001", round(R2_marginal, 3)),
    lm_label = glue::glue(
      "Linear: {coef_linear} ({p_linear_label})\n",
      "Quadratic: {coef_quad} ({p_quad_label})\n",
      "R²: {r2_label}"))

# Keep sample-level data for the selected KOs
biofilm_long_pathwayKO_env_top_soilTemp_selected <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_soilTemp_avg_sig_abundant_pathway,
    by = c("biofilm_pathways_clean","biofilm_KOs"))

# Arrange pathways and create facet labels
ko_soilTemp_avg_sig_abundant_pathway <- ko_soilTemp_avg_sig_abundant_pathway %>%
  arrange(biofilm_pathways_clean) %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_soilTemp_selected <- biofilm_long_pathwayKO_env_top_soilTemp_selected %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_soilTemp_avg_sig_abundant_pathway$pathway_KO_label)))


###### ...FigS16 @ #####
KO_soilTemp_avg_selected_plot <- ggplot(biofilm_long_pathwayKO_env_top_soilTemp_selected,
                                        aes(x = soilTemp_avg, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean), method = "lm", formula = y ~ x + I(x^2),
              se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE,
              color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50", size = 1.5, stroke = 0.35) +
  geom_label(data = ko_soilTemp_avg_sig_abundant_pathway,
             aes(x = Inf, y = -Inf, label = lm_label),
             inherit.aes = FALSE, hjust = 1, vjust = -4, size = 2.3,
             label.size = 0, fill = scales::alpha("white",0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(x = "Soil Temperature (°C)", y = "Total KO Abundance (counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text  = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

plot(KO_soilTemp_avg_selected_plot)

ggsave("Figure_S16.png", plot = KO_soilTemp_avg_selected_plot, units = "in", height = 10, width = 7, dpi = 1000)





###### 7.3.4. Selecting two representative KOs ######

# Select the two representative KOs with the highest marginal R² for each predictor (from the 8-pathway),
# retaining unique KO ids within each environmental predictor

# Soil moisture
ko_soilMoist_top2_from_pathways <- ko_soilMoist_avg_dayMin_sig_abundant_pathway %>%
  arrange(desc(R2_marginal)) %>%
  distinct(biofilm_KOs, .keep_all = TRUE) %>%
  slice_head(n = 2)

# Precipitation CV
ko_secPrecip_top2_from_pathways <- ko_secPrecip_cv_dayTot_sig_abundant_pathway %>%
  arrange(desc(R2_marginal)) %>%
  distinct(biofilm_KOs, .keep_all = TRUE) %>%
  slice_head(n = 2)

# Soil temperature
ko_soilTemp_top2_from_pathways <- ko_soilTemp_avg_sig_abundant_pathway %>%
  arrange(desc(R2_marginal)) %>%
  distinct(biofilm_KOs, .keep_all = TRUE) %>%
  slice_head(n = 2)

# Retain sample-level data for the selected KOs

# Soil moisture
biofilm_long_pathwayKO_env_top_soilMoist_top2 <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_soilMoist_top2_from_pathways,
    by = c("biofilm_pathways_clean", "biofilm_KOs"))

# Precipitation variability
biofilm_long_pathwayKO_env_top_secPrecip_top2 <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_secPrecip_top2_from_pathways,
    by = c("biofilm_pathways_clean", "biofilm_KOs"))

# Soil temperature
biofilm_long_pathwayKO_env_top_soilTemp_top2 <- biofilm_long_pathwayKO_env_top %>%
  semi_join(
    ko_soilTemp_top2_from_pathways,
    by = c("biofilm_pathways_clean", "biofilm_KOs"))

# Create facet labels for the selected KOs

# Soil moisture
ko_soilMoist_top2_from_pathways <- ko_soilMoist_top2_from_pathways %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_soilMoist_top2 <- biofilm_long_pathwayKO_env_top_soilMoist_top2 %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_soilMoist_top2_from_pathways$pathway_KO_label)))

# Precipitation variability
ko_secPrecip_top2_from_pathways <- ko_secPrecip_top2_from_pathways %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_secPrecip_top2 <- biofilm_long_pathwayKO_env_top_secPrecip_top2 %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_secPrecip_top2_from_pathways$pathway_KO_label)))

# Soil temperature
ko_soilTemp_top2_from_pathways <- ko_soilTemp_top2_from_pathways %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = paste0(biofilm_pathways_clean, "\n", biofilm_KOs)))

biofilm_long_pathwayKO_env_top_soilTemp_top2 <- biofilm_long_pathwayKO_env_top_soilTemp_top2 %>%
  mutate(
    pathway_KO_label = factor(
      paste0(biofilm_pathways_clean, "\n", biofilm_KOs),
      levels = levels(ko_soilTemp_top2_from_pathways$pathway_KO_label)))


# Generate regression plots for each environmental predictor (2 panels per variable)

# Soil moisture
p_ko_soilMoist_top2 <- ggplot(biofilm_long_pathwayKO_env_top_soilMoist_top2,
  aes(x = soilMoist_avg_dayMin, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean),
    method = "lm",formula = y ~ x + I(x^2),
    se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2),
    se = FALSE, color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50", size = 1.5, stroke = 0.35) +
  geom_label(data = ko_soilMoist_top2_from_pathways,
    aes(x = Inf, y = -Inf, label = lm_label),
    inherit.aes = FALSE, hjust = 1, vjust = -4, size = 2.3,
    label.size = 0, fill = scales::alpha("white", 0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(x = "Soil Moisture Day Min (%)", y = "Total KO Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

# Precipitation variability
p_ko_secPrecip_top2 <- ggplot(biofilm_long_pathwayKO_env_top_secPrecip_top2, 
                              aes(x = secPrecip_cv_dayTot, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean), method = "lm",
    formula = y ~ x + I(x^2), se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE,
    color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50", size = 1.5, stroke = 0.35) +
  geom_label(data = ko_secPrecip_top2_from_pathways,
    aes(x = Inf, y = -Inf, label = lm_label), inherit.aes = FALSE,
    hjust = 1, vjust = -4, size = 2.3, label.size = 0,
    fill = scales::alpha("white", 0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(x = "Precipitation Day Total CV (%)", y = "Total KO Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

# Soil temperature
p_ko_soilTemp_top2 <- ggplot(biofilm_long_pathwayKO_env_top_soilTemp_top2, aes(x = soilTemp_avg, y = count)) +
  geom_smooth(aes(fill = biofilm_pathways_clean),
    method = "lm", formula = y ~ x + I(x^2),
    se = TRUE, color = NA, alpha = 0.40) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = "black", linewidth = 0.7) +
  geom_point(shape = 21, color = "black", fill = "grey50",size = 1.5, stroke = 0.35) +
  geom_label(data = ko_soilTemp_top2_from_pathways, aes(x = Inf, y = -Inf, label = lm_label),
    inherit.aes = FALSE,hjust = 1, vjust = -4, size = 2.3, label.size = 0,
    fill = scales::alpha("white", 0), color = "black") +
  facet_wrap(~ pathway_KO_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pathway_colors) +
  labs(
    x = "Soil Temperature (°C)",
    y = "Total KO Abundance\n(counts / kb / million reads)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    strip.background = element_blank())

###### ...Fig4 @ #####
# Combine regression plots into Figure 4
KO_top2R2_from_pathway_selected_plot <- ggarrange(
  p_ko_soilMoist_top2, p_ko_secPrecip_top2, p_ko_soilTemp_top2,
  ncol = 1, nrow = 3,
  labels = c("A", "B", "C"),
  font.label = list(size = 11, face = "bold"))

plot(KO_top2R2_from_pathway_selected_plot)

ggsave("Figure_4.png", plot = KO_top2R2_from_pathway_selected_plot,
  units = "in", height = 8.4, width = 7, dpi = 1000)


