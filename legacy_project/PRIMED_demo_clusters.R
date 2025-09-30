remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)
library(readr)

cohort <- "ARIC"
analysis_name <- "WBC"

# each of these files should have sample identifier as the first column
# all other columns should follow naming convention here:
# https://docs.google.com/spreadsheets/d/1B_AUoAlq9dY4GDxv2YgeE_tfDspqn2wJBD5N_P_fZYg/edit?usp=sharing
# example files are here:
# https://github.com/UW-GAC/prsmixsumstats/tree/main/inst/extdata
trait <- read_tsv(paste(cohort, analysis_name, "trait.tsv", sep="_"))
covariates <- read_tsv(paste(cohort, analysis_name, "covariates.tsv", sep="_"))
scores <- read_tsv(paste0(cohort, "_adjusted_scores.tsv"))
clusters <- read_tsv(paste0(cohort, "_cluster_definitions.tsv"))

# create summary statistics
make_sumstats_clusters(trait = trait, covariates = covariates, 
                       scores = scores, clusters = clusters, 
                       analysis_name = analysis_name, cohort_name = cohort)

# check
sumst <- readRDS(paste(analysis_name, cohort, "sumstats.rds", sep="_"))
str(sumst)
# repeat for cluster files

# copy files to google bucket
library(AnVIL)
library(AnVILGCP)
bucket <- "gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/summary_stats/"
avcopy(list.files(pattern="_sumstats.rds$"), bucket)
# copy overlap file
avcopy("gs://fc-secure-d5a0ec8b-663b-4a69-8e57-1f62ee70e694/ARIC_overlap.tsv", bucket)
