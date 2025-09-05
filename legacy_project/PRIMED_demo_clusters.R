library(AnVIL)
library(readr)
remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

cohort <- "ARIC"
trait <- "WBC"
bucket <- "gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/summary_stats/"

# each of these files should have sample identifier as the first column
# all other columns should follow naming convention here:
# https://docs.google.com/spreadsheets/d/1B_AUoAlq9dY4GDxv2YgeE_tfDspqn2wJBD5N_P_fZYg/edit?usp=sharing
trait <- read_tsv(paste(cohort, trait, "trait.tsv", sep="_"))
covariates <- read_tsv(paste(cohort, trait, "covariates.tsv", sep="_"))
scores <- read_tsv(paste0(cohort, "_adjusted_scores.tsv"))
clusters <- read_tsv(paste0(cohort, "_cluster_definitions.tsv"))

make_sumstats_clusters(trait, covariates, scores, clusters, trait_name=trait, cohort_name=cohort)
avcopy(list.files(pattern="_sumstats.rds$"), bucket)

# check
sumst <- readRDS(paste(trait, cohort, "sumstats.rds", sep="_"))
str(sumst)


