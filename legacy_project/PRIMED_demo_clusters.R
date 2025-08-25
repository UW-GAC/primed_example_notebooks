library(readr)
source("primed_example_notebooks/legacy_project/make_sumstats_V2.R")
source("primed_example_notebooks/legacy_project/PRIMED_summary_stats.R")
source("primed_example_notebooks/legacy_project/PRIMED_stats_clusters.R")

cohort <- "ARIC"

# each of these files should have sample identifier as the first column
trait <- read_tsv(paste0(cohort, "_WBC_trait.tsv"))
covariates <- read_tsv(paste0(cohort, "_WBC_covariates.tsv"))
scores <- read_tsv(paste0(cohort, "_adjusted_scores_filtered.tsv"))
clusters <- read_tsv(paste0(cohort, "_cluster_definitions.tsv"))

sumst <- make_sumstats_clusters(trait, covariates, scores, clusters)
saveRDS(sumst, paste0(cohort, "_WBC_XX_XY.rds"))
avcopy(paste0(cohort, "_WBC_XX_XY.rds"), avstorage())

lapply(sumst, function(x) dim(x$xx))
lapply(sumst, function(x) head(x$xy))
lapply(sumst, function(x) attr(x$xy, "nsubj"))
lapply(sumst, function(x) attr(x$xy, "nmiss"))
lapply(sumst, function(x) attr(x$xy, "ysum"))
lapply(sumst, function(x) attr(x$xy, "yssq"))


