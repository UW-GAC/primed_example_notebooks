library(readr)
source("primed_example_notebooks/legacy_project/PRIMED_summary_stats.R")
source("primed_example_notebooks/legacy_project/PRIMED_stats_clusters.R")

# each of these files should have sample identifier as the first column
trait <- read_tsv("JHS_WBC_trait.tsv")
covariates <- read_tsv("JHS_WBC_covariates.tsv")
scores <- read_tsv("JHS_adjusted_scores_filtered.tsv")
clusters <- read_tsv("JHS_cluster_definitions.tsv")

sumst <- make_sumstats_clusters(trait, covariates, scores, clusters)
saveRDS(sumst, "JHS_WBC_XX_XY.rds")
avcopy("JHS_WBC_XX_XY.rds", avstorage())

lapply(sumst, function(x) dim(x$xx))
lapply(sumst, function(x) head(x$xy))
lapply(sumst, function(x) attr(x$xy, "nsubj"))
lapply(sumst, function(x) attr(x$xy, "nmiss"))
lapply(sumst, function(x) attr(x$xy, "ysum"))
lapply(sumst, function(x) attr(x$xy, "yssq"))


