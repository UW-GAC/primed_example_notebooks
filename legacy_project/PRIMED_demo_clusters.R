library(AnVIL)
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

make_sumstats_clusters(trait, covariates, scores, clusters, analysis_name="WBC", cohort_name=cohort)
avcopy(list.files(pattern="_sumstats.rds$"), avstorage())

sumst <- readRDS("WBC_ARIC_sumstats.rds")
dim(sumst$xx)
head(sumst$xy)
attr(sumst$xy, "nsubj")
attr(sumst$xy, "nmiss")
attr(sumst$xy, "ysum")
attr(sumst$xy, "yssq")


