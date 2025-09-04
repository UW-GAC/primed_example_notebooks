library(AnVIL)
library(readr)
remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

cohort <- "ARIC"

# each of these files should have sample identifier as the first column
trait <- read_tsv(paste0(cohort, "_WBC_trait.tsv"))
covariates <- read_tsv(paste0(cohort, "_WBC_covariates.tsv"))
scores <- read_tsv(paste0(cohort, "_adjusted_scores_filtered.tsv"))
clusters <- read_tsv(paste0(cohort, "_cluster_definitions.tsv"))

make_sumstats_clusters(trait, covariates, scores, clusters, trait_name="WBC", cohort_name=cohort)
avcopy(list.files(pattern="_sumstats.rds$"), avstorage())
sumst <- readRDS("WBC_ARIC_sumstats.rds")
dim(sumst$xx)
head(sumst$xy)
attr(sumst, "nsubj")
attr(sumst, "nmiss")
attr(sumst, "ysum")
attr(sumst, "yssq")


