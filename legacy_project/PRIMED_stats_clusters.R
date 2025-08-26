# make summary stats for each cluster

# each input is expected to be a data table where the first column is sample ID
make_sumstats_clusters <- function(trait, covariates, scores, clusters, analysis_name, cohort_name) {
  ids <- intersect(trait[[1]], covariates[[1]])
  ids <- intersect(ids, scores[[1]])
  ids <- intersect(ids, clusters[[1]])
  trait <- trait[match(ids, trait[[1]]),]
  covariates <- covariates[match(ids, covariates[[1]]),]
  scores <- scores[match(ids, scores[[1]]),]
  clusters <- clusters[match(ids, clusters[[1]]),]
  stopifnot(all(trait[[1]] == covariates[[1]]))
  stopifnot(all(trait[[1]] == scores[[1]]))
  stopifnot(all(trait[[1]] == clusters[[1]]))
  
  cov_scores <- cbind(covariates, scores[,-1])
  rm(scores)
  rm(covariates)
  
  res <- make_sumstats_V2(x=as.matrix(cov_scores[,-1]), y=unlist(trait[,-1]))
  saveRDS(res, paste0(analysis_name, "_", cohort_name, "_sumstats.rds"))
  
  cluster_names <- unique(clusters[[2]])
  for (c in cluster_names) {
    index <- which(clusters[[2]] %in% c)
    res <- make_sumstats_V2(x=as.matrix(cov_scores[index,-1]), y=unlist(trait[index,-1]))
    saveRDS(res, paste0(analysis_name, "_", cohort_name, "_cluster", c, "_sumstats.rds"))
  }
}
