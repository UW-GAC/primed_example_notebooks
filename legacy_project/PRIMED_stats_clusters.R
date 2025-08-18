# make summary stats for each cluster

# each input is expected to be a data table where the first column is sample ID
make_sumstats_clusters <- function(trait, covariates, scores, clusters) {
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
  
  res <- list()
  res[["all"]] <- make_sumstats(x=cov_scores[,-1], y=trait[,-1])
  
  cluster_names <- unique(clusters[[2]])
  for (c in cluster_names) {
    index <- which(clusters[[2]] %in% c)
    res[[paste0("cluster", c)]] <- make_sumstats(x=cov_scores[index,-1], y=trait[index,-1])
  }
  
  return(res)
}
