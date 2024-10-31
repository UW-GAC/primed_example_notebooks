define_clusters <- function(pcs, cluster_centers) {
  n_pcs <- ncol(cluster_centers)
  sample_ids <- pcs[,2]
  pcs <- pcs[,3:ncol(pcs)]
  if (ncol(pcs) < n_pcs) stop(paste("only", ncol(pcs), "PCs in pcs but need", n_pcs, "to define clusters"))
  
  n_samples <- nrow(pcs)
  n_clusters <- nrow(cluster_centers)
  cluster_distances <- matrix(nrow=n_samples, ncol=n_clusters, 
                              dimnames=list(NULL, paste0("cluster", 1:n_clusters, "_dist")))
  for (n in 1:n_samples) {
    for (c in 1:n_clusters) {
      cluster_distances[n,c] <- dist(rbind(cluster_centers[c,], pcs[n,]), method = "euclidean")
    }
  }
  best_cluster <- apply(cluster_distances, 1, which.min)
  cluster_defs <- cbind(sample_ids, best_cluster, cluster_distances)
  return(cluster_defs)
}
