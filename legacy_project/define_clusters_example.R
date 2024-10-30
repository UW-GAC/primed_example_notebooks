library(AnVIL)
library(readr)
source("primed_example_notebooks/legacy_project/define_clusters.R")

gsutil_cp("gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/pca/kmeans_14clusters_centers.txt", ".")
cluster_centers <- read_tsv("kmeans_14clusters_centers.txt")

gsutil_cp("gs://fc-43902fa4-c6a5-42fa-acfc-0f3200f3162a/submissions/bef9de29-9f1a-4b7f-a21c-dfb0922518cf/projected_PCA/668cd970-d853-4287-8cd7-ff37a70f047f/call-ProjectArray/projections.txt", ".")
pcs <- read_tsv("projections.txt")

cluster_defs <- define_clusters(pcs, cluster_centers)

head(cluster_defs)
table(cluster_defs$best_cluster)
