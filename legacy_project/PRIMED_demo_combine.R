library(AnVIL)
library(readr)
remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

this_trait <- "WBC"
analysis <- avtable("example_analysis") %>%
  filter(trait == this_trait)

clusters <- unique(analysis$cluster)
for (clust in clusters) {
  this <- analysis %>%
    filter(cluster == clust) %>%
    select(sumstat_file, overlap_file)
  avcopy(unlist(this), "~/pgs_sumstats/")
  
  sumst_list <- list()
  for (i in 1:nrow(this)) {
    sumst <- readRDS(file.path("~/pgs_sumstats", basename(this$sumstat_file[i])))
    overlap <- read_tsv(file.path("~/pgs_sumstats", basename(this$overlap_file[i])))
    sumst_list[[i]] <- filter_sumstats(sumst, overlap, name_col="score", 
                             filter_col="beta_fraction", 
                             filter_threshold=0.8)
  }
  sumst_comb <- combine_sumstats(sumst_list)
  
  #length(sumst_comb$incomplete_cols)
  #str(sumst_comb$sumstats)
  saveRDS(sumst_comb, paste(this_trait, clust, "sumstats.rds", sep="_"))
}


