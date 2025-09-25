library(AnVIL)
library(readr)
remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

this_trait <- "WBC"
analysis <- avtable("example_analysis") %>%
  filter(analysis == this_trait)

aou_file <- "gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/AoU_pgs_id.csv"
#avcopy(aou_file, ".")
aou_scores <- stringr::str_trim(readLines(basename(aou_file)))

data_dir <- "~/pgs_sumstats"
clusters <- unique(analysis$cluster)
for (clust in clusters) {
  this <- analysis %>%
    filter(cluster == clust) %>%
    select(sumstat_file, overlap_file)
  for (f in unlist(this)) {
    local_file <- file.path(data_dir, basename(f))
    if (!file.exists(local_file)) suppressWarnings(avcopy(f, local_file))
  }
  
  sumst_list <- list()
  for (i in 1:nrow(this)) {
    sumst <- readRDS(file.path(data_dir, basename(this$sumstat_file[i])))
    overlap <- read_tsv(file.path(data_dir, basename(this$overlap_file[i])), show_col_types=FALSE)
    sumst_list[[i]] <- filter_sumstats(sumst, overlap, name_col="score", 
                             filter_col="beta_fraction", 
                             filter_threshold=0.8)
  }
  sumst_comb <- combine_sumstats(sumst_list)
  
  # remove AoU scores
  sumst_comb$sumstats <- drop_cols_sumstats(sumst_comb$sumstats, aou_scores)
  
  print(length(sumst_comb$incomplete_cols))
  print(str(sumst_comb$sumstats))
  saveRDS(sumst_comb, file.path(data_dir, paste(this_trait, clust, "sumstats.rds", sep="_")))
}


