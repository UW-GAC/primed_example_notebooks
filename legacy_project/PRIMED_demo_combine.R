library(AnVIL)
library(readr)
#remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

#this_trait <- "WBC"
this_trait <- "CAD_subcovar"
sumst_tbl <- avtable("sumstats") %>%
  filter(analysis == this_trait)

keep_scores <- readLines("harmonized_scores_noAoU_noexcl.txt")

dest_bucket <- file.path(avstorage(), "pgs_sumstats/")
data_dir <- "~/pgs_sumstats"
date_str <- format(Sys.Date(), "%Y%m%d")
clusters <- unique(sumst_tbl$cluster)
#clusters <- "all"

tbl_list <- list()
for (clust in clusters) {
  this <- sumst_tbl %>%
    filter(cluster == clust) 
  files <- this %>%
    select(sumstat_file, overlap_file)
  for (f in unique(unlist(files))) {
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
  
  # remove AoU and excluded scores
  cols_to_drop <- setdiff(colnames(sumst_comb$sumstats$xx), keep_scores)
  sumst_comb$sumstats <- drop_cols_sumstats(sumst_comb$sumstats, cols_to_drop)
  
  # check that we are not missing any covariates
  stopifnot(all(str_detect(sumst_comb$incomplete_cols, "^PGS")))
  print(paste("number of incomplete PRS:", length(sumst_comb$incomplete_cols)))
  
  print(str(sumst_comb$sumstats))
  outfile <- file.path(data_dir, paste(this_trait, clust, date_str, "sumstats.rds", sep="_"))
  saveRDS(sumst_comb, outfile)
  avcopy(outfile, dest_bucket)
  tbl_list[[clust]] <- tibble(
    analysis = this_trait,
    cluster = clust,
    biobanks = paste(sort(unique(this$biobank)), collapse="|"),
    combined_sumstats = paste0(dest_bucket, basename(outfile)),
    analysis_id = paste(this_trait, cluster, date_str, sep="_")
  )
}

tbl <- bind_rows(tbl_list)
avtable_import(tbl, entity="analysis_id")
