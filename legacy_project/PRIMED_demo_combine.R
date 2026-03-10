library(AnVIL)
library(readr)
#remotes::install_github("UW-GAC/prsmixsumstats", ref="rename")
library(prsmixsumstats)

adjusted <- FALSE
this_trait <- "TG"
validation <- "GeneSTAR"
drop_cohorts <- c("eMERGE-Geisinger", "eMERGE-KPWUW", "eMERGE-MountSinai")

if (adjusted) id_str <- "analysis_adjusted_id" else file_str <- "analysis_unadjusted_id"
if (adjusted) file_str <- "sumstats_adjusted" else file_str <- "sumstats_unadjusted"

sumst_tbl <- avtable(file_str) %>%
  filter(analysis == this_trait) %>%
  filter(biobank != validation) %>%
  filter(!(cohort %in% drop_cohorts))

drop_scores <- readLines("drop_scores_AoU_excl.txt")

dest_bucket <- paste0(avstorage(), "/pgs_", file_str, "/")
data_dir <- paste0("~/pgs_", file_str)
date_str <- format(Sys.Date(), "%Y%m%d")
clusters <- sort(unique(sumst_tbl$cluster))

tbl_list <- list()
for (clust in clusters) {
  print(clust)
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
    print(this[[paste0(file_str, "_id")]][i])
    sumst <- readRDS(file.path(data_dir, basename(this$sumstat_file[i])))
    covars <- colnames(sumst$xx)[!grepl("^PGS", colnames(sumst$xx))]
    print(sort(covars))
    # rename columns with incorrect names
    if ("bmi" %in% colnames(sumst$xx)) {
      sumst <- rename_col_sumstats(sumst, "bmi", "BMI")
    }
    overlap <- read_tsv(file.path(data_dir, basename(this$overlap_file[i])), show_col_types=FALSE)
    sumst <- filter_sumstats(sumst, overlap, name_col="score", 
                             filter_col="beta_fraction", 
                             filter_threshold=0.75)
    sumst <- filter_sumstats(sumst, overlap, name_col="score", 
                             filter_col="overlap", 
                             filter_threshold=0.75)
    sumst_list[[i]] <- sumst
  }
  sumst_comb <- combine_sumstats(sumst_list)
  
  # remove AoU and excluded scores
  sumst_comb$sumstats <- drop_cols_sumstats(sumst_comb$sumstats, drop_scores)
  sumst_comb$beta_multiplier <- sumst_comb$beta_multiplier[colnames(sumst_comb$sumstats$xx)]
  stopifnot(all(colnames(sumst_comb$sumstats$xx) == names(sumst_comb$beta_multiplier)))
  
  # check that we are not missing any covariates
  drop <- sumst_comb$incomplete_cols[!grepl("^PGS", sumst_comb$incomplete_cols)]
  if (length(drop) > 0) {
    print("drop:")
    print(drop)
    sumst_comb$sumstats <- drop_cols_sumstats(sumst_comb$sumstats, drop)
    sumst_comb$beta_multiplier <- sumst_comb$beta_multiplier[colnames(sumst_comb$sumstats$xx)]
    stopifnot(all(colnames(sumst_comb$sumstats$xx) == names(sumst_comb$beta_multiplier)))
  }
  
  covars <- colnames(sumst_comb$sumstats$xx)[!grepl("^PGS", colnames(sumst_comb$sumstats$xx))]
  print("covars:")
  print(sort(covars))
  
  # penalty factor should be 0 for covariates and 1 for PGS
  penalty_factor <- rep(0, ncol(sumst_comb$sumstats$xx))
  names(penalty_factor) <- colnames(sumst_comb$sumstats$xx)
  penalty_factor[grepl("^PGS", names(penalty_factor))] <- 1
  sumst_comb$penalty_factor <- penalty_factor
  
  print(paste("number of incomplete PRS:", length(sumst_comb$incomplete_cols)))
  print(paste("number of PRS with no variation:", length(sumst_comb$diag_zero_cols)))
  
  print(str(sumst_comb$sumstats))
  print(head(colnames(sumst_comb$sumstats$xx)))
  print(tail(colnames(sumst_comb$sumstats$xx)))
  outfile <- paste0(data_dir, "/", paste(this_trait, clust, date_str, file_str, sep="_"), ".rds")
  saveRDS(sumst_comb, outfile)
  avcopy(outfile, dest_bucket)
  tbl_list[[clust]] <- tibble(
    analysis = this_trait,
    cluster = clust,
    biobanks = paste(sort(unique(this$biobank)), collapse="|"),
    combined_sumstats = paste0(dest_bucket, basename(outfile))
  )
  tbl_list[[clust]][[id_str]] <- paste(this_trait, clust, date_str, sep="_")
}

tbl <- bind_rows(tbl_list)
avtable_import(tbl, entity=id_str)
