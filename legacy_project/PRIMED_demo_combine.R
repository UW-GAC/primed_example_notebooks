library(AnVIL)
library(readr)
remotes::install_github("UW-GAC/prsmixsumstats")
library(prsmixsumstats)

analysis <- avtable("example_analysis") %>%
  filter(trait == "WBC")

all <- analysis %>%
  filter(cluster == "all") %>%
  select(sumstat_file, overlap_file)
avcopy(all, "~/pgs_sumstats/")

sumst_list <- lapply(all$sumstat_file, function(x) readRDS(file.path("~/pgs_sumstats", basename(x))))
sumst_comb <- combine_sumstats(sumst_list)

length(sumst_comb$incomplete_cols)
str(sumst_comb$sumstats)


