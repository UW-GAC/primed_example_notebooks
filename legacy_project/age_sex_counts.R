library(dplyr)
library(readr)

age_cat_fn <- function(age) {
  cut(age,
      breaks = c(-Inf, 25, seq(30, 75, by = 5), Inf),
      labels = c(
        "<=25",
        ">25 to <=30",
        ">30 to <=35",
        ">35 to <=40",
        ">40 to <=45",
        ">45 to <=50",
        ">50 to <=55",
        ">55 to <=60",
        ">60 to <=65",
        ">65 to <=70",
        ">70 to <=75",
        ">75"
      ),
      right = TRUE,      # intervals are (a, b]
      include.lowest = TRUE,
      ordered_result = TRUE
  )
}


# binary outcome
trait <- read_tsv("ARIC_T2D_trait.tsv")
covs <- read_tsv("ARIC_T2D_covariates.tsv")

tbl <- trait %>%
  left_join(covs) %>%
  mutate(age_cat = age_cat_fn(age)) %>%
  mutate(T2D=as.factor(T2D), sex=as.factor(sex)) %>%
  count(age_cat, T2D, sex, .drop=FALSE)

write_tsv(tbl, "ARIC_T2D_age_sex_counts.tsv")


# quantitative outcome
trait <- read_tsv("ARIC_WBC_trait.tsv")
covs <- read_tsv("ARIC_WBC_covariates.tsv")

tbl <- trait %>%
  left_join(covs) %>%
  mutate(age_cat = age_cat_fn(age)) %>%
  mutate(sex=as.factor(sex)) %>%
  count(age_cat, sex, .drop=FALSE)

write_tsv(tbl, "ARIC_WBC_age_sex_counts.tsv")
