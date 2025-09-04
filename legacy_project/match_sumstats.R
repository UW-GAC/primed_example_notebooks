match_sumstats <- function(sumstats) {
  cols <- lapply(sumstats, function(x) names(x$xy))
  cols_intersect <- cols[[1]]
  cols_union <- cols[[1]]
  for (i in 2:length(cols)) {
    cols_intersect <- intersect(cols_intersect, cols[[i]])
    cols_union <- union(cols_union, cols[[i]])
  }
  writeLines(setdiff(cols_union, cols_intersect), "incomplete_scores.txt")
  tmp <- list()
  for (i in 1:length(sumstats)) {
    missing_cols <- setdiff(cols_union, cols[[i]])
    tmp[[i]] <- add_cols_sumstats(sumstats[[i]], missing_cols)
    tmp[[i]] <- match_cols_sumstats(tmp[[i]], cols_union)
  }
  return(tmp)
}


match_cols_sumstats <- function(ss, col_names) {
  tmp <- list()
  tmp$xx <- ss$xx[col_names, col_names]
  tmp$xy <- ss$xy[col_names]
  return(tmp)
}


add_cols_sumstats <- function(ss, col_names) {
  tmp <- list()
  tmp$xx <- add_cols_square_matrix(ss$xx, col_names)
  tmp$xy <- c(ss$xy, setNames(rep(0, length(col_names)), col_names))
  return(tmp)
}


add_cols_square_matrix <- function(x, col_names) {
  xc <- matrix(0, ncol=length(col_names), nrow=nrow(x), dimnames=list(NULL, col_names))
  x <- cbind(x, xc)
  xr <- matrix(0, ncol=ncol(x), nrow=length(col_names), dimnames=list(col_names, NULL))
  x <- rbind(x, xr)
  return(x)
}

# testing
x <- matrix(1, nrow=3, ncol=3, dimnames=list(letters[1:3], letters[1:3]))
ss <- list(xx=x, xy=setNames(1:3, letters[1:3]))
col_names <- c("d", "e")
add_cols_square_matrix(x, col_names)
add_cols_sumstats(ss, col_names)
match_cols_sumstats(ss, c("c", "b", "a"))
