#' This function takes a data.frame. then calculates the procentage of values
#' between 0.2 and 0.8
mhic_vec <- function(col) {
  unlist(col) %>%
    between(0.2, 0.8) %>%
    mean(rm.na = TRUE)
}

#' This function takes a data.frame. then performs acf on the columns one by
#' one where after it is taking its mean.
acf_df_lag <- function(data, lag) {
  mean(map_dbl(data, ~ acf(.x, lag.max = lag, plot = FALSE,
                           na.action = na.pass)$acf[lag + 1]))
}

#' This function takes a data.frame. then it calculates pairwise distance using
#' manhattan distances and then takes the mean
pwd_df <- function(data) {
  as.matrix(data) %>%
    apply(2, function(x) dist(x, method = "manhattan") %>% mean()) %>%
    mean()
}

#' This function takes a data.frame. then calculates the variance within each
#' column and takes the average.
var_df <- function(data) {
  map_dbl(data, var) %>%
    mean()
}

#' This function calculates the variance on a cpg level basis. Then takes the mean.
cpg_var_df <- function(data) {
  apply(X = data, MARGIN = 1, FUN = var) %>%
    mean()
}
