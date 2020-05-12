#' Clean summerised results
#'
#' - Renames first column of data.frame to "region"
#' - Returns a tibble
#' - Removes acf and n_obs metrics
#' @param x A data.frame
#'
#' @return A tibble
#' @export
clean_datasets <- function(x) {
  x %>%
    dplyr::rename(region = 1) %>%
    select(-n_obs)
}

