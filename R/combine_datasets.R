#' Combine tumor and normal datasets together
#'
#' @param tumor Tibble with tumor regions, output from `clean_datasets()`
#' @param normal Tibble with normal regions, output from `clean_datasets()`
#'
#' @return Tibble
#' @export
combine_datasets <- function(tumor, normal) {
  left_join(
    by = c("region", "n_obs"),
    suffix = c("_tumor", "_normal"),
    tumor,
    normal
  )
}
