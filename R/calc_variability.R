#' Calculate variabilty within tumors and normals
#'
#' @param x tibble
#'
#' @return tibble
#' @export
calc_variability <- function(x) {
  x %>%
    mutate(tumor_variability = mhic_tumor + pwd_tumor * 2 +
                               var_sample_tumor * 4 + var_cpg_tumor * 4,
           normal_variability = mhic_normal + pwd_normal * 2 +
                                var_sample_normal * 4 + var_cpg_normal * 4)
}
