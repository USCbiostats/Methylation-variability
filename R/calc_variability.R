#' Calculate variabilty within tumors and normals
#'
#' @param x tibble
#'
#' @return tibble
#' @export
calc_variability <- function(x) {
  x %>%
    mutate(tumor_variability = scale(mhic_tumor)[, 1] + scale(pwd_tumor)[, 1] +
                               scale(sample_var_tumor)[, 1] ,
           normal_variability = scale(mhic_normal)[, 1] + scale(pwd_normal)[, 1] +
                                scale(sample_var_normal)[, 1])
}
