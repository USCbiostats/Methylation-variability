get_metrics_cpgs <- function(cpg_vec, basename, idat_folder) {

  target_df <- tibble(Basename = basename)

  betas <- read_idat(base = idat_folder,
                     targets = target_df,
                     force = TRUE)

  location_ref <- get_location_ref()

  betas %>%
    right_join(tibble(rowname = cpg_vec), by = "rowname") %>%
    left_join(location_ref, by = "rowname") %>%
    arrange(chr, pos) %>%
    select(chr, pos, everything()) %>%
    select(-rowname)
}
