get_metrics_enhancers <- function(basename, idat_folder, min_obs) {

  ## Loading the data
  target_df <- tibble(Basename = basename)

  betas <- read_idat(base = idat_folder,
                     targets = target_df,
                     force = TRUE) %>%
    idat_to_tibble()

  ## Preproccess

  enhancer_ref <- get_enhancer_ref()
  location_ref <- get_location_ref()

  beta_nested <- betas %>%
    left_join(location_ref, by = "rowname") %>%
    arrange(chr, pos) %>%
    select(-chr, -pos) %>%
    inner_join(enhancer_ref, by = "rowname") %>%
    dplyr::filter(Phantom4_Enhancers != "") %>%
    select(rowname, Phantom4_Enhancers, everything()) %>%
    select(-c(rowname)) %>%
    nest(data = -c(Phantom4_Enhancers))

  ## Calculate matrics

  beta_metrics <- beta_nested %>%
    mutate(n_obs = map_dbl(data, nrow)) %>%
    dplyr::filter(n_obs >= min_obs) %>%
    mutate(mhic = map_dbl(data, mhic_vec),
           pwd = map_dbl(data, pwd_df),
           var_sample = map_dbl(data, var_df),
           var_cpg = map_dbl(data, cpg_var_df)) %>%
    select(-data)
}
