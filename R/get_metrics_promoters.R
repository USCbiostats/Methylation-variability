get_metrics_promoters <- function(basename, idat_folder, min_obs) {

  ## Loading the data
  target_df <- tibble(Basename = basename)

  betas <- read_idat(base = idat_folder,
                     targets = target_df,
                     force = TRUE)

  ## Preproccess

  promoter_ref <- get_promoter_ref()
  location_ref <- get_location_ref()

  beta_nested <- betas %>%
    left_join(location_ref, by = "rowname") %>%
    arrange(chr, pos) %>%
    select(-chr, -pos) %>%
    inner_join(promoter_ref, by = "rowname") %>%
    dplyr::filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
    select(rowname, Regulatory_Feature_Name, everything()) %>%
    select(-c(rowname, Regulatory_Feature_Group)) %>%
    nest(data = -c(Regulatory_Feature_Name))

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
