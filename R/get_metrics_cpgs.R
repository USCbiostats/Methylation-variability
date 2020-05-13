get_metrics_cpgs <- function(cpg_vec, basename, idat_folder, min_obs) {

target_df <- tibble(Basename = basename)

betas <- read_idat(base = idat_folder,
                   targets = target_df,
                   force = TRUE) %>%
  idat_to_tibble()

location_ref <- get_location_ref()

beta_nested <- betas %>%
  right_join(tibble(rowname = cpg_vec), by = "rowname") %>%
  left_join(location_ref, by = "rowname") %>%
  arrange(chr, pos) %>%
  select(-chr, -pos) %>%
  mutate(group = ceiling(row_number() / 20)) %>%
  select(rowname, group, everything()) %>%
  select(-c(rowname)) %>%
  nest(data = -c(group)) %>%
  ungroup()

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
