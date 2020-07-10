get_island150_data <- function(basename, idat_folder) {
  ## Loading the data
  target_df <- tibble(Basename = basename)

  betas <- read_idat(base = idat_folder,
                     targets = target_df,
                     force = TRUE)

  ## Preproccess

  island_ref <- get_island_ref()
  location_ref <- get_location_ref()

  sorted_betas <- betas %>%
    left_join(location_ref, by = "rowname") %>%
    arrange(chr, pos) %>%
    full_join(island_ref, by = "rowname")

  sorted_betas %>%
    dplyr::filter(!is.na(pos)) %>%
    group_by(chr) %>%
    mutate(n_obs = slide_index_dbl(pos, pos, ~ length(.x), .before = 75, .after = 75),
           n_obs_in_island = slide_index_dbl(Relation_to_Island, .i =  pos,
                                             ~ sum(.x == "Island", na.rm = TRUE),
                                             .before = 75, .after = 75)) %>%
    dplyr::filter(n_obs_in_island > 0) %>%
    select(chr, pos, Relation_to_Island, n_obs, n_obs_in_island, everything()) %>%
    ungroup()
}
