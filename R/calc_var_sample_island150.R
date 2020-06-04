calc_var_sample_island150 <- function(data) {

  ref <- data %>%
    ungroup() %>%
    mutate(row_number = row_number()) %>%
    dplyr::select(row_number, everything()) %>%
    group_by(chr)

  pwds <- ref %>%
    mutate(groups = slide_index_chr(row_number, .i = pos, ~paste(.x, sep = "-", collapse = "-"),
                                    .before = 75, .after = 75)) %>%
    ungroup() %>%
    separate_rows(groups, sep = "-", convert = TRUE) %>%
    select(-n_obs, -Relation_to_Island, -n_obs_in_island, -rowname, -row_number, -Islands_Name, -pos) %>%
    nest(matches("[0-9]{12}")) %>%
    mutate(var_sample = map_dbl(data, var_df)) %>%
    select(groups, var_sample)

  ref %>%
    select(chr, pos, row_number) %>%
    left_join(pwds, by = c("row_number" = "groups")) %>%
    select(-row_number)
}
