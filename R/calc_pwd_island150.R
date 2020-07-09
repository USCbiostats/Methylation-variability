calc_pwd_island150 <- function(data) {
  pwd_df <- function(data) {

    if (nrow(data) == 1) return(NaN)
    as.matrix(data) %>%
      apply(2, function(x) dist(x, method = "manhattan") %>% mean()) %>%
      mean()
  }

  ref <- data %>%
    ungroup() %>%
    mutate(row_number = row_number()) %>%
    dplyr::select(row_number, everything()) %>%
    group_by(chr)

  pwds <- ref %>%
    mutate(groups = slide_index_chr(row_number, .i = pos, ~glue_collapse(.x, sep = " "),
                                    .before = 75, .after = 75)) %>%
    separate_rows(groups, sep = " ", convert = TRUE) %>%
    select(-n_obs, -Relation_to_Island, -n_obs_in_island, -rowname, -row_number, -Islands_Name, -pos) %>%
    nest(matches("[0-9]{12}")) %>%
    mutate(pwd = map_dbl(data, pwd_df)) %>%
    select(groups, pwd)

  ref %>%
    select(chr, pos, row_number) %>%
    left_join(pwds, by = c("row_number" = "groups")) %>%
    select(-row_number)
}
