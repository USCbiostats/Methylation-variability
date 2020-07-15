extract_most_variable <- function(data, n_most, min_obs) {
  data %>%
    dplyr::filter(n_obs >= min_obs) %>%
    mutate(variability = tumor_variability * normal_variability) %>%
    arrange(dplyr::desc(variability)) %>%
    dplyr::slice(seq_len(n_most)) %>%
    select(group, ends_with("variability"))
}
