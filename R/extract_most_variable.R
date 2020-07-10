extract_most_variable <- function(data, n_most) {
  data %>%
    mutate(variability = tumor_variability * normal_variability) %>%
    arrange(dplyr::desc(variability)) %>%
    dplyr::slice(seq_len(n_most)) %>%
    select(group, ends_with("variability"))
}
