calc_mhic_island150 <- function(data) {
  mhic <- function(x) {
    between(x, 0.2, 0.8) %>% mean(rm.na = TRUE)
  }

  mhic_slide <- function(x, pos) {
    slide_index_dbl(x, pos, .f = mhic, .before = 75, .after = 75)
  }

  data %>%
    mutate_at(vars(matches("[0-9]{12}")), ~mhic_slide(.x, pos)) %>%
    ungroup() %>%
    mutate(mhic_combined = select(., matches("[0-9]{12}")) %>% rowMeans(na.rm = TRUE))
}
