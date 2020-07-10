variability_chart <- function(data, type, min_obs) {
  data %>%
    dplyr::filter(n_obs >= min_obs) %>%
    ggplot(aes(tumor_variability, normal_variability,
               color = n_obs,
               text = group)) +
    geom_point(alpha = 0.3) +
    scale_color_viridis_c(trans = "log", labels = scales::comma) +
    theme_minimal() +
    geom_abline(slope = 1, intercept = 0) +
    coord_fixed() +
    labs(title = glue("Variability for {type}"),
         x = "tumor variability",
         y = "normal variability",
         color = "# Cpgs")
}
