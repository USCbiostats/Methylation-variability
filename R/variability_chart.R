variability_chart <- function(data, type) {
  data %>%
    select(region, tumor_variability, normal_variability) %>%
    ggplot(aes(tumor_variability, normal_variability,
               color = tumor_variability * normal_variability,
               text = region)) +
    geom_point(alpha = 0.3) +
    scale_color_viridis_c() +
    theme_minimal() +
    geom_abline(slope = 1, intercept = 0) +
    coord_fixed() +
    labs(title = glue("Variability for {type}"),
         x = "tumor variability",
         y = "normal variability") +
    guides(color = "none")
}
