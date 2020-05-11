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
    guides(color = "none") +
    annotate("segment", x = 0.5, xend = 1.5, y = -0.1, yend = -0.1, arrow = arrow()) +
    annotate("text", x = 1, y = -0.05, label = "Most variable") +
    annotate("segment", y = 0.5, yend = 1.5, x = -0.1, xend = -0.1, arrow = arrow()) +
    annotate("text", y = 1, x = -0.05, label = "Most variable", angle = 90)
}
