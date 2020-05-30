scaler <- function(data) {
  recipe(~ ., data = data) %>%
    step_normalize(matches("tumor|normal")) %>%
    prep(strings_as_factors = FALSE)
}

apply_scaler <- function(scaler_recipe, data) {
  bake(scaler_recipe, data)
}
