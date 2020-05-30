scaler <- function(data) {
recipe(~ ., data = data) %>%
  step_normalize(matches("tumor|normal")) %>%
  prep()
}
