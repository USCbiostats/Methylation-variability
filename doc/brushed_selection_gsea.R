
library(ggplot2)
library(dplyr)
library(shiny)
library(ReactomePA)

most_variable <- drake::readd(bulk_genes_variability) %>%
  dplyr::filter(n_obs >= drake::readd(min_obs)) %>%
  left_join(by = c("group" = "symbol"),
    data.frame(drake::readd(genes_annotations)) %>%
      dplyr::select(gene_id, symbol) %>%
      tidyr::drop_na() %>%
      distinct()
  ) %>%
   dplyr::select(gene_id, n_obs, tumor_variability, normal_variability)

ui <- basicPage(
  plotOutput("plot1", brush = "plot_brush", height = 400),
  plotOutput("info", height = 400)
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    ggplot(most_variable, aes(tumor_variability, normal_variability,
                              color = n_obs)) +
      geom_point(alpha = 0.1) +
      scale_color_viridis_c(trans = "log", labels = scales::comma) +
      theme_minimal()
  })

  output$info <- renderPlot({
    brushedPoints(most_variable, input$plot_brush) %>%
      pull(gene_id) %>%
      enrichPathway() %>%
      emapplot()
  }
  )

}

shinyApp(ui, server)

