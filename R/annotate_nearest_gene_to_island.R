annotate_nearest_gene_to_island <- function(data, genes_annotations, island_annotations) {
  genes_annotations <- genes_annotations[!is.na(genes_annotations$symbol)]

  island_gene <- nearest(island_annotations, genes_annotations)

  ccc <- data %>%
    mutate(nearest_gene = genes_annotations$symbol[island_gene[readr::parse_number(annot.id)]])

  left_join(by = c("nearest_gene" = "symbol"),
            ccc,
            data.frame(genes_annotations) %>%
              select(gene_id, symbol) %>%
              distinct()
  )
}
