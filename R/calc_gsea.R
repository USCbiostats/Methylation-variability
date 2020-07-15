calc_gsea <- function(data, order, n_genes, min_obs, genes_annotations) {
  most_variable <- data %>%
    dplyr::filter(n_obs >= min_obs) %>%
    arrange( {{order}} ) %>%
    dplyr::slice(seq_len(n_genes)) %>%
    pull(group)

  most_vairable_gene_id <- data.frame(genes_annotations) %>%
    dplyr::select(gene_id, symbol) %>%
    tidyr::drop_na() %>%
    distinct() %>%
    dplyr::filter(symbol %in% most_variable) %>%
    pull(gene_id)

  enrichPathway(gene = most_vairable_gene_id,
                pvalueCutoff = 0.05,
                readable = TRUE)
}
