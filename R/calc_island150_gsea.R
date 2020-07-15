calc_island150_gsea <- function(data, order,
                                island150_min_obs, island150_n_genes,
                                genes_annotations, island_annotations) {

data %>%
  dplyr::filter(n_obs >= island150_min_obs) %>%
  arrange( {{order}} ) %>%
  dplyr::slice(seq_len(island150_n_genes)) %>%
  annotate_nearest_gene_to_island(genes_annotations, island_annotations) %>%
  pull(gene_id) %>%
  enrichPathway()
}
