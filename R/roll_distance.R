roll_distance <- function(data, basenames, distance) {
  gene_roll_annon <- build_annotations(genome = "hg19", annotations = "hg19_genes_cds")


  gene_regions <- gene_roll_annon[!is.na(gene_roll_annon$symbol)] %>%
    data.frame() %>%
    group_by(symbol, seqnames) %>%
    summarise(start = min(start),
              end = max(end)) %>%
    ungroup()

  roll_1 <- function(data, n) {
    data %>%
      mutate(start = start + n * distance,
             end = end + n * distance,
             symbol = paste0(symbol, "_", n))
  }

  roll_gene_regions <- makeGRangesFromDataFrame(map_dfr(0:5, ~ roll_1(gene_regions, .x)),
                                                keep.extra.columns = TRUE)

  roll_gene_annotated <- annotate_dataset(data, roll_gene_regions)

  get_metrics_by_group(roll_gene_annotated, paste0("X", basenames), annot.symbol) %>%
    separate(group, c("gene", "iter"), sep = "_")
}
