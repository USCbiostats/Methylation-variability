annotate_dataset <- function(data, annotations) {

  regions <- makeGRangesFromDataFrame(data,
                                      seqnames.field = "chr",
                                      start.field = "pos",
                                      end.field = "pos",
                                      keep.extra.columns = TRUE)

  dm_annotated <- annotate_regions(
    regions = regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  data.frame(dm_annotated)
}


slice_1_genes <- function(data) {
  data %>%
    group_by(seqnames, start) %>%
    dplyr::slice(1) %>%
    ungroup()
}
