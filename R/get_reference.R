get_island_ref <- function() {
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(Relation_to_Island == "Island")
}

get_location_ref <- function() {
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::select(-strand)
}

get_promoter_ref <- function() {
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::select(rowname, Regulatory_Feature_Name, Regulatory_Feature_Group)
}

get_gene_ref <- function() {
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    select(rowname, UCSC_RefGene_Name)
}

get_enhancer_ref <- function() {
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    select(rowname, Phantom4_Enhancers)
}
