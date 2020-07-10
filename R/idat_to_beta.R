#' Reads idat files and calculate beta values with Noob
#'
#' idat_to_beta_values_noob
#'
#' @param base String, base directory. Folder where idat files are stored.
#' @param targets Data.frame with a column named 'Basename' in format
#' 'sentrix_id'_'terminus', for example 200511490070_R01C01.
#' @param force Should reading different size IDAT files be forced? Defaults to
#' TRUE.
#'
#' @return tibble
#' @export
read_idat <- function(base = NULL, targets = NULL, force = TRUE) {
  targets <- as.data.frame(targets)
  RGset <- read.metharray.exp(base = base, targets = targets, force = force)
  RGset@annotation = c(array = "IlluminaHumanMethylationEPIC",
                       annotation = "ilm10b4.hg19")

  MSet.noob <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE,
                              verbose = TRUE)

  ratioSet.noob <- ratioConvert(MSet.noob, what =  "both", keepCN = TRUE)
  beta.noob <- getBeta(ratioSet.noob)

  idat_to_tibble(beta.noob)
}

#' Turns idat matrix to tibble
#'
#' @param x idat matrix.
#'
#' @return tibble.
#' @export
idat_to_tibble <- function(x) {
  as.data.frame(x) %>%
    tibble::rownames_to_column()
}

new_read_idat <- function(base = NULL, targets = NULL, force = TRUE) {
  targets <- as.data.frame(targets)
  RGset <- read.metharray.exp(base = base, targets = targets, force = force)
  RGset@annotation = c(array = "IlluminaHumanMethylationEPIC",
                       annotation = "ilm10b4.hg19")

  MSet.noob <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE,
                              verbose = TRUE)

  ratioSet.noob <- ratioConvert(MSet.noob, what =  "both", keepCN = TRUE)
  beta.noob <- getBeta(ratioSet.noob)

  location_ref <- get_location_ref()

  location_ref %>%
    left_join(idat_to_tibble(beta.noob), by = "rowname") %>%
    arrange(chr, pos) %>%
    select(-rowname)
}
