## library() calls go here

suppressMessages({
library(conflicted)
library(dotenv)
library(drake)
library(fs)
library(dplyr)
library(tibble)
library(ggplot2)
library(glue)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(recipes)
library(knitr)
library(slider)

library(bedslider)

library(annotatr)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ReactomePA)
library(org.Hs.eg.db)
library(AnnotationDbi)
})

conflict_prefer("select", "dplyr")
