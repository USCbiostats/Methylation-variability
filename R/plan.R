the_plan <-
  drake_plan(

    data_path = "~/Documents/USC/P4/JohnHopkins/John-Hopkins/data/",

    ## Load in datasets
    bulk_tumor_islands =   read.csv(path(data_path, "bulk-tumor_islands.csv")),
    bulk_tumor_genes =     read.csv(path(data_path, "bulk-tumor_genes.csv")),
    bulk_tumor_promoters = read.csv(path(data_path, "bulk-tumor_promoters.csv")),

    bulk_normal_islands =   read.csv(path(data_path, "bulk-normal_islands.csv")),
    bulk_normal_genes =     read.csv(path(data_path, "bulk-normal_genes.csv")),
    bulk_normal_promoters = read.csv(path(data_path, "bulk-normal_promoters.csv")),

    # Clean all datasets
    bulk_tumor_islands_region =  clean_datasets(bulk_tumor_islands),
    bulk_tumor_genes_region =    clean_datasets(bulk_tumor_genes),
    bulk_tumor_promoters_region = clean_datasets(bulk_tumor_promoters),

    bulk_normal_islands_region =  clean_datasets(bulk_normal_islands),
    bulk_normal_genes_region =    clean_datasets(bulk_normal_genes),
    bulk_normal_promoters_region = clean_datasets(bulk_normal_promoters),

    # Combine
    bulk_islands =   combine_datasets(bulk_tumor_islands_region,
                                      bulk_normal_islands_region),
    bulk_genes =     combine_datasets(bulk_tumor_genes_region,
                                      bulk_normal_genes_region),
    bulk_promoters = combine_datasets(bulk_tumor_promoters_region,
                                      bulk_normal_promoters_region),

    # Calculate variability
    bulk_islands_variability =   calc_variability(bulk_islands),
    bulk_genes_variability =     calc_variability(bulk_genes),
    bulk_promoters_variability = calc_variability(bulk_promoters),

    # Variability charts
    islands_variability_chart = variability_chart(bulk_islands_variability, "islands"),
    genes_variability_chart = variability_chart(bulk_genes_variability, "genes"),
    promoters_variability_chart = variability_chart(bulk_promoters_variability, "promoters"),

    report = target(
      command = {
        rmarkdown::render(knitr_in("doc/analysis.Rmd"))
        file_out("doc/analysis.html")
      }
    )
)
