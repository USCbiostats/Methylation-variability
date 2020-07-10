the_plan <-
  drake_plan(

    min_obs = 10,
    idat_folder = "~/Data/methylation/idat",
    data_directory = "~/Data/methylation/directory.csv",

    n_most = 500,
    data_path = "~/Documents/USC/P4/JohnHopkins/John-Hopkins/data/",

    # extract basenames from data_directory
    bulk_normal_basenames = read_csv(data_directory) %>%
      dplyr::filter(jh_11_19) %>%
      dplyr::filter(str_detect(sample_id, ".{1}N")) %>%
      pull(basename),

    bulk_tumor_basenames = read_csv(data_directory) %>%
      dplyr::filter(jh_11_19) %>%
      dplyr::filter(!str_detect(sample_id, ".{1}N")) %>%
      pull(basename),

    ## Read Idat files
    bulk_tumor_idat = new_read_idat(idat_folder, tibble(Basename = bulk_tumor_basenames)),
    bulk_normal_idat = new_read_idat(idat_folder, tibble(Basename = bulk_normal_basenames)),

    ## Annotations
    island_annotations    = build_annotations(genome = 'hg19', annotations = "hg19_cpg_islands"),
    genes_annotations     = build_annotations(genome = "hg19", annotations = "hg19_basicgenes"),
    promoters_annotations = build_annotations(genome = "hg19", annotations = "hg19_genes_promoters"),
    enhancers_annotations = build_annotations(genome = "hg19", annotations = "hg19_enhancers_fantom"),

    ## Annotated datasets
    bulk_tumor_islands   = annotate_dataset(bulk_tumor_idat, island_annotations),
    bulk_tumor_genes     = annotate_dataset(bulk_tumor_idat, genes_annotations) %>%
      slice_1_genes() %>%
      dplyr::filter(!is.na(annot.symbol)),
    bulk_tumor_promoters = annotate_dataset(bulk_tumor_idat, promoters_annotations) %>%
      dplyr::filter(!is.na(annot.symbol)),
    bulk_tumor_enhancers = annotate_dataset(bulk_tumor_idat, enhancers_annotations),

    tumor_list_a = get_metrics_cpgs(readr::read_rds("data/A.rds"), bulk_tumor_basenames, idat_folder) %>%
      mutate(group = ceiling(row_number() / 10)),
    tumor_list_b = get_metrics_cpgs(readr::read_rds("data/B.rds"), bulk_tumor_basenames, idat_folder) %>%
      mutate(group = ceiling(row_number() / 10)),


    bulk_normal_islands   = annotate_dataset(bulk_normal_idat, island_annotations),
    bulk_normal_genes     = annotate_dataset(bulk_normal_idat, genes_annotations) %>%
      slice_1_genes() %>%
      dplyr::filter(!is.na(annot.symbol)),
    bulk_normal_promoters = annotate_dataset(bulk_normal_idat, promoters_annotations) %>%
      dplyr::filter(!is.na(annot.symbol)),
    bulk_normal_enhancers = annotate_dataset(bulk_normal_idat, enhancers_annotations),

    normal_list_a = get_metrics_cpgs(readr::read_rds("data/A.rds"), bulk_normal_basenames, idat_folder) %>%
      mutate(group = ceiling(row_number() / 10)),
    normal_list_b = get_metrics_cpgs(readr::read_rds("data/B.rds"), bulk_normal_basenames, idat_folder) %>%
      mutate(group = ceiling(row_number() / 10)),

    ## Calculate metrics
    bulk_tumor_islands_metrics   = get_metrics_by_group(bulk_tumor_islands,   paste0("X", bulk_tumor_basenames), annot.id),
    bulk_tumor_genes_metrics     = get_metrics_by_group(bulk_tumor_genes,     paste0("X", bulk_tumor_basenames), annot.symbol),
    bulk_tumor_promoters_metrics = get_metrics_by_group(bulk_tumor_promoters, paste0("X", bulk_tumor_basenames), annot.symbol),
    bulk_tumor_enhancers_metrics = get_metrics_by_group(bulk_tumor_enhancers, paste0("X", bulk_tumor_basenames), annot.id),
    tumor_list_a_metrics         = get_metrics_by_group(tumor_list_a, bulk_tumor_basenames, group),
    tumor_list_b_metrics         = get_metrics_by_group(tumor_list_b, bulk_tumor_basenames, group),

    bulk_normal_islands_metrics   = get_metrics_by_group(bulk_normal_islands,   paste0("X", bulk_normal_basenames), annot.id),
    bulk_normal_genes_metrics     = get_metrics_by_group(bulk_normal_genes,     paste0("X", bulk_normal_basenames), annot.symbol),
    bulk_normal_promoters_metrics = get_metrics_by_group(bulk_normal_promoters, paste0("X", bulk_normal_basenames), annot.symbol),
    bulk_normal_enhancers_metrics = get_metrics_by_group(bulk_normal_enhancers, paste0("X", bulk_normal_basenames), annot.id),
    normal_list_a_metrics         = get_metrics_by_group(normal_list_a, bulk_normal_basenames, group),
    normal_list_b_metrics         = get_metrics_by_group(normal_list_b, bulk_normal_basenames, group),

    ## Combine
    bulk_islands =   combine_datasets(bulk_tumor_islands_metrics,
                                      bulk_normal_islands_metrics),
    bulk_genes =     combine_datasets(bulk_tumor_genes_metrics,
                                      bulk_normal_genes_metrics),
    bulk_promoters = combine_datasets(bulk_tumor_promoters_metrics,
                                      bulk_normal_promoters_metrics),
    bulk_enhancers = combine_datasets(bulk_tumor_enhancers_metrics,
                                     bulk_normal_enhancers_metrics),
    list_a = combine_datasets(tumor_list_a_metrics, normal_list_a_metrics),
    list_b = combine_datasets(tumor_list_b_metrics, normal_list_b_metrics),

    ## Normalize all variables to common scale
    scaler_recipe = scaler(bulk_islands),

    normalized_islands =   apply_scaler(scaler_recipe, bulk_islands),
    normalized_genes =     apply_scaler(scaler_recipe, bulk_genes),
    normalized_promoters = apply_scaler(scaler_recipe, bulk_promoters),
    normalized_enhancers = apply_scaler(scaler_recipe, bulk_enhancers),
    normalized_list_a =    apply_scaler(scaler_recipe, list_a),
    normalized_list_b =    apply_scaler(scaler_recipe, list_b),

    ## Calculate variability
    bulk_islands_variability =   calc_variability(normalized_islands),
    bulk_genes_variability =     calc_variability(normalized_genes),
    bulk_promoters_variability = calc_variability(normalized_promoters),
    bulk_enhancers_variability = calc_variability(normalized_enhancers),
    list_a_variability = calc_variability(normalized_list_a),
    list_b_variability = calc_variability(normalized_list_b),

    ## Variability charts
    islands_variability_chart =   variability_chart(bulk_islands_variability,   "islands", min_obs),
    genes_variability_chart =     variability_chart(bulk_genes_variability,     "genes", min_obs),
    promoters_variability_chart = variability_chart(bulk_promoters_variability, "promoters", min_obs),
    enhancers_variability_chart = variability_chart(bulk_enhancers_variability, "enhancers", min_obs),
    list_a_variability_chart = variability_chart(list_a_variability, "List A", min_obs),
    list_b_variability_chart = variability_chart(list_b_variability, "List B", min_obs),

    ## Generate report with charts
    report = target(
      command = {
        rmarkdown::render(knitr_in("doc/analysis.Rmd"), output_format = "all")
        file_out("doc/analysis.html")
        file_out("doc/analysis.pdf")
      }
    ),

    ## Extract data resuls
    most_variable_islands = write_csv(
      extract_most_variable(bulk_islands_variability, n_most),
      file_out("doc/most_variable_islands.csv")
    ),
    most_variable_genes = write_csv(
      extract_most_variable(bulk_genes_variability, n_most),
      file_out("doc/most_variable_genes.csv")
    ),
    most_variable_promoters = write_csv(
      extract_most_variable(bulk_promoters_variability, n_most),
      file_out("doc/most_variable_promoters.csv")
    ),

    # --------------------------------------------------------------------------
    # Islands 150bp ranges

    tumor_island150 = bulk_tumor_idat %>%
      group_by(chr) %>%
      bed_slide(vars = list(`200511490070_R02C01`, `200511490070_R03C01`,
                            `200511490070_R05C01`, `200511490070_R06C01`,
                            `200511490070_R07C01`, `200511490070_R08C01`,
                            `200511490073_R02C01`, `200511490073_R03C01`,
                            `200511490073_R04C01`, `200511490073_R05C01`,
                            `200511490073_R06C01`, `200511490073_R07C01`,
                            `200511490073_R08C01`, `200360140022_R01C01`,
                            `200360140022_R02C01`, `200360140022_R03C01`,
                            `200360140022_R04C01`, `200514040139_R01C01`,
                            `200514040139_R02C01`, `200514040139_R03C01`,
                            `200514040139_R04C01`, `200514030126_R02C01`,
                            `200514030126_R03C01`, `200705860031_R02C01`,
                            `200705860031_R03C01`, `200705860031_R04C01`,
                            `200705860031_R05C01`, `200705860031_R06C01`,
                            `200705860031_R07C01`, `200705860031_R08C01`,
                            `200357150201_R07C01`, `200357150201_R08C01`),
                funs = list(mhic = metric_mhic,
                            pwd = metric_pwd,
                            sample_var = metric_sample_var,
                            n_obs = nrow),
                .i = pos,
                size = 150),
)
