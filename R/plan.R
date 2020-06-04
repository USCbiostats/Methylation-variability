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

    ## Load in datasets
    bulk_tumor_islands = get_metrics_islands(bulk_tumor_basenames, idat_folder, min_obs),
    bulk_tumor_genes = get_metrics_genes(bulk_tumor_basenames, idat_folder, min_obs),
    bulk_tumor_promoters = get_metrics_promoters(bulk_tumor_basenames, idat_folder, min_obs),
    bulk_tumor_enhancers = get_metrics_enhancers(bulk_tumor_basenames, idat_folder, min_obs),
    tumor_list_a = get_metrics_cpgs(readr::read_rds("data/A.rds"), bulk_tumor_basenames, idat_folder, min_obs),
    tumor_list_b = get_metrics_cpgs(readr::read_rds("data/B.rds"), bulk_tumor_basenames, idat_folder, min_obs),

    bulk_normal_islands = get_metrics_islands(bulk_normal_basenames, idat_folder, min_obs),
    bulk_normal_genes = get_metrics_genes(bulk_normal_basenames, idat_folder, min_obs),
    bulk_normal_promoters = get_metrics_promoters(bulk_normal_basenames, idat_folder, min_obs),
    bulk_normal_enhancers = get_metrics_enhancers(bulk_normal_basenames, idat_folder, min_obs),
    normal_list_a = get_metrics_cpgs(readr::read_rds("data/A.rds"), bulk_normal_basenames, idat_folder, min_obs),
    normal_list_b = get_metrics_cpgs(readr::read_rds("data/B.rds"), bulk_normal_basenames, idat_folder, min_obs),

    # Clean all datasets
    bulk_tumor_islands_region =   clean_datasets(bulk_tumor_islands),
    bulk_tumor_genes_region =     clean_datasets(bulk_tumor_genes),
    bulk_tumor_promoters_region = clean_datasets(bulk_tumor_promoters),
    bulk_tumor_enhancers_region = clean_datasets(bulk_tumor_enhancers),
    tumor_list_a_region =         clean_datasets(tumor_list_a),
    tumor_list_b_region =         clean_datasets(tumor_list_b),

    bulk_normal_islands_region =   clean_datasets(bulk_normal_islands),
    bulk_normal_genes_region =     clean_datasets(bulk_normal_genes),
    bulk_normal_promoters_region = clean_datasets(bulk_normal_promoters),
    bulk_normal_enhancers_region = clean_datasets(bulk_normal_enhancers),
    normal_list_a_region =         clean_datasets(normal_list_a),
    normal_list_b_region =         clean_datasets(normal_list_b),

    # Combine
    bulk_islands =   combine_datasets(bulk_tumor_islands_region,
                                      bulk_normal_islands_region),
    bulk_genes =     combine_datasets(bulk_tumor_genes_region,
                                      bulk_normal_genes_region),
    bulk_promoters = combine_datasets(bulk_tumor_promoters_region,
                                      bulk_normal_promoters_region),
    bulk_enhancers = combine_datasets(bulk_tumor_enhancers_region,
                                      bulk_normal_enhancers_region),
    list_a = combine_datasets(tumor_list_a_region, normal_list_a_region),
    list_b = combine_datasets(tumor_list_b_region, normal_list_b_region),

    # Normalize all variables to common scale
    scaler_recipe = scaler(bulk_islands),

    normalized_islands =   apply_scaler(scaler_recipe, bulk_islands),
    normalized_genes =     apply_scaler(scaler_recipe, bulk_genes),
    normalized_promoters = apply_scaler(scaler_recipe, bulk_promoters),
    normalized_enhancers = apply_scaler(scaler_recipe, bulk_enhancers),
    normalized_list_a =    apply_scaler(scaler_recipe, list_a),
    normalized_list_b =    apply_scaler(scaler_recipe, list_b),

    # Calculate variability
    bulk_islands_variability =   calc_variability(normalized_islands),
    bulk_genes_variability =     calc_variability(normalized_genes),
    bulk_promoters_variability = calc_variability(normalized_promoters),
    bulk_enhancers_variability = calc_variability(normalized_enhancers),
    list_a_variability = calc_variability(normalized_list_a),
    list_b_variability = calc_variability(normalized_list_b),

    # Variability charts
    islands_variability_chart =   variability_chart(bulk_islands_variability,   "islands"),
    genes_variability_chart =     variability_chart(bulk_genes_variability,     "genes"),
    promoters_variability_chart = variability_chart(bulk_promoters_variability, "promoters"),
    enhancers_variability_chart = variability_chart(bulk_enhancers_variability, "enhancers"),
    list_a_variability_chart = variability_chart(list_a_variability, "List A"),
    list_b_variability_chart = variability_chart(list_b_variability, "List B"),

    # Generate report with charts
    report = target(
      command = {
        rmarkdown::render(knitr_in("doc/analysis.Rmd"), output_format = "all")
        file_out("doc/analysis.html")
        file_out("doc/analysis.pdf")
      }
    ),

    # Extract data resuls
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

    tumor_island150 = get_island150_data(bulk_tumor_basenames, idat_folder),

    tumor_island150_mhic = calc_mhic_island150(tumor_island150),
    tumor_island150_pwd = calc_pwd_island150(tumor_island150),
    tumor_island150_var_sample = calc_var_sample_island150(tumor_island150),
    tumor_island150_var_cpg = calc_var_cpg_island150(tumor_island150),

)
