get_metrics_by_group <- function(data, basenames, group) {
  bed_group(data = data,
            vars = basenames,
            funs = list(mhic = metric_mhic,
                        pwd = metric_pwd,
                        sample_var = metric_sample_var,
                        n_obs = nrow),
            group = {{group}}) %>%
    dplyr::rename(group = {{group}})
}
