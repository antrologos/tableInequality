#' @export

calc_mean_MGBE <- function(data_pnad, groups = NULL){

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by(ID, min_faixa) %>%
                        summarise(
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split = split(data_pnad, data_pnad$ID)

        run_GB_family_i = function(data_i){
                function_fitted <- run_GB_family(ID = as.character(data_i$ID),
                                                 hb = data_i$n,
                                                 bin_min = data_i$min_faixa,
                                                 bin_max = data_i$max_faixa,
                                                 obs_mean = NA,
                                                 ID_name = 'ID',
                                                 modelsToFit = c("GB2", "GG", "BETA2",
                                                                 "DAGUM", "SINGMAD",
                                                                 "LOGNO", "WEI", "GA",
                                                                 "LOGLOG", "PARETO2"))

                function_fitted$best_model.filter$bic
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        mean_result <- future_map_dfr(data_split, run_GB_family_i, .progress = T)

        if(is.null(groups)){
                mean_result <- mean_result %>%
                        dplyr::select(estMean) %>%
                        rename(mean = estMean)
        }else{
                mean_result <- mean_result %>%
                        dplyr::select(ID, estMean) %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        rename(mean = estMean)
        }

        mean_result
}

