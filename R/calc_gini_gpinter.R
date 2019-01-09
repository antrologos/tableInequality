#' @export

calc_gini_gpinter <- function(data_pnad, groups = NULL, known_groupMeans = NULL){

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by(ID, faixas_renda) %>%
                        summarise(min_faixa = min(min_faixa),
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        known_groupMeans = tableInequality:::check_known_groupMeans_DF(data_pnad, groups, known_groupMeans)

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[1]]

        gini_paretoGeneralizada <- function(data_i){

                ID_i = data_i$ID %>% unique()
                p    = data_i$n/sum(data_i$n)

                prob_quantis <- tibble(p_cum = c(0, cumsum(p[-length(p)])),
                                       q     =  data_i$min_faixa)

                prob_quantis <- prob_quantis %>%
                        filter(p_cum < 1)
                min_p <- last(which(prob_quantis$p_cum == 0))
                prob_quantis <- prob_quantis[min_p:nrow(prob_quantis),]

                if(sum(prob_quantis$p_cum > 0) < 3){
                        return(NA)
                }

                if(!is.null(known_groupMeans)){
                        known_groupMean_i = known_groupMeans %>%
                                filter(ID == ID_i) %>%
                                .$mean

                        pareto_threshold_fitted <- try(
                                thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q,
                                               average = known_groupMean_i),
                                silent = T)
                }else{
                        pareto_threshold_fitted <- try(
                                thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q),
                                silent = T)
                }

                if("try-error" %in% class(pareto_threshold_fitted)){
                        return(NA)
                }

                gini(pareto_threshold_fitted)

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        gini_result = future_map_dbl(.x = data_split,
                                     .f = gini_paretoGeneralizada,
                                     .progress = T) %>%
                tibble(ID = names(.), gini = .)


        if(is.null(groups)){
                gini_result <- gini_result %>%
                        dplyr::select(gini)
        }else{
                gini_result <- gini_result %>%
                        dplyr::select(ID, gini) %>%
                separate(col = ID, into = groups, sep = "_")
        }

        gini_result
}
