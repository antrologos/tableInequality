#' @export

calc_mean_gpinter <- function(data_pnad, groups = NULL){

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

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[1]]
        mean_paretoGeneralizada <- function(data_i){

                prob_quantis <- with(data_i, {
                        p     = n/sum(n)

                        tibble(p_cum = c(0, cumsum(p[-length(p)])),
                               q     =  min_faixa)
                })


                prob_quantis <- prob_quantis %>%
                        filter(p_cum < 1)
                min_p <- last(which(prob_quantis$p_cum == 0))
                prob_quantis <- prob_quantis[min_p:nrow(prob_quantis),]

                if(sum(prob_quantis$p_cum > 0) < 3){
                        return(NA)
                }

                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q),
                        silent = T)

                if("try-error" %in% class(pareto_threshold_fitted)){
                        return(NA)
                }

                bottom_average(pareto_threshold_fitted, p = 1)

        }

        mean_result <- map(data_split, mean_paretoGeneralizada) %>%
                tibble(ID = names(.),  mean = unlist(.))

        if(is.null(groups)){
                mean_result <- mean_result %>%
                        dplyr::select(mean)
        }else{
                mean_result <- mean_result %>%
                        dplyr::select(ID, mean) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        mean_result
}
