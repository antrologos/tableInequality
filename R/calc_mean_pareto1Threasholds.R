#' @export

calc_mean_pareto1Threasholds <- function(data_pnad, groups = NULL){

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

        data_split <- split(data_pnad, f = data_pnad$ID)

        mean_pareto = function(data_i){

                #i=10
                #for(i in 1:length(data_split)){
                #        print(i)
                #        data_i <- data_split[[i]]

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                        #next
                }

                data_i <- data_i %>%
                        filter(n > 0) %>%
                        arrange(min_faixa) %>%
                        mutate(p     = n/sum(n),
                               p_sup = cumsum(p),
                               p_inf = c(0, p_sup[-length(p_sup)]))

                nrows_max <- first(with(data_i, which(round(p_sup, 15) == 1)))

                data_i <- data_i[1:nrows_max,]

                pareto_parameters <- with(data_i, {
                        theta = (log(1 - p_inf) - log(1 - p_sup))/(log(max_faixa) - log(min_faixa))
                        k     = ( (p_sup - p_inf)/( (1/min_faixa)^theta - (1/max_faixa)^theta )  )^(1/theta)

                        if(is.na(last(theta))|is.nan(last(theta))|!is.finite(last(theta))){
                                theta[length(theta)] = theta[length(theta)-1]
                                k[length(k)] = k[length(k)-1]
                        }
                        tibble(theta, k, grupos = 1:length(k))
                })

                interpolation_data <- tibble(p_cum = seq(0, 1, length.out = 1000),
                                             grupos  = as.numeric(cut(p_cum,
                                                                      breaks = c(-1, data_i$p_sup)))
                ) %>%
                        left_join(y = pareto_parameters, by = "grupos") %>%
                        mutate(quantil     = k*((1-p_cum)^(-1/theta)),
                               quantil     = ifelse(p_cum == 1, NA, quantil),
                               quantil_cum = cumsum(quantil),
                               prop_quantil_cum = (quantil_cum/max(quantil_cum, na.rm = T)),
                               diff_p_cum  = c(0,diff(p_cum)))

                mean = with(interpolation_data, {

                        sum(quantil*diff_p_cum, na.rm = T)
                })
                #}

                mean
        }

        mean_result <- map_df(data_split, mean_pareto)

        mean_result <- tibble(ID   = rownames(t(mean_result)),
                              mean = t(mean_result)[,1])

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
