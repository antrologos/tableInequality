#' @export

calc_quantile_MCIB <- function(p, data_pnad, groups = NULL, known_groupMeans = NULL, topBracket_method = c("gpinter","RPME")){

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

        known_groupMeans_checked <- tableInequality:::check_known_groupMeans_DF(data_pnad = data_pnad,
                                                                                groups = groups,
                                                                                known_groupMeans = known_groupMeans)

        data_split <- split(data_pnad, f = data_pnad$ID)

        topBracket_method_chosen = topBracket_method[1]

        #data_i = data_split[[1]]

        quantile_MCIB = function(p, data_i, topBracket_method_chosen, known_groupMeans_checked){

                ID_i = data_i$ID %>% unique()
                lower_i = data_i$min_faixa
                upper_i = data_i$max_faixa
                n_i     = data_i$n

                N = sum(n_i)

                lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
                upper_i = c(lower_i[-1], NA)

                slope_intercept <- tableInequality:::estimate_slope_and_intercept_MCIB(lower_i = lower_i,
                                                                                       upper_i = upper_i,
                                                                                       n_i = n_i)
                m = slope_intercept$m
                c = slope_intercept$c


                quantileFunctionPareto_lastBracket = tableInequality:::make_quantileFunctionPareto_lastBracket(data_i = data_i,
                                                                                                               topBracket_method_chosen = topBracket_method_chosen,
                                                                                                               known_groupMeans_checked = known_groupMeans_checked,
                                                                                                               m = m,
                                                                                                               c = c)

                pareto_upper_bound = tableInequality:::get_pareto_upper_bound(data_i = data_i,
                                                                              topBracket_method_chosen = topBracket_method_chosen,
                                                                              known_groupMeans_checked = known_groupMeans_checked,
                                                                              m = m,
                                                                              c = c)

                quantile_function_MCIB = function(p){

                        N = sum(n_i)

                        cum_p = cumsum(n_i)/N
                        cum_p = c(0,cum_p[-length(cum_p)])

                        quantile_data = tibble(lower_i,
                                               upper_i,
                                               cum_p_min = cum_p,
                                               cum_p_max = c(cum_p[-1],1))

                        i = sapply(p, function(x){
                                which((x >= quantile_data$cum_p_min) & (x < quantile_data$cum_p_max))
                        })

                        i = as.numeric(i)

                        # Quantiles for closed brackets with non-zero slopes
                        m_i = m[i]
                        c_i = c[i]
                        min_i = lower_i[i]
                        diff_p = p - quantile_data$cum_p_min[i]

                        a1 = (m_i/(2*N))
                        a2 = (c_i/N)
                        a3 = -( a1*(min_i^2) + a2*min_i + diff_p)

                        y0 = (-a2 + sqrt(a2^2 - 4*a1*a3))/(2*a1)
                        y1 = (-a2 - sqrt(a2^2 - 4*a1*a3))/(2*a1)

                        test_y0 <- round(y0,9) >= round(quantile_data$lower_i[i],9) & round(y0,9) < round(quantile_data$upper_i[i],9)
                        test_y1 <- round(y1,9) >= round(quantile_data$lower_i[i],9) & round(y1,9) < round(quantile_data$upper_i[i],9)

                        quantile_closedBracket = ifelse(test_y0 == T, y0, NA)
                        quantile_closedBracket = ifelse(test_y1 == T, y1, quantile_closedBracket)

                        # Quantiles for closed brackets with uniform densities
                        y_uniform = (p - quantile_data$cum_p_min[i])*(N/c_i) + min_i
                        quantile_closedBracket = ifelse(m_i == 0 & c_i >0,
                                                        y_uniform,
                                                        quantile_closedBracket)

                        if(!is.null(quantileFunctionPareto_lastBracket)){
                                quantile_openBracket = quantileFunctionPareto_lastBracket(p)
                        }else{
                                quantile_openBracket <- rep(NA, length(p))
                        }

                        quantile = ifelse(i < nrow(quantile_data),
                                          quantile_closedBracket,
                                          quantile_openBracket)

                        if(last(n_i) == 0){
                                last_valid <- last(which(n_i > 0))
                                quantile[p==1] = upper_i[last_valid]
                        }

                        quantile = ifelse(p < 0 | p > 1, NA, quantile)

                        quantile
                }

                quantiles_i = quantile_function_MCIB(p)
                as_tibble(matrix(quantiles_i, nrow = 1)) %>% setNames(p)

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        quantile_result <- future_map_dfr(.x = data_split,
                                          .f = quantile_MCIB,
                                          .progress = T,
                                          p = p,
                                          topBracket_method_chosen = topBracket_method_chosen,
                                          known_groupMeans_checked = known_groupMeans_checked,
                                          .options = future_options(globals = c("known_groupMeans_checked",
                                                                                "topBracket_method_chosen"),
                                                                    packages = c("tableInequality", "data.table"))) %>%
                bind_cols(ID = names(data_split), .)

        if(is.null(groups)){
                quantile_result <- quantile_result %>%
                        dplyr::select(-ID)
        }else{
                quantile_result <- quantile_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        quantile_result
}



