#' @export

calc_mean_MCIB <- function(data_pnad, groups = NULL){

        data_pnad <- data_pnad %>%
                unite(col = ID, groups) %>%
                group_by(ID, faixas_renda) %>%
                summarise(min_faixa = min(min_faixa),
                          max_faixa = max(max_faixa),
                          n         = sum(n)) %>%
                ungroup() %>%
                arrange(ID, min_faixa)

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[1]]

        mean_MCIB = function(data_i){

                ID = data_i$ID %>% unique()
                lower_i = data_i$min_faixa
                upper_i = data_i$max_faixa
                n_i     = data_i$n

                N = sum(n_i)

                lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
                upper_i = c(lower_i[-1], NA)

                slope_intercept <- estimate_slope_and_intercept_MCIB(lower_i = lower_i,
                                                                     upper_i = upper_i,
                                                                     n_i = n_i)

                m = slope_intercept$m
                c = slope_intercept$c


                alpha_pareto <- getMids(ID = ID,
                                        hb = n_i,
                                        lb = lower_i,
                                        ub = upper_i,
                                        alpha_bound = 2)$alpha

                beta_pareto = last(data_i$min_faixa)
                pareto_upper_bound = exp( log(beta_pareto) - log(1 - 0.9995)/alpha_pareto)

                #y = seq(6.5, 5000, 100)
                pdf_MCIB = function(y){

                        n_brackets = length(lower_i)
                        i = sapply(y, function(x){
                                which((x >= lower_i) & (x < upper_i))
                        })

                        i = as.numeric(i)
                        i[y >= lower_i[n_brackets]] <- n_brackets

                        m_i = m[i]
                        c_i = c[i]

                        probDensity_closedBrackets = (m_i*y + c_i)/N

                        f_pareto_lastBracket = function(y){
                                beta_pareto = lower_i[n_brackets]
                                log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                                exp(log_PDF)
                        }

                        probDensity_openBracket = (n_i[n_brackets]/N)*f_pareto_lastBracket(y)

                        probDensity = ifelse(i < n_brackets,
                                             probDensity_closedBrackets,
                                             probDensity_openBracket)

                        probDensity = ifelse(y < lower_i[1], 0, probDensity)

                        probDensity
                }


                mean = pracma::integral(function(x) x*pdf_MCIB(x),
                                        xmin = first(lower_i),
                                        xmax = pareto_upper_bound)

                mean

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }


        mean_result <- future_map_dbl(.x = data_split,
                                      .f = mean_MCIB,
                                      .progress = T) %>%
                tibble(ID = names(.), mean = .)

        if(is.null(groups)){
                mean_result <- mean_result %>%
                        dplyr::select(-ID)
        }else{
                mean_result <- mean_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        mean_result

}


