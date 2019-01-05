#' @export

calc_gini_MCIB <- function(data_pnad, groups = NULL, known_groupMeans = NULL){

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

        known_groupMeans <- check_known_groupMeans_DF(data_pnad = data_pnad,
                                                      groups = groups,
                                                      known_groupMeans = known_groupMeans)

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[1]]

        gini_MCIB = function(data_i){

                ID_i = data_i$ID %>% unique()
                lower_i = data_i$min_faixa
                upper_i = data_i$max_faixa
                n_i     = data_i$n

                N = sum(n_i)

                lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
                upper_i = c(lower_i[-1], NA)

                slope_intercept <-estimate_slope_and_intercept_MCIB(lower_i = lower_i,
                                                                    upper_i = upper_i,
                                                                    n_i = n_i)
                m = slope_intercept$m
                c = slope_intercept$c

                beta_pareto = last(data_i$min_faixa)
                if(!is.null(known_groupMeans)){

                        knownMean_i = known_groupMeans %>%
                                filter(ID == ID_i) %>%
                                .$mean

                        group_means_by_integral = {
                                ((m/(3*n_i))*(upper_i^3) + (c/(2*n_i))*(upper_i^2)) - ((m/(3*n_i))*(lower_i^3) + (c/(2*n_i))*(lower_i^2))
                        }

                        mean_last_group = (1/last(n_i))*(N*knownMean_i - sum(n_i*group_means_by_integral, na.rm = T))


                        if(mean_last_group > beta_pareto){
                                alpha_pareto = mean_last_group/(mean_last_group - beta_pareto)
                        }else{
                                alpha_pareto <- getMids(ID = ID_i,
                                                        hb = n_i,
                                                        lb = lower_i,
                                                        ub = upper_i,
                                                        alpha_bound = 2)$alpha
                        }

                }else{
                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha

                }

                pareto_upper_bound = exp( log(beta_pareto) - log(1 - 0.995)/alpha_pareto)

                if(is.na(pareto_upper_bound)){
                        pareto_upper_bound = Inf
                }

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
                                if(!is.na(alpha_pareto)){
                                        beta_pareto = lower_i[n_brackets]
                                        log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                                        exp(log_PDF)
                                }else{
                                        0
                                }
                        }

                        probDensity_openBracket = (n_i[n_brackets]/N)*f_pareto_lastBracket(y)

                        probDensity = ifelse(i < n_brackets,
                                             probDensity_closedBrackets,
                                             probDensity_openBracket)

                        probDensity = ifelse(y < lower_i[1], 0, probDensity)

                        if(is.na(alpha_pareto)){
                                probDensity[y > last(lower_i)] <- 0
                        }

                        probDensity
                }


                cdf_MCIB = function(y){

                        n_brackets = length(lower_i)
                        min = lower_i[1]

                        i = sapply(y, function(x){
                                which((x >= lower_i) & (x < upper_i))
                        })

                        i = as.numeric(i)
                        i[y >= lower_i[n_brackets]] <- n_brackets

                        cum_p = cumsum(n_i)/N
                        cum_p = c(0,cum_p[-n_brackets])

                        m_i = m[i]
                        c_i = c[i]
                        min_i = lower_i[i]
                        cum_p_i = cum_p[i]

                        prob_cum_closed = cum_p_i + (m_i/(2*N))*(y^2) + (c_i/N)*y - (m_i*(min_i^2))/(2*N) - (c_i/N)*min_i

                        f_pareto_lastBracket = function(y){

                                if(!is.na(alpha_pareto)){
                                        CDF_pareto = 1 - (beta_pareto/y)^alpha_pareto
                                        last(cum_p) + (n_i[n_brackets]/N)*CDF_pareto
                                }else{
                                        last(cum_p)
                                }
                        }

                        prob_cum_open = f_pareto_lastBracket(y)

                        prob_cum = ifelse(i < n_brackets,
                                          prob_cum_closed,
                                          prob_cum_open)

                        prob_cum = ifelse(y < lower_i[1], 0, prob_cum)

                        prob_cum

                }


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
                        quantile_closedBracket = ifelse(m_i == 0, y_uniform, quantile_closedBracket)


                        f_pareto_quantile = function(p, i){
                                beta_pareto = last(quantile_data$lower_i)
                                y =  beta_pareto/((1 - (p - last(cum_p))/(last(n_i)/N))^(1/alpha_pareto)) #quantile function for pareto
                                y
                        }

                        quantile_openBracket = f_pareto_quantile(p, i)

                        quantile = ifelse(i < nrow(quantile_data),
                                          quantile_closedBracket,
                                          quantile_openBracket)

                        quantile = ifelse(p < 0 | p > 1, NA, quantile)

                        if(is.na(alpha_pareto)){
                                quantile[p==1] <- last(lower_i)
                        }

                        quantile
                }

                if(!is.null(knownMean_i)){
                        grand_mean = knownMean_i
                }else{
                        if(is.na(alpha_pareto)){
                                pareto_upper_bound = last(lower_i)
                        }

                        grand_mean = pracma::integral(function(x) x*pdf_MCIB(x),
                                                      xmin = first(lower_i),
                                                      xmax = pareto_upper_bound)
                }

                lorenz_MCIB = function(x) {

                        f = function(y) (1/grand_mean)*y*pdf_MCIB(y)

                        # lorenz value for one observation
                        lorenz_i = function(z) integrate(f = f,
                                                         lower = first(lower_i),
                                                         upper = z,
                                                         rel.tol = .Machine$double.eps^0.275,
                                                         stop.on.error = FALSE,
                                                         subdivisions = 2000)$value

                        # lorenz value for a vector
                        future_map_dbl(x, lorenz_i)
                }

                p_max = cdf_MCIB(pareto_upper_bound)

                lorenz_integral = integrate(f = function(x) lorenz_MCIB(quantile_function_MCIB(x)),
                                            lower = 0,
                                            upper = p_max,
                                            rel.tol = .Machine$double.eps^0.25,
                                            stop.on.error = FALSE,
                                            subdivisions = 2000)$value

                gini = 1 - 2*lorenz_integral

                gini

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        gini_result <- future_map_dbl(.x = data_split,
                                          .f = gini_MCIB,
                                          .progress = T) %>%
                tibble(ID = names(.), gini = .)

        if(is.null(groups)){
                gini_result <- gini_result %>%
                        dplyr::select(-ID)
        }else{
                gini_result <- gini_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        gini_result

}



