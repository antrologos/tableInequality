#' @export


calc_mean_MCIB <- function(data_pnad, groups = NULL,
                           topBracket_method = c("gpinter","RPME"),
                           firstBracket_flat = TRUE){

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

        topBracket_method_chosen = topBracket_method[1]

        #data_i = data_split[[304]]

        known_groupMeans_checked = NULL

        mean_MCIB = function(data_i, topBracket_method_chosen, grid_mean){

                ID_i = data_i$ID %>% unique()
                lower_i = data_i$min_faixa
                upper_i = data_i$max_faixa
                n_i     = data_i$n

                N = sum(n_i, na.rm = T)

                if(N == 0){
                        return(NA)
                }


                lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
                upper_i = c(lower_i[-1], NA)

                slope_intercept <- tableInequality:::estimate_slope_and_intercept_MCIB(
                        lower_i           = lower_i,
                        upper_i           = upper_i,
                        n_i               = n_i,
                        firstBracket_flat = firstBracket_flat)

                m = slope_intercept$m
                c = slope_intercept$c

                PDFpareto_lastBracket =  tableInequality:::make_PDFpareto_lastBracket(data_i = data_i,
                                                                                      topBracket_method_chosen = topBracket_method_chosen,
                                                                                      known_groupMeans_checked = known_groupMeans_checked,
                                                                                      m = m,
                                                                                      c = c)

                CDFpareto_lastBracket =  tableInequality:::make_CDFpareto_lastBracket(data_i = data_i,
                                                                                      topBracket_method_chosen = topBracket_method_chosen,
                                                                                      known_groupMeans_checked = known_groupMeans_checked,
                                                                                      m = m,
                                                                                      c = c)

                quantileFunctionPareto_lastBracket =  tableInequality:::make_quantileFunctionPareto_lastBracket(data_i = data_i,
                                                                                                                topBracket_method_chosen = topBracket_method_chosen,
                                                                                                                known_groupMeans_checked = known_groupMeans_checked,
                                                                                                                m = m,
                                                                                                                c = c)


                pareto_upper_bound =  tableInequality:::get_pareto_upper_bound(data_i = data_i,
                                                                               topBracket_method_chosen = topBracket_method_chosen,
                                                                               known_groupMeans_checked = known_groupMeans_checked,
                                                                               m = m,
                                                                               c = c)

                pdf_MCIB = function(y){

                        n_brackets = length(lower_i)
                        i = sapply(y, function(x){
                                which((x >= lower_i) & (x < upper_i))
                        })

                        i = as.numeric(i)
                        i[y >= lower_i[n_brackets]] <- n_brackets

                        m_i = m[i]
                        c_i = c[i]

                        probDensity_closedBrackets <- (m_i*y + c_i)/N

                        if(!is.null(PDFpareto_lastBracket)){
                                probDensity_openBracket <- PDFpareto_lastBracket(y)
                        }else{
                                probDensity_openBracket <- rep(0, length(y))
                        }

                        probDensity = ifelse(i < n_brackets,
                                             probDensity_closedBrackets,
                                             probDensity_openBracket)

                        probDensity = ifelse(y < lower_i[1], NA, probDensity)
                        probDensity = ifelse(y > pareto_upper_bound, NA, probDensity)

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

                        if(!is.null(CDFpareto_lastBracket)){
                                prob_cum_open = CDFpareto_lastBracket(y)
                        }else{
                                prob_cum_open <- rep(1, length(y))
                        }

                        prob_cum = ifelse(i < n_brackets,
                                          prob_cum_closed,
                                          prob_cum_open)

                        prob_cum = ifelse(y < lower_i[1], NA, prob_cum)
                        prob_cum = ifelse(y > pareto_upper_bound, NA, prob_cum)

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

                        if(!is.finite(pareto_upper_bound)){
                                quantile = ifelse(p == 1, NA, quantile)
                        }

                        quantile
                }

                grid_meanCopy <- grid_mean

                if(!is.finite(pareto_upper_bound)){
                        upper_bound = quantileFunctionPareto_lastBracket(1 - .Machine$double.eps^.4)
                }else{
                        upper_bound = pareto_upper_bound
                }

                rescale(grid_meanCopy, domain = c(first(lower_i), upper_bound))
                mean <- mvQuad::quadrature(f = function(y) y*pdf_MCIB(y),
                                           grid = grid_meanCopy)
                mean
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        grid_mean = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 10000)


        mean_result <- future_map_dbl(.x = data_split,
                                      .f = mean_MCIB,
                                      topBracket_method_chosen = topBracket_method_chosen,
                                      grid_mean = grid_mean,
                                      .progress = T,
                                      .options = future_options(globals = c("topBracket_method_chosen",
                                                                            "known_groupMeans_checked",
                                                                            "firstBracket_flat",
                                                                            "grid_mean"),
                                                                packages = c("tableInequality", "data.table"))) %>%
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


