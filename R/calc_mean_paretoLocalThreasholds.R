#' @export

calc_mean_paretoLocalThreasholds <- function(data_pnad, groups = NULL){

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

        #data_i = data_split[[1]]

        mean_paretoLocal = function(data_i,  grid_mean){

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                        #next
                }

                data_i <- data_i %>%
                        arrange(min_faixa) %>%
                        mutate(log_min = log(min_faixa),
                               log_max = log(max_faixa),
                               p     = n/sum(n),
                               p_inf = c(0, cumsum(p)[-length(cumsum(p))]),
                               p_sup = cumsum(p))

                max_value <- data_i$max_faixa[first(which(data_i$p_sup == 1))]
                max_value <- ifelse(is.na(max_value), Inf, max_value)

                #===========================================================
                # Passo 1 - Interpolação de Pareto Local

                #data_ii <- data_i %>%
                #        filter(n > 0)

                nrows_max <- first(which(round(data_i$p_sup, 15) == 1))
                data_ii <- data_i[1:nrows_max,] %>%
                        filter(n > 0)

                pareto_parameters <- with(data_ii, {
                        theta = (log(1 - p_inf) - log(1 - p_sup))/(log(max_faixa) - log(min_faixa))
                        k     = ( (p_sup - p_inf)/( (1/min_faixa)^theta - (1/max_faixa)^theta )  )^(1/theta)

                        if(is.na(last(theta))|is.nan(last(theta))|!is.finite(last(theta))){
                                theta[length(theta)] = theta[length(theta)-1]
                                k[length(k)] = k[length(k)-1]
                        }
                        tibble(theta, k, groups = 1:length(k))
                })

                data_pareto <- bind_cols(pareto_parameters,
                                         data_ii %>% dplyr::select(min_faixa, max_faixa)) %>%
                        arrange(min_faixa)

                two_point_theta <- binequality::getMids(ID = "1",
                                                        hb = data_ii$n,
                                                        lb = data_ii$min_faixa,
                                                        ub = data_ii$max_faixa)$alpha

                data_pareto$theta[is.na(data_pareto$max_faixa)] <- two_point_theta
                data_pareto$k[is.na(data_pareto$max_faixa)]     <- data_pareto$min_faixa[is.na(data_pareto$max_faixa)]

                #y = 0:10000

                theta_test = data_pareto$theta[is.na(data_pareto$max_faixa)]
                k_test = data_pareto$k[is.na(data_pareto$max_faixa)]

                pdf_lastBracket = function(y) (theta_test*(k_test^theta_test))/(y^(theta_test+1))
                #constant <- integrate(pdf_lastBracket, lower = k_test, upper = Inf)$value

                correction_factor = last(data_i$p) #/constant

                data_pareto$correction_factor <- 1
                data_pareto$correction_factor[is.na(data_pareto$max_faixa)] <- correction_factor

                data_pareto$cumulative_factor <- 0
                data_pareto$cumulative_factor[is.na(data_pareto$max_faixa)] <- 1 - last(data_i$p)

                pdf_pareto  <- function(y){

                        group_y = map_dbl(y , function(x) {
                                group = last(which(x >= data_pareto$min_faixa))
                                ifelse(length(group) == 0, NA, group)
                        })

                        alpha = data_pareto$theta[group_y]
                        y_min = data_pareto$k[group_y]
                        correction = data_pareto$correction_factor[group_y]

                        density = (alpha*(y_min^alpha))/(y^(alpha+1))
                        density = density * correction

                        ifelse(is.na(density), 0, density)
                }

                cdf_pareto = function(y){
                        group_y = map_dbl(y , function(x) {
                                group = last(which(x >= data_pareto$min_faixa))
                                ifelse(length(group) == 0, NA, group)
                        })

                        alpha = data_pareto$theta[group_y]
                        y_min = data_pareto$k[group_y]
                        correction1 = data_pareto$correction_factor[group_y]
                        correction2 = data_pareto$cumulative_factor[group_y]

                        p = 1 - (y_min/y)^alpha
                        p = p*correction1 + correction2

                        ifelse(is.na(p), 0, p)
                }


                if(is.finite(max_value)){
                        p_cum_maxValue_pareto = cdf_pareto(max_value)
                }else{
                        k_paretoLast     = last(data_pareto$k)
                        theta_paretoLast = last(data_pareto$theta)

                        # New maximum value
                        max_value        = exp( log(k_paretoLast) - log(1 - 0.99)/theta_paretoLast)

                        p_cum_maxValue_pareto = cdf_pareto(max_value)
                }

                pdf_pareto_adj <- function(y){
                        density = pdf_pareto(y)/p_cum_maxValue_pareto
                        density = ifelse(y > max_value, 0, density)
                        density
                }

                cdf_pareto_adj <- function(y){
                        p = cdf_pareto(y)/p_cum_maxValue_pareto
                        p = ifelse(y > max_value, 1, p)
                        p
                }

                quantile_pareto_adj = Vectorize(tableInequality:::inverse(cdf_pareto_adj, lower = 0, upper = max_value, extendInt = "yes"))


                #=====================================================================

                #system.time({
                grid_meanCopy <- grid_mean
                rescale(grid_meanCopy, domain = c(0, max_value))
                mean <- mvQuad::quadrature(f = function(y) y*pdf_pareto_adj(y),
                                                 grid = grid_meanCopy)
                #})

                mean

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        grid_mean = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 1500)
        mean_result <- future_map_dfr(.x = data_split,
                                      .f = mean_paretoLocal,
                                      grid_mean = grid_mean,
                                      .progress = T)

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

