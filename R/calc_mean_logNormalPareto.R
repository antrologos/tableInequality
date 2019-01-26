#' @export


calc_mean_logNormalPareto <- function(data_pnad,
                                      groups = NULL,
                                      limite_distribuicoes = .9){

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


        #problems
        #data_i =  data_split[[ 25 ]]
        #data_i =  data_split[[2229]]

        #which(names(data_split) == 316590)

        mean_loglinPareto = function(data_i, grid_mean){

                #for(i in 1:length(data_split)){

                #       print(i)

                #data_i = data_split[[i]]

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


                #===========================================================
                # PASSO 2 - Interpolação Log-normal


                likelihood <- function(logNormalParameters){

                        mu     <- logNormalParameters[1]
                        sigma2 <- exp(logNormalParameters[2])

                        sigma  <- sqrt(sigma2)

                        with(data_i, {
                                probs <- pnorm(log_max - mu, sd = sigma) - pnorm(log_min - mu, sd = sigma)
                                probs[length(probs)] = 1 - pnorm(last(log_min) - mu, sd = sigma)

                                -sum(n*log(probs)) #negativo porque o nlm minimiza
                        })
                }

                parameters <- try( maxLik(logLik = function(x) -likelihood(x),
                                          start = c(1,1)),
                                   silent = TRUE)

                if("try-error" %in% class(parameters)){
                        parameters <- nlm(f = likelihood, p = c(1,1))
                }else{
                        if(parameters$code == 3){
                                parameters <- maxLik::maxLik(logLik = function(x) -likelihood(x),
                                                             start = c(1,1), method = "BFGS")
                        }
                }

                mu     = parameters$estimate[1]
                sigma2 = exp(parameters$estimate[2])
                sigma4 = sigma2^2

                #correction factor
                #https://stats.stackexchange.com/questions/221465/why-is-the-arithmetic-mean-smaller-than-the-distribution-mean-in-a-log-normal-di
                cf = ((exp(sigma2) - 1)/(sigma2 + sigma4/2))

                sigma2_corrected = cf*sigma2
                sigma  = sqrt(sigma2_corrected)


                pdf_lognormal      <- function(y) dlnorm(x = y,meanlog = mu, sdlog = sigma)
                cdf_lognormal      <- function(y) plnorm(q = y,meanlog = mu, sdlog = sigma)
                quantile_lognormal <- function(p) qlnorm(p = p,meanlog = mu, sdlog = sigma)

                p_cum_maxValue_lognormal <- cdf_lognormal(max_value)

                pdf_lognormal_adj <- function(y){
                        density = pdf_lognormal(y)/p_cum_maxValue_lognormal
                        density = ifelse(y > max_value, 0, density)
                        density
                }

                cdf_lognormal_adj <- function(y){
                        p = cdf_lognormal(y)/p_cum_maxValue_lognormal
                        p = ifelse(y > max_value, 1, p)
                        p
                }

                if(is.finite(max_value)){
                        quantile_lognormal_adj = Vectorize(tableInequality:::inverse(f = cdf_lognormal_adj,
                                                                                     lower = 0,
                                                                                     upper = max_value,
                                                                                     extendInt = "yes"))
                }else{
                        quantile_lognormal_adj = quantile_lognormal
                }

              #=============================================================================================================
                # Passo 3 - Combinando distribuições

                #y = seq(0, 10000, 100)

                pdf_combined <- function(y){
                        threashold        <- quantile_lognormal_adj(limite_distribuicoes)
                        survival_pareto   <- 1 - cdf_pareto_adj(threashold)
                        correction_factor <- (1 - limite_distribuicoes)/survival_pareto

                        density <- pdf_lognormal_adj(y)
                        density[y > threashold] <- pdf_pareto_adj(y[y > threashold]) * correction_factor

                        density
                }

                #system.time({
                #        mean <- integrate(f = function(y) y*pdf_combined(y),
                #                          lower = 0,
                #                          upper = max_value,
                #                          subdivisions = 2000,
                #                          stop.on.error = F,
                #                          rel.tol = .Machine$double.eps^.5
                #                          )$value
                #})

                #system.time({
                        grid_meanCopy <- grid_mean
                        rescale(grid_meanCopy, domain = c(0, max_value))
                        mean <- mvQuad::quadrature(f = function(y) y*pdf_combined(y),
                                           grid = grid_meanCopy)
                #})
                        #if(str_detect(mean$message, "divergent")){
                        #        mean <- integral(function(y) y*pdf_combined(y), xmin = 0, xmax = max_value)
                        #}else{
                        #        mean <- mean$value
                        #}


                mean
        }

        #teste = NULL
        #for(i in 1:length(data_split)){
        #        print(i)
        #        teste[[i]] = mean_loglinPareto(data_split[[i]], grid_mean )
        #}


        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        grid_mean = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 1500)
        mean_result <- future_map_dfr(.x = data_split,
                              .f = mean_loglinPareto,
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


