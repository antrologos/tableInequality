#' @export

calc_mean_logNormalPareto <- function(data_pnad,
                                      groups = NULL,
                                      limite_distribuicoes = .9,
                                      n_quantiles = 100000){

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
        mean_loglinPareto = function(data_i){

                #for(i in 1:length(data_split)){

                #       print(i)

                #data_i = data_split[[i]]

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                        #next
                }

                # PASSO 1 - Interpolação Log-normal

                data_i <- data_i %>%
                        mutate(log_min = log(min_faixa),
                               log_max = log(max_faixa))

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
                }else if(parameters$code == 3){
                        parameters <- nlm(f = likelihood, p = c(1,1))
                }

                mu    = parameters$estimate[1]
                sigma2 = exp(parameters$estimate[2])
                sigma4 = sigma2^2

                #correction factor
                #https://stats.stackexchange.com/questions/221465/why-is-the-arithmetic-mean-smaller-than-the-distribution-mean-in-a-log-normal-di
                cf = ((exp(sigma2) - 1)/(sigma2 + sigma4/2))

                sigma2_corrected = cf*sigma2

                sigma  = sqrt(sigma2_corrected)

                quantis_logNormal <- tibble(p_cum             = seq(0, 1, length.out = n_quantiles),
                                            quantil_logNormal = exp(qnorm(p_cum, mean = mu, sd = sigma))) %>%
                        mutate(quantil_logNormal     = ifelse(p_cum == 1, NA, quantil_logNormal))


                #===========================================================
                # Passo 2 - Interpolação de Pareto Local

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
                        tibble(theta, k, groups = 1:length(k))
                })


                quantis_Pareto <- tibble(p_cum = seq(0, 1, length.out = n_quantiles),
                                         groups  = as.numeric(cut(p_cum,
                                                                  breaks = c(-1, data_i$p_sup)))) %>%
                        left_join(y = pareto_parameters, by = "groups") %>%
                        mutate(quantil_Pareto     = k*((1-p_cum)^(-1/theta)),
                               quantil_Pareto     = ifelse(p_cum == 1, NA, quantil_Pareto)) %>%
                        dplyr::select(p_cum, quantil_Pareto)



                interpolation_data <- bind_rows(quantis_logNormal %>%
                                                        filter(p_cum < limite_distribuicoes) %>%
                                                        rename(quantil = quantil_logNormal),
                                                quantis_Pareto %>%
                                                        filter(p_cum >= limite_distribuicoes) %>%
                                                        rename(quantil = quantil_Pareto)) %>%
                        mutate(quantil_cum = cumsum(quantil),
                               prop_quantil_cum = (quantil_cum/max(quantil_cum, na.rm = T)),
                               diff_p_cum  = c(0,diff(p_cum)))


                mean <- with(interpolation_data, {
                        sum(quantil*diff_p_cum, na.rm=T)
                })

                mean
        }

        mean_result <- map_df(data_split, mean_loglinPareto)

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

        gc()

        mean_result

}
