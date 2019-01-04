#' @export

calc_mean_logNormal <- function(data_pnad,
                                groups = NULL){

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

        #data_i = data_split[[8]]
        mean_logNormal = function(data_i){

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

                        sigma <- sqrt(sigma2)

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

                mu     = parameters$estimate[1]
                sigma2 = exp(parameters$estimate[2])
                sigma4 = sigma2^2

                #correction factor
                #https://stats.stackexchange.com/questions/221465/why-is-the-arithmetic-mean-smaller-than-the-distribution-mean-in-a-log-normal-di
                cf = ((exp(sigma2) - 1)/(sigma2 + sigma4/2))

                sigma2_corrected = cf*sigma2

                mean = exp(mu + sigma2_corrected/2)

                mean
        }

        mean_result <- map_df(data_split, mean_logNormal)

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
