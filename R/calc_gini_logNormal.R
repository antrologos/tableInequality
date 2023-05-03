#' @export


calc_gini_logNormal <- function(data_pnad,
                                groups = NULL){

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
        gini_logNormal = function(data_i){

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

                likelihood <- function(logNormalParameters, log_min, log_max, n){

                        mu     <- logNormalParameters[1]
                        sigma2 <- exp(logNormalParameters[2])

                        sigma <- sqrt(sigma2)

                        probs <- pnorm(log_max - mu, sd = sigma) - pnorm(log_min - mu, sd = sigma)
                        probs[length(probs)] = 1 - pnorm(last(log_min) - mu, sd = sigma)

                        sum(n*log(probs)) #negativo porque o nlm minimiza
                }

                neg_likelihood <- function(logNormalParameters, log_min, log_max, n){
                        -likelihood(logNormalParameters, log_min, log_max, n)
                }

                parameters <- try( maxLik::maxLik(logLik = likelihood,
                                                  start = c(1,1),
                                                  log_min = data_i$log_min,
                                                  log_max = data_i$log_max,
                                                  n = data_i$n),
                                   silent = TRUE)

                if("try-error" %in% class(parameters)){
                        parameters <- nlm(f = neg_likelihood,
                                          p = c(1,1),
                                          log_min = data_i$log_min,
                                          log_max = data_i$log_max,
                                          n = data_i$n)
                }else{
                        if(parameters$code == 3){
                                parameters <- maxLik::maxLik(logLik = likelihood,
                                                             start = c(1,1),
                                                             log_min = data_i$log_min,
                                                             log_max = data_i$log_max,
                                                             n = data_i$n,
                                                             method = "BFGS")
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

                gini_lognormal = erf(x = sigma/2)
                gini_lognormal
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        gini_result <- map_dfr(data_split, gini_logNormal)

        gini_result <- tibble(ID   = rownames(t(gini_result)),
                              gini = t(gini_result)[,1])

        if(is.null(groups)){
                gini_result <- gini_result %>%
                        dplyr::select(gini)
        }else{
                gini_result <- gini_result %>%
                        dplyr::select(ID, gini) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        gini_result
}
