#' @export

model_ordProbit_defCutsHomoskedastic <- function(formula,
                            data_pnad,
                            groups = NULL){

        dep   = all.vars(formula[[2]])
        indep = all.vars(formula[[3]])

        if(any(indep %in% groups) | any(groups %in% indep)){
                stop("Group variables cannot be used as independent variables in the model")
        }

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <-
                        data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by_at(c("ID",dep,indep)) %>%
                        summarise(max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[24]]

        reg_loglin = function(data_i){

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                }

                ID_i = data_i$ID %>% unique()

                # PASSO 1 -
                data_i <- data_i %>%
                        mutate(log_min = log(min_faixa),
                               log_max = log(max_faixa))

                # Removing empty categories
                for(i in 1:length(indep)){
                        test_i <- data_i %>%
                                group_by_at(indep[i]) %>%
                                summarise(n = sum(n))
                        cat_problem = test_i[[1]][which(test_i$n == 0)]
                        if(length(cat_problem > 0)){
                                data_i <- data_i %>%
                                        filter_at(vars(indep[i]), any_vars(. != cat_problem))
                        }
                }

                X = model.matrix(formula, data = data_i)
                w = data_i$n
                log_min = data_i$log_min
                log_max = data_i$log_max


                #beta = rep(1, ncol(X))
                #parameters_to_estimate = c(beta,1)

                likelihood <- function(parameters_to_estimate){
                        beta  <- parameters_to_estimate[1:ncol(X)]
                        sigma2 <- exp(parameters_to_estimate[(ncol(X)+1):length(parameters_to_estimate)])

                        sigma = sqrt(sigma2)

                        mu = as.numeric(X%*%beta)
                        is_max = which(is.na(log_max))

                        probs <- pnorm(log_max - mu, sd = sigma) - pnorm(log_min - mu, sd = sigma)
                        probs[is_max] = 1 - pnorm(log_min[is_max] - mu[is_max], sd = sigma)

                        sum(w*log(probs)) #positivo porque o maxLik maximiza
                }

                start_values = rep(1, ncol(X)+1)
                parameters <- try(maxLik(logLik = likelihood, start  = start_values), silent = TRUE)

                if("try-error" %in% class(parameters)){
                        trial_number = 1

                        print("Convergence error. Trying again")

                        while(("try-error" %in% class(parameters)) & (trial_number <= 10) ){
                                parameters <- try(
                                        maxLik(logLik = likelihood,
                                               start  = start_values*jitter(sample(start_values))), silent = T)
                        }
                }

                #if("try-error" %in% class(parameters)){
                #        parameters <- nlm(f = function(x) -likelihood(x),
                #                          p = rep(1, ncol(X)+1) )
                #}

                beta  = parameters$estimate[1:ncol(X)]
                sigma2 = exp(parameters$estimate[(ncol(X)+1):length(parameters$estimate)])

                Hessian = pracma::hessian(likelihood, x0 = parameters$estimate)

                varCov = try(solve(-Hessian), silent = T)
                if("try-error" %in% class(varCov)){
                        beta = sd = rep(NA, ncol(X))

                }else{
                        sd = sqrt(diag(varCov))[1:ncol(X)] # the last value is "sigma" (the lonormal dispersion parameter)
                }

                names(beta) = colnames(X)
                names(sd) = paste0("sdError_",colnames(X))

                results = as_tibble(matrix(c(beta,sd, sigma2), nrow = 1))
                names(results) = c(colnames(X), paste0("sdError_",colnames(X)), "sigma2")

                results = bind_cols(tibble(ID = ID_i), results)
                results
        }

        plan(multiprocess)

        #parameters = list()
        #for(i in 1:length(data_split)){
                #print(i)
                #parameters[[i]] = reg_loglin(data_i = data_split[[i]])
        #}
        #
        #regression_result <- do.call(what = bind_rows, parameters)

        regression_result <- future_map_dfr(data_split, reg_loglin, .progress = T)


        if(is.null(groups)){
                regression_result <- regression_result %>%
                        dplyr::select(-ID)
        }else{
                regression_result <- regression_result %>%
                        separate(col = ID, into = groups, sep =  "_")
        }

        regression_result
}


