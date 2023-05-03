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
                        mutate(ID = 1)

                IDs <- unique(data_pnad$ID)
        }else{
                data_pnad <-
                        data_pnad %>%
                        unite(col = ID, groups)

                IDs <- unique(data_pnad$ID)

                data_pnad <- data_pnad %>%
                        group_by_at(c("ID",dep,indep)) %>%
                        summarise(max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_pnad <- prepare_to_regression(data_pnad, indep)
        gc()

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[153]]

        reg_loglin = function(data_i){

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                }

                ID_i = data_i$ID %>% unique()

                # PASSO 1 -
                data_i <- data_i %>%
                        mutate(log_min = log(min_faixa),
                               log_max = log(max_faixa))

                X = model.matrix(formula, data = data_i)
                w = data_i$n
                log_min = data_i$log_min
                log_max = data_i$log_max

                # All categories of the independent variables must have valid observations
                # (otherwise the number of estimated coeficients may not be same in each group)
                test_indeps <- test_for_regression(data_i, indep)

                if(any(test_indeps)){
                        coefs = matrix(rep(as.numeric(NA), ncol(X)*2 + 1) , nrow = 1)
                        colnames(coefs) = c(colnames(X), paste0("sdError_",colnames(X)), "sigma2")
                        coefs = as_tibble(coefs)
                        coefs = bind_cols(tibble(ID = ID_i), coefs)
                        return(coefs)
                }

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

                                trial_number = trial_number + 1
                        }

                        if(trial_number > 10 & ("try-error" %in% class(parameters)) ){
                                coefs = matrix(rep(as.numeric(NA), ncol(X)*2 + 1) , nrow = 1)
                                colnames(coefs) = c(colnames(X), paste0("sdError_",colnames(X)), "sigma2")
                                coefs = as_tibble(coefs)
                                coefs = bind_cols(tibble(ID = ID_i), coefs)
                                return(coefs)
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

        plan(multisession)

        #parameters = list()
        #for(i in 1:length(data_split)){
        #        print(i)
        #        parameters[[i]] = reg_loglin(data_i = data_split[[i]])
        #}

        #regression_result <- do.call(what = bind_rows, parameters)

        regression_result <- future_map_dfr(data_split, reg_loglin, .progress = T)

        regression_result <- left_join(tibble(ID = IDs),
                                       regression_result)


        if(is.null(groups)){
                regression_result <- regression_result %>%
                        dplyr::select(-ID)
        }else{
                regression_result <- regression_result %>%
                        separate(col = ID, into = groups, sep =  "_")
        }

        regression_result
}


