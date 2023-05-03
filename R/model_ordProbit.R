#' @export


model_ordProbit <- function(formula_rhs,
                            data_pnad,
                            groups = NULL,
                            max_trials = 20){

        if(length(formula_rhs) != 2){
                stop("formula_rhs must cointain only the right hand side (independent variables)")
        }

        indep = all.vars(formula_rhs[[2]])

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
                        group_by_at(c("ID","min_faixa",indep)) %>%
                        summarise(max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_pnad <- tableInequality:::prepare_to_regression(data_pnad, indep)
        gc()

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[5]]

        reg_loglin = function(data_i, formula_rhs, max_trials){

                if(sum(data_i$n) == 0){
                        return(as.numeric(NA))
                }

                ID_i = data_i$ID %>% unique()

                # PASSO 1 -
                data_i <- data_i %>%
                        mutate(log_min = log(min_faixa),
                               log_max = log(max_faixa))

                independentVariablesValues <- NULL
                for(indep_i in indep){
                        independentVariablesValues[[indep_i]] <- unique(data_i[[indep_i]])
                }

                X = model.matrix(formula_rhs, data = data_i)
                w = data_i$n

                N = sum(w)
                k = 2*ncol(X)

                log_min = data_i$log_min
                log_max = data_i$log_max

                # All categories of the independent variables must have valid observations
                # (otherwise the number of estimated coeficients may not be same in each group)
                test_indeps <- tableInequality:::test_for_regression(data_i, indep)

                if(any(test_indeps)){
                        stop(cat("One or more independent variables have empty categories in this group",
                                 "\nGroup:",ID_i,".",
                                 "\nVariables:", paste0(names(which(test_indeps==T)),collapse = ", "))
                             )
                }

                #beta   = rep(1, ncol(X))
                #lambda = rep(1, ncol(X))

                #parameters_to_estimate = c(beta,lambda)

                loglikelihood_function_null <- function(parameters_null){
                        mu     <- parameters_null[1]
                        sigma2 <- exp(parameters_null[2])

                        sigma = sqrt(sigma2)

                        is_max = which(is.na(log_max))

                        probs <- pnorm(log_max - mu, sd = sigma) - pnorm(log_min - mu, sd = sigma)
                        probs[is_max] = 1 - pnorm(log_min[is_max] - mu, sd = sigma)

                        sum(w*log(probs)) #positivo porque o maxLik maximiza
                }
                parameters_null <- try(maxLik(logLik = loglikelihood_function_null,
                                              start  = c(1,1),
                                              method = "CG"),
                                       silent = TRUE)


                if("try-error" %in% class(parameters_null)){
                        trial_number = 1

                        print("Convergence error. Trying again")

                        while(("try-error" %in% class(parameters_null)) & (trial_number <= max_trials) ){
                                parameters_null <- try(
                                        maxLik(logLik = loglikelihood_function_null,
                                               start  = jitter(c(1,1)),
                                               method = "CG"),
                                        silent = T)

                                trial_number = trial_number + 1
                        }

                        if(trial_number > max_trials & ("try-error" %in% class(parameters_null)) ){
                                stop(paste0("Null model failed to converge. Group: ",ID_i))
                        }
                }
                loglikelihood_null = logLik(parameters_null)


                loglikelihood_function <- function(parameters_to_estimate){
                        beta   <- parameters_to_estimate[1:ncol(X)]
                        lambda <- parameters_to_estimate[(ncol(X)+1):length(parameters_to_estimate)]

                        mu     = as.numeric(X%*%beta)
                        sigma2 = exp(as.numeric(X%*%lambda))

                        sigma = sqrt(sigma2)

                        is_max = which(is.na(log_max))

                        probs <- pnorm(log_max - mu, sd = sigma) - pnorm(log_min - mu, sd = sigma)
                        probs[is_max] = 1 - pnorm(log_min[is_max] - mu[is_max], sd = sigma[is_max])

                        sum(w*log(probs)) #positivo porque o maxLik maximiza
                }

                start_values = c(parameters_null$estimate[1], rep(0, ncol(X)-1),
                                 parameters_null$estimate[2], rep(0, ncol(X)-1))

                parameters <- try(maxLik(logLik = loglikelihood_function,
                                         start  = start_values,
                                         method = "CG"),
                                  silent = TRUE)

                if("try-error" %in% class(parameters)){
                        trial_number = 1

                        print("Convergence error. Trying again")

                        new_start_values = jitter(start_values)
                        while(("try-error" %in% class(parameters)) & (trial_number <= max_trials) ){
                                parameters <- try(
                                        maxLik(logLik = loglikelihood_function,
                                               start  = new_start_values,
                                               method = "CG"),
                                        silent = T)

                                if("maxLik" %in% class(parameters)){
                                        if(parameters$code == 100){
                                                parameters <- list()
                                                class(parameters) = "try-error"
                                        }
                                }

                                new_start_values <- jitter(new_start_values)
                                trial_number = trial_number + 1
                        }
                }

                if(("try-error" %in% class(parameters))){
                        stop(paste0("Model failed to converge. Group: ",ID_i))
                }

                loglikelihood = logLik(parameters)



                #if("try-error" %in% class(parameters)){
                #        parameters <- nlm(f = function(x) -likelihood(x),
                #                          p = rep(1, ncol(X)+1) )
                #}

                # -- Main Parameters  --

                beta   = parameters$estimate[1:ncol(X)]
                lambda = parameters$estimate[(ncol(X)+1):length(parameters$estimate)]

                # -- Statistical Inference  --
                #varCov = vcov(parameters)

                sd = stdEr(parameters)

                sd_beta   = sd[1:ncol(X)]
                sd_lambda = sd[(ncol(X)+1):length(parameters$estimate)]

                t_statistic_beta   = beta/sd_beta
                t_statistic_lambda = beta/sd_lambda

                p_value_beta   = 2*(1 - pt(abs(t_statistic_beta),df = N-k))
                p_value_lambda = 2*(1 - pt(abs(t_statistic_lambda),df = N-k))


                # -- Classification Matrix --

                lower_b <- unique(log_min) %>% sort()
                upper_b <- c(unique(log_max) %>% sort(),Inf)
                pred    <- as.numeric(X%*%beta)

                lower_b_matrix <- matrix(lower_b, nrow = length(pred), ncol = length(lower_b), byrow = T)
                upper_b_matrix <- matrix(upper_b, nrow = length(pred), ncol = length(upper_b), byrow = T)

                test_category <- (lower_b_matrix <= pred) & (upper_b_matrix > pred)
                pred_cateogry <- lower_b[apply(test_category, 1, function(x) which(x == T))]

                classification_table <- questionr::wtd.table(log_min, pred_cateogry, w)

                classification_df <- classification_table %>%
                        as_tibble()

                corr_Spearman <- wCorr::weightedCorr(x = classification_df[[1]],
                                                     y = classification_df[[2]],
                                                     method = c("Spearman"),
                                                     weights = classification_df[[3]])

                parameters_beta <- tibble(parameter    = colnames(X),
                                          coefficients = beta,
                                          sd           = sd_beta,
                                          t_statistic  = t_statistic_beta,
                                          p_value      = p_value_beta)

                parameters_lambda <- tibble(parameter    = colnames(X),
                                            coefficients = lambda,
                                            sd           = sd_lambda,
                                            t_statistic  = t_statistic_lambda,
                                            p_value      = p_value_lambda)


                list(ID                      = tibble(ID = ID_i) %>% separate(col = ID, into = groups, sep = "_"),
                     parameters_beta         = parameters_beta,
                     parameters_lambda       = parameters_lambda,
                     model_matrix            = X,
                     independentVariablesValues = independentVariablesValues,
                     weights                 = w,
                     `-2loglik`              = -2*loglikelihood,
                     `-2loglik_null`         = -2*loglikelihood_null,
                     AIC                     = AIC(parameters),

                     pseudoR2_mcfadden       = 1 - loglikelihood/loglikelihood_null,
                     pseudoR2_mcfadden_adj   = 1 - (loglikelihood - k)/loglikelihood_null,
                     pseudoR2_CoxSnell       = 1 - exp(loglikelihood_null - loglikelihood)^(2/N),
                     pseudoR2_Nagelkerke     = (1 - exp(loglikelihood_null - loglikelihood)^(2/N))/(1 - exp(loglikelihood_null)^(2/N)),

                     classification_table    = classification_table,
                     prop_correct_classified = sum(diag(classification_table))/sum(classification_table),
                     R_Spearman              = corr_Spearman,
                     R_Spearman2             = (corr_Spearman)^2,
                     N                       = N,
                     formula                 = formula_rhs)
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        #parameters = list()
        #for(i in 1:length(data_split)){
        #        print(i)
        #        parameters[[i]] = reg_loglin(data_i = data_split[[i]])
        #}

        #regression_result <- do.call(what = bind_rows, parameters)

        regression_result <- future_map(.x          = data_split,
                                        .f          = reg_loglin,
                                        formula_rhs = formula_rhs,
                                        max_trials  = max_trials,
                                        .progress   = T)

        class(regression_result) <- "ordProbitModel"

        #if(is.null(groups)){
        #        regression_result <- regression_result %>%
        #                dplyr::select(-ID)
        #}else{
        #        regression_result <- regression_result %>%
        #                separate(col = ID, into = groups, sep =  "_")
        #}

        regression_result
}



#' @export
model.matrix.ordProbitModel <- function(ordProbitModel, weights = F){
        model_matrix_list = purrr::map(ordProbitModel, function(x) x$model_matrix)

        if(weights == T){
                weights_list = purrr::map(ordProbitModel, function(x) x$weights)

                for(i in 1:length(model_matrix_list)){
                        model_matrix_list[[i]] <- cbind(model_matrix_list[[i]], weights_list[[i]])
                }

        }

        model_matrix_list
}



#' @export
coef.ordProbitModel <- function(ordProbitModel){

        beta_df <- purrr::map_dfr(ordProbitModel, function(x) as_tibble(matrix(x$parameters_beta$coefficients, nrow = 1)))
        lambda_df <- purrr::map_dfr(ordProbitModel, function(x) as_tibble(matrix(x$parameters_lambda$coefficients, nrow = 1)))

        parameter_names <- ordProbitModel[[1]]$parameters_beta$parameter

        beta_matrix   <- beta_df   %>% as.matrix() %>% t()
        lambda_matrix <- lambda_df %>% as.matrix() %>% t()

        colnames(beta_matrix)   <- names(ordProbitModel)
        colnames(lambda_matrix) <- names(ordProbitModel)

        rownames(beta_matrix)   <- parameter_names
        rownames(lambda_matrix) <- parameter_names

        list(beta_matrix   = beta_matrix,
             lambda_matrix = lambda_matrix)
}


#' @export
summary.ordProbitModel <- function(ordProbitModel, digits = 3){

        extract_summary <- function(type){

                coef_matrix <- purrr::map(ordProbitModel, function(x) as_tibble(matrix(x[[paste0("parameters_",type)]]$coefficients, nrow = 1)) %>% setNames(x[[paste0("parameters_",type)]]$parameter))
                sd_matrix   <- purrr::map(ordProbitModel, function(x) as_tibble(matrix(x[[paste0("parameters_",type)]]$sd,           nrow = 1)) %>% setNames(x[[paste0("parameters_",type)]]$parameter))
                p_matrix    <- purrr::map(ordProbitModel, function(x) as_tibble(matrix(x[[paste0("parameters_",type)]]$p_value,      nrow = 1)) %>% setNames(x[[paste0("parameters_",type)]]$parameter))

                coef_matrix <- do.call(bind_rows, coef_matrix) %>% as.matrix() %>% t()
                sd_matrix   <- do.call(bind_rows, sd_matrix)   %>% as.matrix() %>% t()
                p_matrix    <- do.call(bind_rows, p_matrix)    %>% as.matrix() %>% t()

                var_names   <- rownames(coef_matrix)

                AIC                     <- purrr::map_dbl(ordProbitModel, function(x) x$AIC)
                `-2loglik`              <- purrr::map_dbl(ordProbitModel, function(x) x$`-2loglik`)
                pseudoR2_mcfadden_adj   <- purrr::map_dbl(ordProbitModel, function(x) x$pseudoR2_mcfadden_adj)
                prop_correct_classified <- purrr::map_dbl(ordProbitModel, function(x) x$prop_correct_classified)
                N                       <- purrr::map_dbl(ordProbitModel, function(x) x$N)
                R_Spearman              <- purrr::map_dbl(ordProbitModel, function(x) x$R_Spearman)
                R_Spearman2             <- purrr::map_dbl(ordProbitModel, function(x) x$R_Spearman2)

                coef_matrix_char <- format(coef_matrix, digits = digits)
                sd_matrix_char   <- format(sd_matrix,   digits = digits)

                coef_matrix_char <- as_tibble(coef_matrix_char) %>%
                        mutate_all(function(x) x %>% str_replace_all(pattern = "NA", "") %>% str_trim()) %>%
                        as.matrix()

                sd_matrix_char <- as_tibble(sd_matrix_char) %>%
                        mutate_all(function(x) x %>% str_replace_all(pattern = "NA", "") %>% str_trim()) %>%
                        as.matrix()

                sd_matrix_char   <- paste0("(",format(sd_matrix_char, digits = digits),")")
                dim(sd_matrix_char) <- dim(coef_matrix_char)
                sd_matrix_char[is.na(sd_matrix)] <- ""

                star1 <- p_matrix < .1
                star2 <- p_matrix < .05
                star3 <- p_matrix < .01

                star1[is.na(star1)] <- FALSE
                star2[is.na(star2)] <- FALSE
                star3[is.na(star3)] <- FALSE

                coef_matrix_char[star1] <- paste0(coef_matrix_char[star1],"*")
                coef_matrix_char[star2] <- paste0(coef_matrix_char[star2],"*")
                coef_matrix_char[star3] <- paste0(coef_matrix_char[star3],"*")

                results = matrix("",
                                 ncol = ncol(coef_matrix_char),
                                 nrow = 2*nrow(coef_matrix_char))

                uneven_lines = seq(1, 2*nrow(coef_matrix_char), 2)
                even_lines   = uneven_lines + 1

                results[uneven_lines,] <- coef_matrix_char
                results[even_lines,]   <- sd_matrix_char

                results <- rbind(results, "",
                                 AIC,
                                 `-2loglik`,
                                 format(pseudoR2_mcfadden_adj,   digits = digits),
                                 format(prop_correct_classified,   digits = digits),
                                 format(R_Spearman,   digits = digits),
                                 format(R_Spearman2,   digits = digits),
                                 format(N,    digits = 0))


                rownames(results) <- c(rep("", 2*length(var_names)), "", "AIC", "-2loglik", "McFadden Pseudo R2 (Adj)",
                                       "Prop Corr. Class","R_Spearman","R_Spearman2", "N")
                rownames(results)[uneven_lines] <- var_names

                results
        }

        beta_summary   <- extract_summary("beta")
        lambda_summary <- extract_summary("lambda")

        list(beta = beta_summary, lambda = lambda_summary)
}

#' @export
print.ordProbitModel <- function(ordProbitModel){
        print(summary(ordProbitModel))
}

