#' @export
varDecomp <- function (model) {
        UseMethod("varDecomp")
}


#' @export
varDecomp.midpointModel <- function(model){
        matrix_beta <- coefficients(model)

        list_X <- model.matrix(model)
        list_e <- residuals(model)
        list_w <- map(.x = names(model), function(x) model[[x]]$weights)

        formula   = model[[1]]$formula
        indep_vars <- all.vars(formula[[3]])

        matrix_lambda = NULL
        lambdas_list = list()
        betas_list = list()

        lambda_start    <- rep(1, nrow(matrix_beta))
        variance_decomp <- NULL
        for(k in 1:length(list_X)){
                print(k)

                X  = list_X[[k]]
                e2 = list_e[[k]]^2
                w  = list_w[[k]]

                model_k  <- glm.fit(x = X, y = e2, weights = w, family = Gamma(link = log),
                                    start = lambda_start)

                lambda <- coefficients(model_k)
                beta   <- matrix_beta[, k]

                p_i = w/sum(w)
                s2_i <- exp(X%*%lambda)

                mu_i <- as.numeric(X%*%beta)
                mu   <- sum(p_i*mu_i)

                between <- sum( p_i*((mu_i - mu)^2) )
                within  <- sum( p_i*s2_i )

                varLog  <- between + within
                between_share <- between/varLog
                within_share  <- within/varLog

                variance_decomp_k <- tibble(ID = names(list_X)[k],
                                            varLog, between, within,
                                            between_share, within_share)
                variance_decomp <- bind_rows(variance_decomp, variance_decomp_k)

                lambda_start = lambda
                matrix_lambda = cbind(matrix_lambda, lambda)

                lambdas_list[[names(list_X)[k]]]$`(Intercept)` <- lambda[1]
                for(indep_var in indep_vars){
                        lambdas_list[[names(list_X)[k]]][[indep_var]] <- lambda[str_detect(names(lambda), indep_var)]
                }

                betas_list[[names(list_X)[k]]]$`"(Intercept)"` <- beta[1]
                for(indep_var in indep_vars){
                        betas_list[[names(list_X)[k]]][[indep_var]] <- beta[str_detect(names(beta), indep_var)]
                }
        }
        colnames(matrix_lambda) = colnames(matrix_beta)

        list_vars = model[[1]]$independentVariablesValues
        data_X <- do.call(expand.grid, list_vars)
        data_X <- data_X %>% mutate_all(as.character)

        formula[2] <- NULL
        X     <- model.matrix(formula, data_X)
        X     <- X[ , colnames(list_X[[1]])]
        IDs_X <- apply(X, 1, function(x) paste(x, collapse = "_"))

        data_weight <- tibble(IDs = IDs_X)
        for(i in 1:length(list_X)){
                data_weight_tmp <- tibble(IDs = apply(list_X[[i]], 1, function(x) paste(x, collapse = "_")),
                                          w   = model[[i]]$weights)

                data_weight_tmp <- data_weight_tmp %>%
                        group_by(IDs) %>%
                        summarise(p  = sum(w)) %>%
                        ungroup() %>%
                        mutate(p = p/sum(p))

                data_weight <- left_join(data_weight,
                                         data_weight_tmp)

                names(data_weight)[i+1] <- names(list_X[i])
        }

        matrix_weight <- as.matrix(data_weight[, -1])
        matrix_weight[is.na(matrix_weight)] <- 0

        #var = names(data_X)[1]
        i = 1

        pContra_list <- list()
        for(var in names(data_X)){
                data_p = data_X[var]

                for(i in 1:ncol(matrix_weight)){
                        ID       <- colnames(matrix_weight)[i]
                        tab      <- colSums(X*matrix_weight[,i])
                        freq_var <- tab[str_detect(names(tab), var)]

                        data_p_var <- tibble(unique(data_X[[var]]), c(1 - sum(freq_var), freq_var)) %>%
                                setNames(c(var, ID))

                        data_p <- left_join(data_p, data_p_var)
                }

                matrix_p <- data_p[-1] %>% as.matrix()
                matrix_p <- (1/(matrix_p/matrix_p[,1]))*matrix_weight # Western & Bloome, p. 313

                pContra_list[[var]] <- matrix_p
        }

        betasContra_list <- NULL
        for(indep_var in indep_vars){
                betasContra_matrix <- matrix(NA,
                                             ncol = length(betas_list),
                                             nrow = length(beta))

                for(i in 1:length(betas_list)){
                        if(i == 1){
                                betasContra1 <- do.call(c, betas_list[[i]])
                                betasContra_matrix[, 1] <- betasContra1
                        }else{

                                betasContra <- do.call(c, betas_list[[i]])
                                betasContra[str_detect(names(betasContra), indep_var)] <- betasContra1[str_detect(names(betasContra), indep_var)]
                                betasContra_matrix[, i] <- betasContra
                        }
                }

                rownames(betasContra_matrix) <- names(beta)
                betasContra_list[[indep_var]] <- betasContra_matrix
        }


        lambdasContra_list <- NULL
        for(indep_var in indep_vars){
                lambdasContra_matrix <- matrix(NA,
                                               ncol = length(lambdas_list),
                                               nrow = length(lambda))

                for(i in 1:length(lambdas_list)){
                        if(i == 1){
                                lambdasContra1 <- do.call(c, lambdas_list[[i]])
                                lambdasContra_matrix[, 1] <- lambdasContra1
                        }else{

                                lambdasContra <- do.call(c, lambdas_list[[i]])
                                lambdasContra[str_detect(names(lambdasContra), indep_var)] <- lambdasContra1[str_detect(names(lambdasContra), indep_var)]
                                lambdasContra_matrix[, i] <- lambdasContra
                        }
                }

                rownames(lambdasContra_matrix) <- names(lambda)
                lambdasContra_list[[indep_var]] <- lambdasContra_matrix
        }


        list_counterfactuals <- NULL
        for(indep_var in indep_vars){

                Xb <- X%*%matrix_beta
                Xl <- X%*%matrix_lambda

                Xb_contra <- X%*%betasContra_list[[indep_var]]
                Xl_contra <- X%*%lambdasContra_list[[indep_var]]

                p        = matrix_weight
                p_contra = pContra_list[[indep_var]]

                mu          <- colSums(Xb        *p)
                mu_contra_m <- colSums(Xb_contra *p)
                mu_contra_p <- colSums(Xb        *p_contra)

                between          <- colSums(p*(t(t(Xb) - mu)^2))
                between_contra_m <- colSums(p*(t(t(Xb_contra) - mu_contra_m)^2))
                between_contra_p <- colSums(p_contra*(t(t(Xb) - mu_contra_p)^2))

                within_contra_p       <- colSums(p_contra*exp(Xl))
                within_contra_lambda  <- colSums(p*exp(Xl_contra))

                list_counterfactuals[[indep_var]] <- tibble(ID           = names(between),
                                                            between_beta = between_contra_m,
                                                            between_p    = between_contra_p,

                                                            within_lambda  = within_contra_lambda,
                                                            within_p       = within_contra_p,

                                                            mean_ef = variance_decomp$varLog - (between_beta + variance_decomp$within),
                                                            var_ef  = variance_decomp$varLog - (variance_decomp$between + within_lambda),
                                                            comp_ef = variance_decomp$varLog - (between_p + within_p)
                )
        }


        list(variance_decomp, list_counterfactuals)
}



#' @export
varDecomp.ordProbitModel <- function(model){
        matrix_coef <- coefficients(model)

        list_X <- model.matrix(model)
        list_w <- map(.x = names(model), function(x) model[[x]]$weights)

        formula   = model[[1]]$formula
        indep_vars <- all.vars(formula)

        k=1

        matrix_beta   = matrix_coef$beta_matrix
        matrix_lambda = matrix_coef$lambda_matrix

        lambdas_list = list()
        betas_list = list()

        variance_decomp <- NULL
        for(k in 1:length(list_X)){
                print(k)

                w = list_w[[k]]
                X = list_X[[k]]

                lambda <- matrix_lambda[, k]
                beta   <- matrix_beta[, k]

                p_i = w/sum(w)
                s2_i <- exp(X%*%lambda)

                mu_i <- as.numeric(X%*%beta)
                mu   <- sum(p_i*mu_i)

                between <- sum( p_i*((mu_i - mu)^2) )
                within  <- sum( p_i*s2_i )

                varLog  <- between + within
                between_share <- between/varLog
                within_share  <- within/varLog

                variance_decomp_k <- tibble(ID = names(list_X)[k],
                                            varLog, between, within,
                                            between_share, within_share)
                variance_decomp <- bind_rows(variance_decomp, variance_decomp_k)

                lambdas_list[[names(list_X)[k]]]$`(Intercept)` <- lambda[1]
                for(indep_var in indep_vars){
                        lambdas_list[[names(list_X)[k]]][[indep_var]] <- lambda[str_detect(names(lambda), indep_var)]
                }

                betas_list[[names(list_X)[k]]]$`"(Intercept)"` <- beta[1]
                for(indep_var in indep_vars){
                        betas_list[[names(list_X)[k]]][[indep_var]] <- beta[str_detect(names(beta), indep_var)]
                }
        }


        list_vars = model[[1]]$independentVariablesValues
        data_X <- do.call(expand.grid, list_vars)
        data_X <- data_X %>% mutate_all(as.character)

        X     <- model.matrix(formula, data_X)
        X     <- X[ , colnames(list_X[[1]])]
        IDs_X <- apply(X, 1, function(x) paste(x, collapse = "_"))

        data_weight <- tibble(IDs = IDs_X)
        for(i in 1:length(list_X)){
                data_weight_tmp <- tibble(IDs = apply(list_X[[i]], 1, function(x) paste(x, collapse = "_")),
                                          w   = model[[i]]$weights)

                data_weight_tmp <- data_weight_tmp %>%
                        group_by(IDs) %>%
                        summarise(p  = sum(w)) %>%
                        ungroup() %>%
                        mutate(p = p/sum(p))

                data_weight <- left_join(data_weight,
                                         data_weight_tmp)

                names(data_weight)[i+1] <- names(list_X[i])
        }

        matrix_weight <- as.matrix(data_weight[, -1])
        matrix_weight[is.na(matrix_weight)] <- 0

        #var = names(data_X)[1]
        i = 1

        pContra_list <- list()
        for(var in names(data_X)){
                data_p = data_X[var]

                for(i in 1:ncol(matrix_weight)){
                        ID       <- colnames(matrix_weight)[i]
                        tab      <- colSums(X*matrix_weight[,i])
                        freq_var <- tab[str_detect(names(tab), var)]

                        data_p_var <- tibble(unique(data_X[[var]]), c(1 - sum(freq_var), freq_var)) %>%
                                setNames(c(var, ID))

                        data_p <- left_join(data_p, data_p_var)
                }

                matrix_p <- data_p[-1] %>% as.matrix()
                matrix_p <- (1/(matrix_p/matrix_p[,1]))*matrix_weight # Western & Bloome, p. 313

                pContra_list[[var]] <- matrix_p
        }

        betasContra_list <- NULL
        for(indep_var in indep_vars){
                betasContra_matrix <- matrix(NA,
                                             ncol = length(betas_list),
                                             nrow = length(beta))

                for(i in 1:length(betas_list)){
                        if(i == 1){
                                betasContra1 <- do.call(c, betas_list[[i]])
                                betasContra_matrix[, 1] <- betasContra1
                        }else{

                                betasContra <- do.call(c, betas_list[[i]])
                                betasContra[str_detect(names(betasContra), indep_var)] <- betasContra1[str_detect(names(betasContra), indep_var)]
                                betasContra_matrix[, i] <- betasContra
                        }
                }

                rownames(betasContra_matrix) <- names(beta)
                betasContra_list[[indep_var]] <- betasContra_matrix
        }


        lambdasContra_list <- NULL
        for(indep_var in indep_vars){
                lambdasContra_matrix <- matrix(NA,
                                               ncol = length(lambdas_list),
                                               nrow = length(lambda))

                for(i in 1:length(lambdas_list)){
                        if(i == 1){
                                lambdasContra1 <- do.call(c, lambdas_list[[i]])
                                lambdasContra_matrix[, 1] <- lambdasContra1
                        }else{

                                lambdasContra <- do.call(c, lambdas_list[[i]])
                                lambdasContra[str_detect(names(lambdasContra), indep_var)] <- lambdasContra1[str_detect(names(lambdasContra), indep_var)]
                                lambdasContra_matrix[, i] <- lambdasContra
                        }
                }

                rownames(lambdasContra_matrix) <- names(lambda)
                lambdasContra_list[[indep_var]] <- lambdasContra_matrix
        }


        list_counterfactuals <- NULL
        for(indep_var in indep_vars){

                Xb <- X%*%matrix_beta
                Xl <- X%*%matrix_lambda

                Xb_contra <- X%*%betasContra_list[[indep_var]]
                Xl_contra <- X%*%lambdasContra_list[[indep_var]]

                p        = matrix_weight
                p_contra = pContra_list[[indep_var]]

                mu          <- colSums(Xb        *p)
                mu_contra_m <- colSums(Xb_contra *p)
                mu_contra_p <- colSums(Xb        *p_contra)

                between          <- colSums(p*(t(t(Xb) - mu)^2))
                between_contra_m <- colSums(p*(t(t(Xb_contra) - mu_contra_m)^2))
                between_contra_p <- colSums(p_contra*(t(t(Xb) - mu_contra_p)^2))

                within_contra_p       <- colSums(p_contra*exp(Xl))
                within_contra_lambda  <- colSums(p*exp(Xl_contra))

                list_counterfactuals[[indep_var]] <- tibble(ID           = names(between),
                                                            between_beta = between_contra_m,
                                                            between_p    = between_contra_p,

                                                            within_lambda  = within_contra_lambda,
                                                            within_p       = within_contra_p,

                                                            mean_ef = variance_decomp$varLog - (between_beta + variance_decomp$within),
                                                            var_ef  = variance_decomp$varLog - (variance_decomp$between + within_lambda),
                                                            comp_ef = variance_decomp$varLog - (between_p + within_p)
                )
        }


        list(variance_decomp, list_counterfactuals)
}



