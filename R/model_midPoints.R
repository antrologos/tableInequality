#' @export
model_midPoints <- function(formula,
                            data_pnad,
                            weight = NULL,
                            groups = NULL){

        dep   = all.vars(formula[[2]])
        indep = all.vars(formula[[3]])

        if(any(indep %in% groups) | any(groups %in% indep)){
                stop("Group variables cannot be used as independent variables in the model")
        }

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange_at(ID, min_faixa)
        }else{
                data_pnad <-
                        data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by_at(c("ID",dep,indep)) %>%
                        summarise(min_faixa = min(min_faixa),
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[1]]

        regression = function(data_i, formula, weight){

                if(is.null(weight)){
                        vars = all.vars(formula)
                }else{
                        vars = c(all.vars(formula), weight)
                }

                ID_i <- data_i$ID %>% unique()

                data_i$case_number = 1:nrow(data_i)

                data_ii = data_i %>%
                        dplyr::select(vars, case_number) %>%
                        filter(complete.cases(.))

                if(is.null(weight)){
                        data_ii$w = rep(1, length(y))
                }else{
                        data_ii$w = data_ii[[weight]]
                }

                data_ii <- data_ii[data_ii$w > 0, ]

                X = model.matrix(formula, data = data_ii)
                y = with(data_ii, eval(formula[[2]]))

                w      = data_ii$w
                sqrt_w = sqrt(w)

                beta = try(
                        crossprod(solve(crossprod(sqrt_w*X)), crossprod(sqrt_w*X, sqrt_w*y)) %>% as.numeric(),
                        silent = T)

                if("try-error" %in% class(beta)){
                        beta = matrix(rep(as.numeric(NA), ncol(X)), nrow = 1)
                        colnames(beta) = colnames(X)
                }

                N     = sum(w)
                k     = ncol(X)
                y_hat = X %*% beta
                res   = as.numeric(y - y_hat)
                var_e = sum(w*(res^2))/(N - k)
                sd_e  = sqrt(var_e)

                var_beta_he <- solve(crossprod(sqrt_w*X)) %*% t(sqrt_w*X) %*% (w*(res^2)*sqrt_w*X) %*% solve(crossprod(sqrt_w*X))
                sd_beta_he  <- sqrt(diag(var_beta_he))

                t_statistic = beta/sd_beta_he
                p_value = 2*(1 - pt(abs(t_statistic),df = N-k))

                mean_y = sum(w*y)/sum(w)
                SSE = sum(w*((y - y_hat)^2))
                SST = sum(w*((y - mean_y)^2))

                r2 = 1 - SSE/SST

                parameters <- tibble(parameter    = colnames(X),
                                     coefficients = beta,
                                     sd           = sd_beta_he,
                                     t_statistic  = t_statistic,
                                     p_value      = p_value)

                list(ID            = tibble(ID = ID_i) %>% separate(col = ID, into = groups, sep = "_"),
                     parameters    = parameters,
                     dependent_var = y,
                     model_matrix  = X,
                     residuals     = res,
                     weights       = w,
                     rmse          = sd_e,
                     r2            = r2,
                     N             = N)

                #bind_cols(tibble(ID = ID_i), as_tibble(beta))

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        #betas = NULL
        #for(i in 1:length(data_split)){
        #        print(i)

        #        betas = bind_rows(betas,
        #                           regression(data_i = data_split[[i]],
        #                                      formula = formula,
        #                                      weight = weight))
        #}

        betas_result <- future_map(.x = data_split,
                                .f = regression,
                                formula = formula,
                                weight = weight,
                                .progress = T)

        class(betas_result) <- "midpointModel"

        betas_result
}

#' @export
residuals.midpointModel <- function(midpointModel){
        residuals_list = purrr::map(midpointModel, function(x) x$residuals)
        residuals_list
}


#' @export
model.matrix.midpointModel <- function(midpointModel){
        model_matrix_list = purrr::map(midpointModel, function(x) x$model_matrix)
        model_matrix_list
}

#' @export
model.matrix.midpointModel <- function(midpointModel, weights = F){
        model_matrix_list = purrr::map(midpointModel, function(x) x$model_matrix)

        if(weights == T){
                weights_list = purrr::map(midpointModel, function(x) x$weights)

                for(i in 1:length(model_matrix_list)){
                        model_matrix_list[[i]] <- cbind(model_matrix_list[[i]], weights_list[[i]])
                }

        }

        model_matrix_list
}

#' @export
coefficients.midpointModel <- function(midpointModel){
        coefficients_list = purrr::map(midpointModel, function(x) x$parameters$coefficients)
        coefficients_matrix <- do.call(cbind, coefficients_list)

        var_names <- midpointModel[[1]]$parameters[[1]] %>% as.character()

        rownames(coefficients_matrix) <- var_names
}


#' @export
summary.midpointModel <- function(midpointModel, digits = 3){

        beta_matrix <- purrr::map(midpointModel, function(x) as_tibble(matrix(x$parameters$coefficients, nrow = 1)) %>% setNames(x$parameters$parameter))
        sd_matrix   <- purrr::map(midpointModel, function(x) as_tibble(matrix(x$parameters$sd,           nrow = 1)) %>% setNames(x$parameters$parameter))
        p_matrix    <- purrr::map(midpointModel, function(x) as_tibble(matrix(x$parameters$p_value,      nrow = 1)) %>% setNames(x$parameters$parameter))

        beta_matrix <- do.call(bind_rows, beta_matrix) %>% as.matrix() %>% t()
        sd_matrix   <- do.call(bind_rows, sd_matrix)   %>% as.matrix() %>% t()
        p_matrix    <- do.call(bind_rows, p_matrix)    %>% as.matrix() %>% t()

        var_names   <- rownames(beta_matrix)

        rmse        <- purrr::map_dbl(midpointModel, function(x) x$rmse)
        r2          <- purrr::map_dbl(midpointModel, function(x) x$r2)
        N           <- purrr::map_dbl(midpointModel, function(x) x$N)

        beta_matrix_char <- format(beta_matrix, digits = digits)
        sd_matrix_char   <- format(sd_matrix,   digits = digits)

        beta_matrix_char <- as_tibble(beta_matrix_char) %>%
                mutate_all(function(x) x %>% str_replace_all(pattern = "NA", "") %>% str_trim()) %>%
                as.matrix()

        sd_matrix_char <- as_tibble(sd_matrix_char) %>%
                mutate_all(function(x) x %>% str_replace_all(pattern = "NA", "") %>% str_trim()) %>%
                as.matrix()

        sd_matrix_char   <- paste0("(",format(sd_matrix_char, digits = digits),")")
        dim(sd_matrix_char) <- dim(beta_matrix_char)
        sd_matrix_char[is.na(sd_matrix)] <- ""

        star1 <- p_matrix < .1
        star2 <- p_matrix < .05
        star3 <- p_matrix < .01

        star1[is.na(star1)] <- FALSE
        star2[is.na(star2)] <- FALSE
        star3[is.na(star3)] <- FALSE

        beta_matrix_char[star1] <- paste0(beta_matrix_char[star1],"*")
        beta_matrix_char[star2] <- paste0(beta_matrix_char[star2],"*")
        beta_matrix_char[star3] <- paste0(beta_matrix_char[star3],"*")

        results = matrix("",
                         ncol = ncol(beta_matrix_char),
                         nrow = 2*nrow(beta_matrix_char))

        uneven_lines = seq(1, 2*nrow(beta_matrix_char), 2)
        even_lines   = uneven_lines + 1

        results[uneven_lines,] <- beta_matrix_char
        results[even_lines,]   <- sd_matrix_char

        results <- rbind(results, "",
                         format(rmse, digits = digits),
                         format(r2,   digits = digits),
                         format(N,    digits = 0))

        rownames(results) <- c(rep("", 2*length(var_names)), "", "RMSE", "R2", "N")
        rownames(results)[uneven_lines] <- var_names

        results
}

#' @export
print.midpointModel <- function(midpointModel){
        print(summary(midpointModel))
}

