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

                data_ii = data_i %>%
                        dplyr::select(vars) %>%
                        filter(complete.cases(.))

                if(is.null(weight)){
                        w = rep(1, length(y))
                }else{
                        w = data_ii[[weight]]
                }

                data_ii <- data_ii[w > 0, ]
                w = data_ii[[weight]]

                X = model.matrix(formula, data = data_ii)
                y = with(data_ii, eval(formula[[2]]))

                sqrt_w = sqrt(w)

                beta = try(
                        crossprod(solve(crossprod(sqrt_w*X)), crossprod(sqrt_w*X, sqrt_w*y)) %>% t(),
                        silent = T)

                if("try-error" %in% class(beta)){
                        beta = matrix(rep(as.numeric(NA), ncol(X)), nrow = 1)
                        colnames(beta) = colnames(X)
                }

                bind_cols(tibble(ID = ID_i), as_tibble(beta))

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

        betas_result <- future_map_dfr(.x = data_split,
                                .f = regression,
                                formula = formula,
                                weight  = weight,
                                .progress = T)

        if(is.null(groups)){
                betas_result <- betas_result %>%
                        dplyr::select(-ID)
        }else{
                betas_result <- betas_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        betas_result
}

