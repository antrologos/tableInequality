#' @export

calc_gini_stepbins <- function(data_pnad, groups = NULL){

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

        gini_stepbins = function(data_i){

                limites  <- c(data_i$min_faixa[1],data_i$max_faixa)
                contagem <- c(0, data_i$n)

                fit <- suppressWarnings(stepbins(bEdges = limites, bCounts = contagem))

                tableInequality:::gini_numerical_integration(PDF_func = fit$stepPDF,
                                           CDF_func = fit$stepCDF,
                                           max_x    = fit$E)
        }

        ginis <- map(.x = data_split, .f = ~gini_stepbins(.x)) %>%
                unlist()

        gini_result <- tibble(ID = names(data_split),
                              gini = ginis)

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
