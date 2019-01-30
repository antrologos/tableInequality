#' @export

calc_gini_splinebins <- function(data_pnad, groups = NULL){

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

        gini_splinebins = function(data_i){

                limites  <- c(data_i$min_faixa[1],data_i$max_faixa)
                contagem <- c(0, data_i$n)

                fit <- splinebins(bEdges = limites, bCounts = contagem)

                gini_numerical_integration(PDF_func = fit$splinePDF,
                                           CDF_func = fit$splineCDF,
                                           max_x    = fit$E)
        }

        gini_result <- future_map_parallel(.x = data_split, .f = ~gini_splinebins(.x), .progress = T) %>%
                tibble(ID = names(.),  gini = unlist(.))

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
