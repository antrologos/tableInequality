#' @export

calc_quantile_rsubbins <- function(p, data_pnad, groups = NULL){

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

        #data_i <- data_split[[7]]
        quantile_rsubbins = function(p, data_i = data_i){

                limites  <- c(data_i$min_faixa[1], data_i$max_faixa)
                contagem <- c(0, data_i$n)
                fit      <- rsubbins(bEdges = limites, bCounts = contagem)

                quantiles_i = tableInequality:::quantile(p , CDF_func = fit$rsubCDF, max_x = fit$E)

                as_tibble(matrix(quantiles_i, nrow = 1)) %>% setNames(p)
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        quantile_result <- future_map_dfr(.x = data_split,
                                          .f = quantile_rsubbins, p = p,
                                          .progress = T) %>%
                bind_cols(ID = names(data_split), .)

        if(is.null(groups)){
                quantile_result <- quantile_result %>%
                        dplyr::select(-ID)
        }else{
                quantile_result <- quantile_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        quantile_result
}

