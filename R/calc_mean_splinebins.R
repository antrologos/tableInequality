#' @export

calc_mean_splinebins <- function(data_pnad, groups = NULL){

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by(ID, faixas_renda) %>%
                        summarise(min_faixa = min(min_faixa),
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split <- split(data_pnad, f = data_pnad$ID)

        mean_splinebins = function(data_i){

                limites  <- c(data_i$min_faixa[1],data_i$max_faixa)
                contagem <- c(0, data_i$n)

                fit <- splinebins(bEdges = limites, bCounts = contagem)

                # mean by integration (quadrature - Gauss-Legendre)
                grid_mean = mvQuad::createNIGrid(dim=1, type="GLe", level=75)
                mvQuad::rescale(grid_mean, domain = matrix(c(0, fit$E), ncol=2))

                mean = mvQuad::quadrature(f = function(y) y*splinePDF(y), grid = grid_mean)

                mean
        }

        mean_result <- future_map(.x = data_split, .f = mean_splinebins) %>%
                tibble(ID = names(.),  mean = unlist(.))

        if(is.null(groups)){
                mean_result <- mean_result %>%
                        dplyr::select(mean)
        }else{
                mean_result <- mean_result %>%
                        dplyr::select(ID, mean) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        mean_result

}
