#' @export

calc_mean_stepbins <- function(data_pnad, groups = NULL){

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

        mean_stepbins = function(data){

                fit <- with(data, {
                        limites  = c(min_faixa[1],max_faixa)
                        contagem = c(0, n)
                        stepbins(bEdges = limites, bCounts = contagem)
                })

                # grand mean
                mean = integrate(f = function(z){z*fit$stepPDF(z)},
                                 lower = 0,
                                 upper =  fit$E,
                                 subdivisions = 2000,
                                 stop.on.error = F)$value

                mean
        }

        mean_result <- map(.x = data_split, .f = mean_stepbins) %>%
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
