#' @export

calc_theil_RPME <- function(data_pnad, groups = NULL){

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

        data_pnad_groups <- get_midPoints(data_pnad, groups = "ID")

        data_split <- split(data_pnad_groups, f = data_pnad_groups$ID)

        #data_i = data_split[[1]]
        calc_theil <- function(data_i){

                x = data_i$midArit_topArit
                w = data_i$n

                x_bar = sum(x*w)/sum(w)

                theil = (1/sum(w))*sum(w*((x/x_bar)*log(x/x_bar)))

                theil
        }

        theil_result <- tableInequality:::future_map_parallel(.x = data_split,
                                                              .f = ~calc_theil(.x),
                                                              .progress = T) %>%
                tibble(ID = names(.),  theil = unlist(.))

        if(is.null(groups)){
                theil_result <- theil_result %>%
                        dplyr::select(theil)
        }else{
                theil_result <- theil_result %>%
                        dplyr::select(ID, theil) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        theil_result

}

