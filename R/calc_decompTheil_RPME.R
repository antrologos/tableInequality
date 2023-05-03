#' @export

calc_decompTheil_RPME <- function(data_pnad, groups_theil, groups = NULL){

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange_at(c("ID", groups_theil, "min_faixa"))
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by_at(c("ID", groups_theil, "min_faixa")) %>%
                        summarise(
                                max_faixa = max(max_faixa),
                                n         = sum(n)) %>%
                        ungroup() %>%
                        arrange_at(c("ID", groups_theil, "min_faixa"))
        }

        data_pnad_groups <- get_midPoints(data_pnad, groups = c("ID", groups_theil))
        data_split       <- split(data_pnad_groups, f = data_pnad_groups$ID)

        calc_theil <- function(data_i){

                x = data_i$midArit_topArit
                w = data_i$n

                x_bar = sum(x*w)/sum(w)

                theil = (1/sum(w))*sum(w*((x/x_bar)*log(x/x_bar)))

                tibble(mean = x_bar, theil)
        }

        #data_i = data_split[[1]]
        calc_decompTheil <- function(data_i){

                data_i_split <- split(data_i, data_i[[groups_theil]])

                group_mean_theil  <- map(.x = data_i_split,
                                         .f = calc_theil) %>%
                        rbindlist(idcol = TRUE)

                group_sizes    <- map(data_i_split, function(x) tibble(n = sum(x$n))) %>%
                        rbindlist(idcol = TRUE)

                group_results <- left_join(group_mean_theil, group_sizes) %>%
                        mutate(grand_mean = sum(n*mean)/sum(n),
                               share      = (n/sum(n))*(mean/grand_mean))

                between <- with(group_results, sum(share*log(mean/grand_mean)) )
                within  <- with(group_results, sum(theil*share))

                tibble(theil = between + within, between, within)

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multisession)
        }

        #calc_decompTheil(data_split[[1]])
        decompTheil_result <- future_map_dfr(.x = data_split,
                                             .f = calc_decompTheil,
                                             .progress = T)

        decompTheil_result <- bind_cols(tibble(ID = names(data_split)), decompTheil_result)

        if(is.null(groups)){
                decompTheil_result <- decompTheil_result %>%
                        dplyr::select(theil, between, within)
        }else{
                decompTheil_result <- decompTheil_result %>%
                        dplyr::select(ID, theil, between, within) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        decompTheil_result

}

