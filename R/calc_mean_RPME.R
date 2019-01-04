#' @export

calc_mean_RPME <- function(data_pnad, groups = NULL){

        data_pnad <- data_pnad %>%
                unite(col = ID, groups) %>%
                group_by(ID, faixas_renda) %>%
                summarise(min_faixa = min(min_faixa),
                          max_faixa = max(max_faixa),
                          n         = sum(n)) %>%
                ungroup() %>%
                arrange(ID, min_faixa)


        bin_mids_unbound <- with(data_pnad, {
                getMids(ID = ID,
                        hb = n,
                        lb = min_faixa,
                        ub = max_faixa,
                        alpha_bound = numeric(0))
        })

        colMax <- function(...){
                apply(cbind(...), 1, max, na.rm=T)
        }

        midPoints_highestBin <- data_pnad %>%
                ungroup() %>%
                filter(is.na(max_faixa)) %>%
                dplyr::select(ID, min_faixa) %>%
                rename(lowerbound = min_faixa) %>%
                left_join(y = tibble(ID              = data_pnad$ID %>% unique() %>% as.character(), #bin_mids_unbound$mids$ID %>% unique() %>% as.character(),
                                     alpha           = bin_mids_unbound$alpha,
                                     alpha_bounded_1 = colMax(1,alpha),
                                     alpha_bounded_2 = colMax(2,alpha)),
                          by = "ID") %>%
                mutate(aritMean = lowerbound * (alpha_bounded_2/(alpha_bounded_2 - 1)),
                       median   = lowerbound*(2^(1/alpha_bounded_1)),
                       geoMean  = lowerbound*(exp(1/alpha_bounded_1)),
                       HarMean  = lowerbound*(1 + (1/alpha_bounded_1)))

        pnads_midpoints <- data_pnad %>%
                ungroup() %>%
                left_join(y  = midPoints_highestBin %>% dplyr::select(ID, aritMean:HarMean),
                          by = "ID") %>%
                mutate(arit_midPoints = (min_faixa + max_faixa)/2,
                       geom_midPoints = (min_faixa*max_faixa)^(1/2),

                       arit_midPoints_aritMean = ifelse(is.na(arit_midPoints), aritMean, arit_midPoints),
                       arit_midPoints_median   = ifelse(is.na(arit_midPoints), median,   arit_midPoints),
                       arit_midPoints_geoMean  = ifelse(is.na(arit_midPoints), geoMean,  arit_midPoints),
                       arit_midPoints_HarMean  = ifelse(is.na(arit_midPoints), HarMean,  arit_midPoints),

                       geom_midPoints_aritMean = ifelse(is.na(geom_midPoints), aritMean, geom_midPoints),
                       geom_midPoints_median   = ifelse(is.na(geom_midPoints), median,   geom_midPoints),
                       geom_midPoints_geoMean  = ifelse(is.na(geom_midPoints), geoMean,  geom_midPoints),
                       geom_midPoints_HarMean  = ifelse(is.na(geom_midPoints), HarMean,  geom_midPoints)
                ) %>%
                dplyr::select(ID, n, arit_midPoints_aritMean:geom_midPoints_HarMean)

        pnads_midpoints_byGroups <- split(pnads_midpoints, pnads_midpoints$ID)

        var_names <- pnads_midpoints_byGroups[[1]] %>%
                dplyr::select(arit_midPoints_aritMean:geom_midPoints_HarMean) %>%
                names()

        mean_result <- map_df(1:length(pnads_midpoints_byGroups), .f = function(i){

                nome <- names(pnads_midpoints_byGroups)[i]
                data <- pnads_midpoints_byGroups[[i]] %>%
                        filter(n > 0)

                if(nrow(data) > 1){
                        map(var_names, .f = function(var){
                                sum(data[[var]]*(data$n/sum(data$n)))
                        }) %>%
                                setNames(var_names) %>%
                                as_tibble() %>%
                                mutate(ID = nome)
                }else{
                        matrix(as.numeric(NA), nrow = 1, ncol = 8) %>%
                                as_tibble() %>%
                                setNames(var_names)%>%
                                mutate(ID = nome)
                }
        })


        if(is.null(groups)){
                mean_result <- mean_result %>%
                        gather(key = `Tipo de Média`,
                               value = "mean",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)

        }else{
                mean_result <- mean_result %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        gather(key = `Tipo de Média`,
                               value = "mean",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)
        }

        mean_result
}
