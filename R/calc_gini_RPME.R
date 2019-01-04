#' @export

calc_gini_RPME <- function(data_pnad, groups = NULL){

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

        data_pnad <- data_pnad %>% mutate(ID = as.character(ID))

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

        gini_result <- map_dfr(1:length(pnads_midpoints_byGroups), .f = function(i){

                nome <- names(pnads_midpoints_byGroups)[i]
                data <- pnads_midpoints_byGroups[[i]] %>%
                        filter(n > 0)

                if(nrow(data) > 1){
                        map(var_names, .f = function(var) calcSGini(data[[var]], data[["n"]])$ineq$index) %>%
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
                gini_result <- gini_result %>%
                        gather(key = `Tipo de Média`,
                               value = "gini",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)

        }else{
                gini_result <- gini_result %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        gather(key = `Tipo de Média`,
                               value = "gini",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)
        }

        gini_result
}
