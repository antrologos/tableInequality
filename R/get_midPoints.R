#' @export

get_midPoints <- function(data_pnad, groups = NULL){

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

                       midArit_topArit   = ifelse(is.na(arit_midPoints), aritMean, arit_midPoints),
                       midArit_topMedian = ifelse(is.na(arit_midPoints), median,   arit_midPoints),
                       midArit_topGeom   = ifelse(is.na(arit_midPoints), geoMean,  arit_midPoints),
                       midArit_topHarm   = ifelse(is.na(arit_midPoints), HarMean,  arit_midPoints),

                       midGeom_topArit   = ifelse(is.na(geom_midPoints), aritMean, geom_midPoints),
                       midGeom_topMedian = ifelse(is.na(geom_midPoints), median,   geom_midPoints),
                       midGeom_topGeom   = ifelse(is.na(geom_midPoints), geoMean,  geom_midPoints),
                       midGeom_topHarm   = ifelse(is.na(geom_midPoints), HarMean,  geom_midPoints)
                ) %>%
                dplyr::select(ID, n, midArit_topArit:midGeom_topHarm)

        pnads_midpoints <- bind_cols(pnads_midpoints,
                                     data_pnad %>% dplyr::select(min_faixa, max_faixa))

        if(is.null(groups)){
                pnads_midpoints
        }else{
                pnads_midpoints %>%
                        separate(col = ID, into = groups, sep = "_")
        }
}


