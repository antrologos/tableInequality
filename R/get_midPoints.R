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

        data_split = split(data_pnad, data_pnad$ID)
        IDs = names(data_split)

        #data_i = data_split[[1]]
        get_alpha = function(data_i){

                calc_alpha <- try(
                        getMids(ID = data_i$ID,
                                hb = data_i$n,
                                lb = data_i$min_faixa,
                                ub = data_i$max_faixa,
                                alpha_bound = numeric(0))$alpha,
                        silent = T)

                alpha = ifelse("try-error" %in% class(calc_alpha),
                               NA,
                               as.numeric(calc_alpha))

                tibble(ID = unique(data_i$ID), alpha)
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        data_alpha <- future_map_dfr(.x = data_split, .f = get_alpha, .progress = T)


        #for(i in 1:length(data_split)){
        #        print(i)
        #        get_alpha(data_split[[i]])
        #}

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
                                     alpha           = data_alpha$alpha,
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
                pnads_midpoints$ID = NULL
                pnads_midpoints
        }else{
                if(length(groups) == 1){
                        pnads_midpoints[[groups]] = pnads_midpoints[["ID"]]

                        if(groups != "ID"){
                                pnads_midpoints$ID = NULL
                        }
                        pnads_midpoints
                }else{
                        pnads_midpoints %>%
                                separate(col = ID, into = groups, sep = "_")
                }
        }
}


