#' @export

#data_pnad <- pnads_1968_1973
#groups = c("ano", "trimestre", "regiao", "educacao")

calc_mean_RPME <- function(data_pnad, groups = NULL){

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

        pnads_midpoints = get_midPoints(data_pnad, groups = "ID")

        pnads_midpoints_byGroups <- split(pnads_midpoints, pnads_midpoints$ID)

        var_names <- pnads_midpoints_byGroups[[1]] %>%
                dplyr::select(midArit_topArit:midGeom_topHarm) %>%
                names()

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        mean_result <- future_map_dfr(1:length(pnads_midpoints_byGroups), .f = function(i){

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
                        gather(key = type,
                               value = "mean",
                               midArit_topArit:midGeom_topHarm)

        }else{
                mean_result <- mean_result %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        gather(key = type,
                               value = "mean",
                               midArit_topArit:midGeom_topHarm)
        }

        mean_result
}

