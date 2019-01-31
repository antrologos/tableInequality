#' @export

calc_gini_RPME <- function(data_pnad, groups = NULL){

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
                        gather(key = type,
                               value = "gini",
                               midArit_topArit:midGeom_topHarm)

        }else{
                gini_result <- gini_result %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        gather(key = type,
                               value = "gini",
                               midArit_topArit:midGeom_topHarm)
        }

        gini_result
}


