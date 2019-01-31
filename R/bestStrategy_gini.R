#' @export

bestStrategy_gini <- function(data_pnad, groups = NULL){

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


        ID_list <- data_pnad$ID %>% unique()
        distributionProfile = verify_distributionProfile(data_pnad, "ID")

        best_strategy_list <- system.file("extdata",
                                     "bestStrategy_for_distributionProfile.csv",
                                     package = "tableInequality") %>%
                read.csv(stringsAsFactors = F)

        bestStrategy <- left_join(distributionProfile,
                                   best_strategy_list,
                                   by = "distributionProfile")
        bestStrategy <-
                bestStrategy %>%
                rename(bestStrategy = tipo) %>%
                dplyr::select(ID, bestStrategy)

        data_pnad <- left_join(x = data_pnad,
                               y = bestStrategy)

        bestStrategy_split = split(x = data_pnad, data_pnad$bestStrategy)

        results <- NULL
        #bestStrategy_name = names(bestStrategy_split)[[1]]
        for(bestStrategy_name in names(bestStrategy_split)){

                print(bestStrategy_name)

                data_split_strategy = bestStrategy_split[[bestStrategy_name]]

                if(str_detect(bestStrategy_name, "RPME: ")){

                        selected_type = str_remove(bestStrategy_name, "RPME: ")

                        result <- calc_gini_RPME(data_split_strategy,
                                                 groups = "ID") %>%
                                filter(type == selected_type)
                }

                if(bestStrategy_name == "Splinebins"){
                        result <- calc_gini_splinebins(data_split_strategy,
                                                       groups = "ID")
                }

                if(bestStrategy_name == "stepbins"){
                        result <- calc_gini_stepbins(data_split_strategy,
                                                     groups = "ID")
                }

                if(bestStrategy_name == "rsubbins"){
                        result <- calc_gini_rsubbins(data_split_strategy,
                                                     groups = "ID")
                }

                if(bestStrategy_name == "MGBE"){
                        result <- calc_gini_MGBE(data_split_strategy,
                                                 groups = "ID")
                }

                if(bestStrategy_name == "logNormal"){
                        result <- calc_gini_logNormal(data_split_strategy,
                                                      groups = "ID")
                }

                if(bestStrategy_name == "logNormalPareto P85"){
                        result <- calc_gini_logNormalPareto(data_split_strategy,
                                                            groups = "ID",
                                                            limite_distribuicoes = .85)
                }

                if(bestStrategy_name == "logNormalPareto P90"){
                        result <- calc_gini_logNormalPareto(data_split_strategy,
                                                            groups = "ID",
                                                            limite_distribuicoes = .9)
                }

                if(bestStrategy_name == "logNormalPareto P95"){
                        result <- calc_gini_logNormalPareto(data_split_strategy,
                                                            groups = "ID",
                                                            limite_distribuicoes = .95)
                }

                if(bestStrategy_name == "paretoLocal"){
                        result <- calc_gini_paretoLocalThreasholds(data_split_strategy,
                                                                   groups = "ID")
                }

                if(bestStrategy_name == "gpinter"){
                        result <- calc_gini_gpinter(data_split_strategy,
                                                    groups = "ID")
                }

                if(bestStrategy_name == "MCIB - RPME"){
                        result <- calc_gini_MCIB(data_split_strategy,
                                                 groups = "ID",
                                                 topBracket_method = "RPME")
                }

                if(bestStrategy_name == "MCIB - gpinter"){
                        result <- calc_gini_MCIB(data_split_strategy,
                                                 groups = "ID",
                                                 topBracket_method = "gpinter")
                }


                if(str_detect(bestStrategy_name, "Parabolic: ")){

                        selected_type = str_remove(bestStrategy_name, "Parabolic: ")

                        midpoints_data = get_midPoints(data_split_strategy, groups = "ID") %>%
                                gather(key = type, value = midpoint, midArit_topArit:midGeom_topHarm) %>%
                                filter(type == selected_type)

                        result <- calc_gini_parabolicInterp(midpoints_data,
                                                            groups = "ID")
                }

                result$strategy = bestStrategy_name


                results <- bind_rows(results, result)

        }


        final_results = left_join(tibble(ID = ID_list),
                                  results) %>%
                dplyr::select(ID, gini, strategy)

        if(is.null(groups)){
                final_results <- final_results %>%
                        dplyr::select(-ID)
        }else{
                final_results <- final_results %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        final_results
}


