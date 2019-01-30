#' @export


verify_distributionProfile <- function(data_pnad, groups){

        data_pnad <- data_pnad %>%
                unite(col = ID, groups) %>%
                group_by(ID, min_faixa) %>%
                summarise(
                        max_faixa = max(max_faixa),
                        n         = sum(n)) %>%
                ungroup() %>%
                arrange(ID, min_faixa)

        verification_step1 <-
                data_pnad %>%
                arrange(ID, min_faixa) %>%
                group_by(ID) %>%
                mutate(p     = n/sum(n),
                       rank  = rank(min_faixa)) %>%
                ungroup() %>%
                dplyr::select(ID, p, rank) %>%
                spread(key = rank, value = p, fill = 0) %>%
                mutate(perc90_or_more_1st_bracket = `1` > .9,
                       perc80_or_more_1st_bracket = `1` > .8,
                       perc70_or_more_1st_bracket = `1` > .7,
                       perc60_or_more_1st_bracket = `1` > .6,
                       perc50_or_more_1st_bracket = `1` > .5,
                       perc90_or_more_1st_e_2nd_brackets    = `1` + `2` > .9,
                       perc80_or_more_1st_e_2nd_brackets    = `1` + `2` > .8,
                       perc70_or_more_1st_e_2nd_brackets    = `1` + `2` > .7,
                       perc60_or_more_1st_e_2nd_brackets    = `1` + `2` > .6,
                       perc50_or_more_1st_e_2nd_brackets    = `1` + `2` > .5,
                       perc90_or_more_1st_2nd_e_3rd_brackets = `1` + `2` + `3` > .9,
                       perc80_or_more_1st_2nd_e_3rd_brackets = `1` + `2` + `3` > .8,
                       perc70_or_more_1st_2nd_e_3rd_brackets = `1` + `2` + `3` > .7,
                       perc60_or_more_1st_2nd_e_3rd_brackets = `1` + `2` + `3` > .6,
                       perc50_or_more_1st_2nd_e_3rd_brackets = `1` + `2` + `3` > .5)


        verification_step2 <-
                verification_step1 %>%
                mutate(distributionProfile = "q. well distributed",
                       distributionProfile = ifelse(perc50_or_more_1st_2nd_e_3rd_brackets, "p. 50%+ - 3rd bracket", distributionProfile),
                       distributionProfile = ifelse(perc60_or_more_1st_2nd_e_3rd_brackets, "o. 60%+ - 3rd bracket", distributionProfile),
                       distributionProfile = ifelse(perc70_or_more_1st_2nd_e_3rd_brackets, "n. 70%+ - 3rd bracket", distributionProfile),
                       distributionProfile = ifelse(perc80_or_more_1st_2nd_e_3rd_brackets, "m. 80%+ - 3rd bracket", distributionProfile),
                       distributionProfile = ifelse(perc90_or_more_1st_2nd_e_3rd_brackets, "l. 90%+ - 3rd bracket", distributionProfile),

                       distributionProfile = ifelse(perc50_or_more_1st_e_2nd_brackets,    "k. 50%+ - 2nd bracket", distributionProfile),
                       distributionProfile = ifelse(perc60_or_more_1st_e_2nd_brackets,    "j. 60%+ - 2nd bracket", distributionProfile),
                       distributionProfile = ifelse(perc70_or_more_1st_e_2nd_brackets,    "h. 70%+ - 2nd bracket", distributionProfile),
                       distributionProfile = ifelse(perc80_or_more_1st_e_2nd_brackets,    "g. 80%+ - 2nd bracket", distributionProfile),
                       distributionProfile = ifelse(perc90_or_more_1st_e_2nd_brackets,    "f. 90%+ - 2nd bracket", distributionProfile),

                       distributionProfile = ifelse(perc50_or_more_1st_bracket,  "e. 50%+ - 1st bracket", distributionProfile),
                       distributionProfile = ifelse(perc60_or_more_1st_bracket,  "d. 60%+ - 1st bracket", distributionProfile),
                       distributionProfile = ifelse(perc70_or_more_1st_bracket,  "c. 70%+ - 1st bracket", distributionProfile),
                       distributionProfile = ifelse(perc80_or_more_1st_bracket,  "b. 80%+ - 1st bracket", distributionProfile),
                       distributionProfile = ifelse(perc90_or_more_1st_bracket,  "a. 90%+ - 1st bracket", distributionProfile))

        verification_step3 <- verification_step2 %>%
                dplyr::select(ID, distributionProfile)


        if(is.null(groups)){
                verification_step3

        }else{
                verification_step3 %>%
                        separate(col = ID, into = groups, sep = "_")

        }
}


