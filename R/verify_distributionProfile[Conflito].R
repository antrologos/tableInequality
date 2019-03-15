#' @export

verify_distributionProfile <- function(data_pnad, groups = NULL){

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
                mutate(perc85_or_more_1st_bracket = `1` > .85,
                       perc85_or_more_1st_e_2nd_brackets    = `1` + `2` > .85)


        verification_step2 <-
                verification_step1 %>%
                mutate(distributionProfile = "c. well distributed",
                       distributionProfile = ifelse(perc85_or_more_1st_e_2nd_brackets, "b. 90%+ - 1st & 2nd bracket", distributionProfile),
                       distributionProfile = ifelse(perc85_or_more_1st_bracket,        "a. 90%+ - 1st bracket", distributionProfile))

        verification_step3 <- verification_step2 %>%
                dplyr::select(ID, distributionProfile)

        if(is.null(groups)){
                verification_step3

        }else{
                verification_step3 %>%
                        separate(col = ID, into = groups, sep = "_")

        }
}
