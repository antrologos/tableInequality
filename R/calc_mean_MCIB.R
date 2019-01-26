#' @export

calc_mean_MCIB <- function(data_pnad, groups = NULL, topBracket_method = c("gpinter", "RMPE")){

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


        data_split <- split(data_pnad, f = data_pnad$ID)

        topBracket_method_chosen = topBracket_method[1]

        #data_i = data_split[[1]]

        mean_MCIB = function(data_i, topBracket_method_chosen){

                ID_i = data_i$ID %>% unique()
                lower_i = data_i$min_faixa
                upper_i = data_i$max_faixa
                n_i     = data_i$n

                N = sum(n_i)

                lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
                upper_i = c(lower_i[-1], NA)

                slope_intercept <- tableInequality:::estimate_slope_and_intercept_MCIB(lower_i = lower_i,
                                                                                       upper_i = upper_i,
                                                                                       n_i = n_i)
                m = slope_intercept$m
                c = slope_intercept$c

                PDFpareto_lastBracket = tableInequality:::make_PDFpareto_lastBracket(data_i = data_i,
                                                                                     topBracket_method_chosen = topBracket_method_chosen,
                                                                                     known_groupMeans_checked = NULL,
                                                                                     m = m,
                                                                                     c = c)

                pareto_upper_bound = tableInequality:::get_pareto_upper_bound(data_i = data_i,
                                                                              topBracket_method_chosen = topBracket_method_chosen,
                                                                              known_groupMeans_checked = NULL,
                                                                              m = m,
                                                                              c = c)

                pdf_MCIB = function(y){

                        n_brackets = length(lower_i)
                        i = sapply(y, function(x){
                                which((x >= lower_i) & (x < upper_i))
                        })

                        i = as.numeric(i)
                        i[y >= lower_i[n_brackets]] <- n_brackets

                        m_i = m[i]
                        c_i = c[i]

                        probDensity_closedBrackets <- (m_i*y + c_i)/N

                        if(!is.null(PDFpareto_lastBracket)){
                                probDensity_openBracket <- PDFpareto_lastBracket(y)
                        }else{
                                probDensity_openBracket <- rep(0, length(y))
                        }

                        probDensity = ifelse(i < n_brackets,
                                             probDensity_closedBrackets,
                                             probDensity_openBracket)

                        probDensity = ifelse(y < lower_i[1], 0, probDensity)

                        probDensity
                }


                mean = integrate(f = function(y) y*pdf_MCIB(y),
                                    lower = first(lower_i),
                                    upper = pareto_upper_bound,
                                    subdivisions = 2000,
                                    rel.tol = 1e-10,
                                    stop.on.error = FALSE)$value

                mean

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        mean_result <- future_map_dbl(.x = data_split,
                                      .f = mean_MCIB,
                                      topBracket_method_chosen = topBracket_method_chosen,
                                      .progress = T,
                                      .options = future_options(globals = c("topBracket_method_chosen"),
                                                                packages = c("tableInequality", "data.table"))) %>%
                tibble(ID = names(.), mean = .)

        if(is.null(groups)){
                mean_result <- mean_result %>%
                        dplyr::select(-ID)
        }else{
                mean_result <- mean_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        mean_result

}


