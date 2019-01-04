#' @export

varDecomp_ordProbit_defCutsHeteroskedastic = function(formula, data_pnad, groups){

        model_df <- model_ordProbit_defCutsHeteroskedastic(formula = formula,
                                                           data_pnad = data_pnad,
                                                           groups = groups)

        dep   = as.character(formula[[2]])
        indep = formula[[3]] %>%
                as.character() %>%
                str_split(" ") %>%
                unlist() %>%
                str_replace_all("[+]", "") %>%
                .[nchar(.)>=1]

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)

                model_df <- model_df %>%
                        mutate(ID = 1)
        }else{
                data_pnad <-
                        data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by_at(c("ID",dep,indep)) %>%
                        summarise(min_faixa = min(min_faixa),
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)

                model_df <-
                        model_df %>%
                        unite(col = ID, groups) %>%
                        arrange(ID)
        }

        data_split  <- split(data_pnad, f = data_pnad$ID)
        model_split <- split(model_df,  f = model_df$ID)

        get_beta = function(model_df_i){
                beta_i_df <- model_df_i %>%
                        dplyr::select(starts_with("beta_")) %>%
                        setNames(str_replace_all(names(.), "beta_", ""))

                beta_i        <- as.numeric(beta_i_df)
                names(beta_i) <- names(beta_i_df)

                beta_i[complete.cases(beta_i)]
        }

        get_lambda = function(model_df_i){
                lambda_i_df <- model_df_i %>%
                        dplyr::select(starts_with("lambda_")) %>%
                        setNames(str_replace_all(names(.), "lambda_", ""))

                lambda_i        <- as.numeric(lambda_i_df)
                names(lambda_i) <- names(lambda_i_df)

                lambda_i[complete.cases(lambda_i)]
        }


        var_decomp_i = function(i){
                #i = 1
                #print(i)

                ID_i       = names(data_split)[i]
                data_i     = data_split[[ID_i]]
                model_df_i = model_split[[ID_i]]

                beta_i   <- get_beta(model_df_i)
                lambda_i <- get_lambda(model_df_i)

                # Removing empty categories
                for(j in 1:length(indep)){
                        test_j <- data_i %>%
                                group_by_at(indep[j]) %>%
                                summarise(n = sum(n))
                        cat_problem = test_j[[1]][which(test_j$n == 0)]
                        if(length(cat_problem > 0)){
                                data_i <- data_i %>%
                                        filter_at(vars(indep[j]), any_vars(. != cat_problem))
                        }
                }

                X_i = model.matrix(formula, data = data_i)
                beta_i   = beta_i[colnames(X_i)]
                lambda_i = lambda_i[colnames(X_i)]

                p_i = data_i$n/sum(data_i$n)

                m_j  = as.numeric(X_i%*%beta_i)
                mu   = sum(p_i*m_j)
                s2_j = as.numeric(exp(X_i%*%lambda_i))

                between = sum(p_i*(m_j - mu)^2)
                within  = sum(p_i*s2_j)

                tibble(ID      = ID_i,
                       var     = between + within,
                       between = between,
                       within  = within,
                       R2      = between/var)

                #i =  i + 1
        }

        var_decomp_result <- map_dfr(.x = 1:length(data_split), var_decomp_i)

        var_decomp_result %>%
                separate(col = ID, into = groups, sep =  "_")
}
