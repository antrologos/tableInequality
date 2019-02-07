#' @export

model_quantreg_rsubbins <- function(formula, data_pnad, tau = .5,
                                    groups = NULL, confidence_intervals = F,
                                    binSampleSize = 1000,
                                    CI_lowerBound = .025, CI_upperBound = .975,
                                    nsim_CI = 2000){

        dep   = all.vars(formula[[2]])
        indep = all.vars(formula[[3]])

        if(any(indep %in% groups) | any(groups %in% indep)){
                stop("Group variables cannot be used as independent variables in the model")
        }

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        arrange(ID, min_faixa)
        }


        data_pnad <- prepare_to_regression(data_pnad, indep)
        gc()

        data_split <- split(data_pnad, f = data_pnad$ID)

        #which(names(data_split) == 1431102)

        #data_i = data_split[[1]]

        regquantile <- function(data_i, formula, tau,
                                groups, confidence_intervals,
                                binSampleSize, dep, indep,
                                CI_lowerBound, CI_upperBound,
                                nsim_CI){

                ID_i <- data_i$ID %>% unique()

                data_grouped <- data_i %>%
                        unite(col = varIndepIDs, indep) %>%
                        group_by(varIndepIDs, min_faixa, max_faixa) %>%
                        summarise(n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(varIndepIDs, min_faixa)

                data_p <- data_grouped %>%
                        group_by(varIndepIDs) %>%
                        summarise(n = sum(n)) %>%
                        ungroup() %>%
                        mutate(p = n/sum(n))

                X_groups = separate(data_p[,"varIndepIDs"],
                                    col = varIndepIDs, into = indep)
                formula_tmp = formula
                formula_tmp[[2]] = NULL
                X_groups = model.matrix(formula_tmp, data = X_groups)

                test_indeps <- test_for_regression(data_i, indep)

                if(any(test_indeps)){
                        if(confidence_intervals == F){
                                coefs = matrix(rep(as.numeric(NA), ncol(X_groups)), nrow = 1)
                                colnames(coefs) = colnames(X_groups)
                                coefs <- as_tibble(coefs)
                                coefs <- bind_cols(tibble(ID = ID_i), coefs)
                                return(coefs)
                        }else{
                                coefs = matrix(rep(as.numeric(NA), ncol(X_groups)*3), nrow = 1)

                                colNames_tmp1 <- rep(colnames(X_groups), each = 3)
                                colNames_tmp2 <- rep(paste0("_P", c(CI_lowerBound, .50, CI_upperBound)), length(colnames(X_groups)))

                                colNames <- paste0(colNames_tmp1, colNames_tmp2)

                                colnames(coefs) = colNames
                                coefs <- as_tibble(coefs)
                                coefs <- bind_cols(tibble(ID = ID_i), coefs)
                                return(coefs)
                        }
                }

                data_split_indep <- split(data_grouped, data_grouped$varIndepIDs)

                rsub_fits <- lapply(data_split_indep, function(data_x){
                        try(rsubbins(bEdges  = c(1, data_x$max_faixa),
                                     bCounts = c(0,data_x$n)), silent = T)
                })

                quantile_functions <- lapply(rsub_fits, function(x) {
                        if("try-error" %in% class(x)){
                                function(p) rep(NA, length(p))
                        }else{
                                Vectorize(tableInequality:::inverse(x$rsubCDF,
                                                                           lower = 1,
                                                                           upper = x$E,
                                                                           extendInt = "yes"))
                        }
                })

                generate_simData_group = function(i, data_p, rsub_fits, quantile_functions, indep, nsims_i = NULL){

                        #print(i)
                        if(is.null(nsims_i)){
                                nsims_i <- data_p$n[i]
                        }

                        if("try-error" %in% class(rsub_fits[[i]])|nsims_i == 0){
                                return(NULL) #################################################TESTAR
                        }

                        X_values = separate(data_p[i,"varIndepIDs"],
                                            col = varIndepIDs, into = indep)

                        data_sim_i <- bind_rows(replicate(nsims_i, X_values , simplify = FALSE)) %>%
                                mutate(y = quantile_functions[[i]](runif(nsims_i)),
                                       wgt := rsub_fits[[i]]$rsubPDF(y) * data_p[i,]$p)

                        data_sim_i
                }

                get_coef <- function(data_p,
                                     formula,
                                     rsub_fits,
                                     tau,
                                     quantile_functions,
                                     indep,
                                     nsims_i){

                        data_sim = map_dfr(.x = 1:length(rsub_fits),
                                           .f  = function(x) generate_simData_group(i = x,
                                                                                    data_p = data_p,
                                                                                    rsub_fits = rsub_fits,
                                                                                    quantile_functions = quantile_functions,
                                                                                    indep = indep, nsims_i = nsims_i))
                        setDT(data_sim)
                        y   = data_sim[ , eval(formula[[2]])]
                        wgt = data_sim$wgt
                        X   = model.matrix(formula, data_sim)

                        obj = function(b, tau){
                                Indicador1 = as.numeric(y >= as.numeric(X%*%b))
                                Indicador2 = as.numeric(y < as.numeric(X%*%b))

                                sum(wgt*Indicador1*tau*abs(y - as.numeric(X%*%b))) + sum(wgt*Indicador2*(1 - tau)*abs(y - as.numeric(X%*%b)))
                        }

                        estimated_values <- try( nlm(f       = obj,
                                                     p       = rep(1, ncol(X)),
                                                     tau     = tau,
                                                     steptol = 1e-10,
                                                     gradtol = 1e-10,
                                                     iterlim = 1000),
                                                 silent = T)

                        if("try-error" %in% class(estimated_values)){
                                return(rep(NA, ncol(X)))
                        }else{
                                if(estimated_values$code == 5){
                                        return(rep(NA, ncol(X)))
                                }else{
                                        return(estimated_values$estimate)
                                }
                        }
                }


                if(confidence_intervals == F){

                        coefs <- get_coef(tau = tau,
                                          data_p    = data_p,
                                          formula   = formula,
                                          rsub_fits = rsub_fits,
                                          quantile_functions = quantile_functions,
                                          indep   = indep,
                                          nsims_i = binSampleSize)

                        coefs <- matrix(coefs, nrow = 1)

                        colnames(coefs) = colnames(X_groups)
                        coefs <- as_tibble(coefs)
                        coefs <- bind_cols(tibble(ID = ID_i), coefs)
                        return(coefs)

                }else{

                        coeficients_distribution <- future_map(.x = 1:nsim_CI,
                                                               .f = function(x){
                                                                       get_coef(tau       = tau,
                                                                                data_p    = data_p,
                                                                                formula   = formula,
                                                                                rsub_fits = rsub_fits,
                                                                                quantile_functions = quantile_functions,
                                                                                indep     = indep,
                                                                                nsims_i   = NULL)},
                                                               .progress = T)

                        coeficients_distribution <- do.call(rbind, coeficients_distribution)

                        coef_ConfInterval <- apply(coeficients_distribution, 2, quantile, probs = c(CI_lowerBound, .5, CI_upperBound))
                        dim(coef_ConfInterval) = NULL
                        coef_ConfInterval = matrix(coef_ConfInterval, nrow = 1)


                        colNames_tmp1 <- rep(colnames(X_groups), each = 3)
                        colNames_tmp2 <- rep(paste0("_P", c(CI_lowerBound, .50, CI_upperBound)), length(colnames(X_groups)))
                        colNames      <- paste0(colNames_tmp1, colNames_tmp2)

                        colnames(coef_ConfInterval) = colNames
                        coef_ConfInterval <- as_tibble(coef_ConfInterval)
                        coef_ConfInterval <- bind_cols(tibble(ID = ID_i), coef_ConfInterval)
                        return(coef_ConfInterval)
                }

        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        #test = NULL
        #for(i in 1:30){
        #        print(i)
        #        test[[i]] <- regquantile(data_i = data_split[[i]],
        #                                 formula = formula,
        #                                 tau = tau,
        #                                 groups = groups,
        #                                 confidence_intervals = confidence_intervals,
        #                                 binSampleSize = binSampleSize,
        #                                 CI_lowerBound = CI_lowerBound,
        #                                 CI_upperBound = CI_upperBound,
        #                                 nsim_CI = nsim_CI)
        #}

        #data_split2 = data_split[1:10]

        regquantile_result <- future_map(.x = data_split,
                                         .f = regquantile,
                                         formula = formula,
                                         tau = tau,
                                         groups = groups,
                                         confidence_intervals = confidence_intervals,
                                         binSampleSize = binSampleSize,
                                         CI_lowerBound = CI_lowerBound,
                                         CI_upperBound = CI_upperBound,
                                         dep = dep, indep = indep,
                                         nsim_CI = nsim_CI,
                                         .progress = T,
                                         .options = future_options(packages = c("tidyr",
                                                                                "dplyr",
                                                                                "binsmooth",
                                                                                "quantreg",
                                                                                "data.table")
                                         ))

        regquantile_result <- do.call(bind_rows, regquantile_result)

        if(is.null(groups)){
                regquantile_result <- regquantile_result %>%
                        dplyr::select(-ID)
        }else{
                regquantile_result <- regquantile_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        regquantile_result

}
