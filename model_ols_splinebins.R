#rm(list=ls())
library(data.table)
library(tidyverse)
library(binsmooth)
library(furrr)

#load("C:/Dropbox/Rogerio/Formacao/Pos-Doc/Exercicios com PNAD 1968-1973/2-Analises/Ambiente - TESTE 20 - testando model_ols_splinebins.RData")
load("E:/Dropbox-Ro/Dropbox/Rogerio/Formacao/Pos-Doc/Exercicios com PNAD 1968-1973/2-Analises/Ambiente - TESTE 20 - testando model_ols_splinebins.RData")

rm(list = setdiff(ls(), c("c1970", "c1970_aggreg")))
gc()

c1970        = setDT(c1970)
c1970_aggreg = setDT(c1970_aggreg)

formula_micro = log(totalIncome2010Values) ~ male + as.factor(educationAttainment_aggreg)
formula = log(y) ~ male + as.factor(educationAttainment_aggreg)

#=============================================================================
# Com microdados

regression = function(formula, data_i, weight = NULL){

        if(is.null(weight)){
                vars = all.vars(formula)
        }else{
                vars = c(all.vars(formula), weight)
        }

        data = data_i %>%
                dplyr::select(vars) %>%
                filter(complete.cases(.))

        X = try(model.matrix(formula, data = data), silent = T)

        if("try-error" %in% class(X)){
                return(matrix(NA, nrow = 1, ncol = 5))
        }

        if(ncol(X) < 5){
                return(matrix(NA, nrow = 1, ncol = 5))
        }

        y = data[[vars[1]]]

        if(is.null(weight)){
                w = rep(1, length(y))
        }else{
                w = data[[weight]]
        }

        sqrt_w = sqrt(w)

        beta = try( {crossprod(solve(crossprod(sqrt_w*X)), crossprod(sqrt_w*X, sqrt_w*log(y))) %>% t()} ,
                    silent = T)

        if("try-error" %in% class(beta)){
                beta = matrix(NA, nrow = 1, ncol = 5)
        }else{

                y_hat = X%*%as.numeric(beta)
                e = as.numeric(log(y) - y_hat)

                var_cov = solve(crossprod(sqrt_w*X)) %*%  t(e*sqrt_w*X) %*% (e*sqrt_w*X) %*% t(solve(crossprod(sqrt_w*X)))

                sd = sqrt(diag(var_cov))

                if(any(is.na(sd))|any(is.nan(sd))|any(!is.finite(sd))){
                        beta = matrix(NA, nrow = 1, ncol = 5)
                }
        }

        beta
}

selected_munic = unique(c1970$municipalityCurrent)

betas = list()
for(i in 1:length(selected_munic)){
        print(i)
        munic = selected_munic[i]
        betas[[i]] <- regression(formula = formula_micro,
                                 data_i = c1970[municipalityCurrent == munic], weight = "wgtperson")

        betas[[i]] = cbind(municipalityCurrent = munic , betas[[i]])
}

betas = do.call(what = rbind, args = betas)
betas = as_tibble(betas)

#=============================================================================

formula = log(y) ~ male + as.factor(educationAttainment_aggreg)
data_pnad = c1970_aggreg
groups = "municipalityCurrent"
confidence_intervals = F
binSampleSize = 1000
CI_lowerBound = .025
CI_upperBound = .975
nsim_CI = 2000

model_ols_splinebins <- function(formula, data_pnad,
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
                        mutate(ID = 1)

                IDs <- unique(data_pnad$ID)

                data_pnad <- data_pnad %>%
                        group_by_at(c("ID", "min_faixa", "max_faixa", indep)) %>%
                        summarise(n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)

        }else{
                data_pnad <-
                        data_pnad %>%
                        unite(col = ID, groups)

                IDs <- unique(data_pnad$ID)

                data_pnad <- data_pnad %>%
                        group_by_at(c("ID", "min_faixa", "max_faixa", indep)) %>%
                        summarise(n = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_pnad <- tableInequality:::prepare_to_regression(data_pnad, indep)
        gc()


        categories <- list()
        for(var in indep){
                categories[[var]] <- unique(data_pnad[[var]])
        }
        X_model <- model.matrix(formula[c(1,3)], data = expand.grid(categories))

        data_split <- split(data_pnad, f = data_pnad$ID)

        #which(names(data_split) == 1431102)

        #data_i = data_split[[1]]

        reg_ols <- function(data_i, formula,
                            groups, confidence_intervals,
                            binSampleSize,
                            CI_lowerBound, CI_upperBound,
                            nsim_CI){

                dep   = all.vars(formula[[2]])
                indep = all.vars(formula[[3]])

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

                test_indeps <- tableInequality:::test_for_regression(data_i, indep)

                if(any(test_indeps)|ncol(X_groups) < ncol(X_model)){
                        if(confidence_intervals == F){
                                coefs = matrix(rep(as.numeric(NA), ncol(X_model)), nrow = 1)
                                colnames(coefs) = colnames(X_model)
                                coefs <- as_tibble(coefs)
                                coefs <- bind_cols(tibble(ID = ID_i), coefs)
                                return(coefs)
                        }else{
                                coefs = matrix(rep(as.numeric(NA), ncol(X_model)*3), nrow = 1)

                                colNames_tmp1 <- rep(colnames(X_model), each = 3)
                                colNames_tmp2 <- rep(paste0("_P", c(CI_lowerBound, .50, CI_upperBound)), length(colnames(X_groups)))

                                colNames <- paste0(colNames_tmp1, colNames_tmp2)

                                colnames(coefs) = colNames
                                coefs <- as_tibble(coefs)
                                coefs <- bind_cols(tibble(ID = ID_i), coefs)
                                return(coefs)
                        }
                }

                data_split_indep <- split(data_grouped, data_grouped$varIndepIDs)

                spline_fits <- lapply(data_split_indep, function(data_x){
                        try(splinebins(bEdges  = c(1, data_x$max_faixa),
                                     bCounts = c(0,data_x$n)), silent = T)
                })

                quantile_functions <- lapply(spline_fits, function(x) {
                        if("try-error" %in% class(x)){
                                function(p) rep(NA, length(p))
                        }else{
                                Vectorize(tableInequality:::inverse(x$splineCDF,
                                                                    lower = 1,
                                                                    upper = x$E,
                                                                    extendInt = "yes"))
                        }
                })

                generate_simData_group = function(i, data_p, spline_fits, quantile_functions, indep, nsims_i = NULL){

                        #print(i)
                        if(is.null(nsims_i)){
                                nsims_i <- data_p$n[i]
                        }

                        if("try-error" %in% class(spline_fits[[i]])|nsims_i == 0){
                                return(NULL) #################################################TESTAR
                        }

                        X_values = separate(data_p[i,"varIndepIDs"],
                                            col = varIndepIDs, into = indep)

                        data_sim_i <- bind_rows(replicate(nsims_i, X_values , simplify = FALSE)) %>%
                                mutate(y   = quantile_functions[[i]](runif(nsims_i)),
                                       wgt = spline_fits[[i]]$splinePDF(y) * data_p[i,]$p)

                        data_sim_i
                }

                get_coef <- function(data_p,
                                     formula,
                                     spline_fits,
                                     quantile_functions,
                                     indep,
                                     nsims_i){

                        data_sim = map_dfr(.x = 1:length(spline_fits),
                                           .f  = function(x) generate_simData_group(i = x,
                                                                                    data_p = data_p,
                                                                                    spline_fits = spline_fits,
                                                                                    quantile_functions = quantile_functions,
                                                                                    indep = indep, nsims_i = nsims_i))

                        y   = model.matrix(formula[1:2], data_sim)[, 2]
                        wgt = data_sim$wgt
                        sqrt_w = sqrt(wgt)
                        X   = model.matrix(formula, data_sim)

                        beta = try( {crossprod(solve(crossprod(sqrt_w*X)), crossprod(sqrt_w*X, sqrt_w*y)) %>% t()} ,
                                    silent = T)

                        if("try-error" %in% class(beta)){
                                return(matrix(NA, nrow = 1, ncol = ncol(X)))
                        }else{

                                y_hat = X%*%as.numeric(beta)
                                e = as.numeric(log(y) - y_hat)

                                var_cov = solve(crossprod(sqrt_w*X)) %*%  t(e*sqrt_w*X) %*% (e*sqrt_w*X) %*% t(solve(crossprod(sqrt_w*X)))

                                sd = sqrt(diag(var_cov))

                                if(any(is.na(sd))|any(is.nan(sd))|any(!is.finite(sd))){
                                        beta = matrix(NA, nrow = 1, ncol = ncol(X))
                                }
                        }

                        beta

                }


                if(confidence_intervals == F){

                        coefs <- get_coef(data_p    = data_p,
                                          formula   = formula,
                                          spline_fits = spline_fits,
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
                                                                       get_coef(data_p    = data_p,
                                                                                formula   = formula,
                                                                                spline_fits = spline_fits,
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
        #        test[[i]] <- reg_ols(data_i = data_split[[i]],
        #                                 formula = formula,
        #                                 groups = groups,
        #                                 confidence_intervals = confidence_intervals,
        #                                 binSampleSize = binSampleSize,
        #                                 CI_lowerBound = CI_lowerBound,
        #                                 CI_upperBound = CI_upperBound,
        #                                 nsim_CI = nsim_CI)
        #}

        #data_split2 = data_split[1:10]

        reg_ols_result <- future_map(.x = data_split,
                                     .f = reg_ols,
                                     formula = formula,
                                     groups = groups,
                                     confidence_intervals = confidence_intervals,
                                     binSampleSize = binSampleSize,
                                     CI_lowerBound = CI_lowerBound,
                                     CI_upperBound = CI_upperBound,
                                     nsim_CI = nsim_CI,
                                     .progress = T,
                                     .options = future_options(globals = c("formula",
                                                                           "groups",
                                                                           "confidence_intervals",
                                                                           "binSampleSize",
                                                                           "CI_lowerBound",
                                                                           "CI_upperBound",
                                                                           "dep",
                                                                           "indep",
                                                                           "nsim_CI"),
                                                               packages = c("tidyr",
                                                                            "dplyr",
                                                                            "binsmooth",
                                                                            "data.table")
                                     ))

        reg_ols_result <- do.call(bind_rows, reg_ols_result)

        if(is.null(groups)){
                reg_ols_result <- reg_ols_result %>%
                        dplyr::select(-ID)
        }else{
                reg_ols_result <- reg_ols_result %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        reg_ols_result

}


formula = log(y) ~ male + as.factor(educationAttainment_aggreg)

betas_table_splinebins_1000 <- model_ols_splinebins(formula = formula,
                                                data_pnad = c1970_aggreg,
                                                groups = "municipalityCurrent",
                                                binSampleSize = 1000)

betas_table_splinebins_2000 <- model_ols_splinebins(formula = formula,
                                                data_pnad = c1970_aggreg,
                                                groups = "municipalityCurrent",
                                                binSampleSize = 2000)

betas_table_splinebins_3000 <- model_ols_splinebins(formula = formula,
                                                data_pnad = c1970_aggreg,
                                                groups = "municipalityCurrent",
                                                binSampleSize = 3000)


banco_betas <- left_join(betas %>%
                                 gather(key = coef, value = beta_micro, -municipalityCurrent),
                         betas_table_splinebins_1000 %>%
                                 rename(male = male1) %>%
                                 gather(key = coef, value = beta_table_1000, -municipalityCurrent) %>%
                                 mutate(municipalityCurrent = as.numeric(municipalityCurrent))) %>%
        left_join(betas_table_splinebins_2000 %>%
                          rename(male = male1) %>%
                          gather(key = coef, value = beta_table_2000, -municipalityCurrent) %>%
                          mutate(municipalityCurrent = as.numeric(municipalityCurrent))) %>%
        left_join(betas_table_splinebins_3000 %>%
                          rename(male = male1) %>%
                          gather(key = coef, value = beta_table_3000, -municipalityCurrent) %>%
                          mutate(municipalityCurrent = as.numeric(municipalityCurrent)))

banco_betas %>%
        ggplot(aes(x = beta_micro, y= beta_table_1000)) +
        geom_point(alpha = .20) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        geom_smooth() +
        facet_wrap(~coef, scales = "free")

Sys.sleep(5)

banco_betas %>%
        ggplot(aes(x = beta_micro, y= beta_table_2000)) +
        geom_point(alpha = .20) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        geom_smooth() +
        facet_wrap(~coef, scales = "free")

Sys.sleep(5)

banco_betas %>%
        ggplot(aes(x = beta_micro, y= beta_table_3000)) +
        geom_point(alpha = .20) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        geom_smooth() +
        facet_wrap(~coef, scales = "free")



