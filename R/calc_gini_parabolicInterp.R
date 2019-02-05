#' @export

calc_gini_parabolicInterp <- function(data_pnad,
                                      groups = NULL,
                                      combine_function = "mean"){

        if(!exists("midpoint", where = data_pnad)){
                stop("This method requires a column named 'midpoint', which should cointain the backets' midpoints")

        }

        if(is.null(groups)){
                data_pnad <- data_pnad %>%
                        mutate(ID = 1) %>%
                        arrange(ID, min_faixa)
        }else{
                data_pnad <- data_pnad %>%
                        unite(col = ID, groups) %>%
                        group_by(ID, min_faixa) %>%
                        summarise(
                                midpoint  = mean(midpoint, na.rm = T),
                                max_faixa = max(max_faixa),
                                n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[which(names(data_split)==241190)]]
        #data_i = data_split[[1]]

        gini_parabolas <- function(data_i){

                data_i <- data_i %>%
                        filter(n > 0 )

                if(nrow(data_i) < 3){
                        return(NA)
                }

                n_classes = data_i$n
                p_medios  = data_i$midpoint

                # Calculando renda total de cada classe
                total_classes <- n_classes * p_medios
                #total_classes/sum(total_classes) # apenas para teste - ver participacao na renda por classe

                # Calculando a renda acumulada
                r = c(0, cumsum(total_classes))
                r = r/last(r)

                # Calculando a populacao acumulada
                p = n_classes/sum(n_classes)
                p = c(0, cumsum(p))

                # Lorenz nao interpolada
                #graph_lorenz_points <- ggplot(data = NULL, aes(x = p, y = r)) +
                #        geom_point(size = 3) +
                #        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), lty = 2, size = 1) +
                #        ylim(0,1) +
                #        xlim(0,1) +
                #        theme_classic()

                # Ajustando parabolas
                betas <- matrix(NA, ncol=3, nrow = (length(p)-2))
                r2 <- NULL
                for(i in 1:(length(p)-2)){
                        reg        <- lm(r[i:(i+2)] ~ p[i:(i+2)] + I((p[i:(i+2)])^2))
                        r2[i]      <- summary(reg)$r.squared
                        betas[i, ] <- coef(reg)
                }
                colnames(betas) = c("intercept", "p", "p^2")


                # Ajustando equações lineares (para quando parábolas geram resultados negativos)
                betas_lin <- matrix(NA, ncol=2, nrow = (length(p)-1))
                r2_lin    <- NULL
                for(i in 1:(length(p)-1)){
                        reg_lin    <- lm(r[i:(i+1)] ~ p[i:(i+1)])
                        r2_lin[i]      <- summary(reg_lin)$r.squared
                        betas_lin[i, ] <- coef(reg_lin)
                }
                colnames(betas_lin) = c("intercept", "p")


                # Ajustando exponenciais
                #betas_exp <- matrix(NA, ncol=2, nrow = (length(p)-1))
                #for(i in 1:(length(p)-1)){
                #
                #        dep   = c(r[i], r[i+1])
                #        indep = c(p[i], p[(i+1)])
                #
                #        betas_exp[i, ] = coef(lm(log1p(dep) ~ log1p(indep)))
                #
                #        # Inverso:
                #        #expm1(betas[1] + betas[2]*log1p(indep))
                #}


                # Interpolando
                make_parabolafunctions <- function(i, x, beta_matrix){
                        env_tmp       <- new.env(parent = emptyenv())
                        env_tmp$linha <- i
                        env_tmp$beta  <- beta_matrix

                        f = function(x) env_tmp$beta[env_tmp$linha,1] + env_tmp$beta[env_tmp$linha,2]*x + env_tmp$beta[env_tmp$linha,3]*(x^2)
                        f
                }

                # Interpolando - LINEAR
                make_Linearfunctions <- function(i, x, beta_matrix){
                        env_tmp       <- new.env(parent = emptyenv())
                        env_tmp$linha <- i
                        env_tmp$beta  <- beta_matrix

                        f = function(x) env_tmp$beta[env_tmp$linha,1] + env_tmp$beta[env_tmp$linha,2]*x
                        f
                }

                #make_exponentialfunctions <- function(i, x, beta_matrix){
                #        env_tmp       <- new.env(parent = emptyenv())
                #        env_tmp$linha <- i
                #        env_tmp$beta  <- beta_matrix
                #
                #       f = function(x) {
                #
                #                beta0 = env_tmp$beta[env_tmp$linha,1]
                #                beta1 = env_tmp$beta[env_tmp$linha,2]
                #
                #                expm1(beta0 + beta1*log1p(x))
                #        }
                #
                #        f
                #}


                functions <- list()
                for(j in 1:nrow(betas)){
                        functions[[j]] <- make_parabolafunctions(i = j,
                                                                 x = x,
                                                                 beta_matrix = betas)
                }


                lin_functions <- list()
                for(j in 1:nrow(betas_lin)){
                        lin_functions[[j]] <- make_Linearfunctions(i = j,
                                                                   x = x,
                                                                   beta_matrix = betas_lin)
                }

                #exp_functions <- list()
                #for(j in 1:nrow(betas_lin)){
                #        exp_functions[[j]] <- make_exponentialfunctions(i = j,
                #                                                   x = x,
                #                                                   beta_matrix = betas_exp)
                #}

                size = length(p)
                threasholds <- cbind(1:(size-2),3:size)
                p_threasholds <- p[threasholds]
                dim(p_threasholds) = dim(threasholds)

                threasholds_linExp <- cbind(1:(size-1),2:size)
                p_threasholds_linExp <- p[threasholds_linExp]
                dim(p_threasholds_linExp) = dim(threasholds_linExp)

                is_between <-function(x){
                        which(x >= p_threasholds[,1] & x <= p_threasholds[,2])
                }

                is_between_linExp <-function(x){
                        which(x >= p_threasholds_linExp[,1] & x <= p_threasholds_linExp[,2])
                }

                #y = .4
                lorenz_i <- function(y){

                        functionsNumbers <- is_between(y)
                        lorenz_value_temp = NULL

                        count <- 1
                        for(functionsNumber in functionsNumbers){
                                lorenz_value_temp_test1 <- functions[[functionsNumber]](y)
                                #lorenz_value_temp_test2 <- functions[[functionsNumber]](y - .Machine$double.eps^.5)

                                #minimum_parabola <- optimize(f =  functions[[functionsNumber]], interval = c(0, 1), maximum = F)$objective
                                #test_negative = minimum_parabola < 0

                                #if(!test_negative & (lorenz_value_temp_test1 > lorenz_value_temp_test2)){
                                        lorenz_value_temp[count] <- lorenz_value_temp_test1
                                #}else{
                                #        lorenz_value_temp[count] <- NA
                                #}
                                count = count + 1
                        }

                        #if(!all(is.na(lorenz_value_temp))){
                                # Não leva em conta as parábolas que geraram valores negativos ou que não não monotonicamente crescentes
                                expression = paste0(combine_function,"( c(",paste(lorenz_value_temp, collapse = ", "),"), na.rm = T)")

                                lorenz_value <- eval(expr = parse(text = expression))

                        #}else{
                                # Usa interpolação linear quando não houve valores válidos gerados pelas parábolas
                                #functionsNumbers_lin <- is_between_linExp(y)

                                #lorenz_value         <- exp_functions[[functionsNumbers_lin]](y)

                                #if(is.nan(lorenz_value)|is.na(lorenz_value)|!is.finite(lorenz_value)){
                                #        lorenz_value         <- lin_functions[[functionsNumbers_lin]](y)
                                #}
                        #}

                        lorenz_value
                }

                lorenz_parabola = Vectorize(lorenz_i)

                # NUMERICAL INTEGRAL - QUADRATURE
                # create grid
                nw = createNIGrid(dim=1, type="nLe", level=25)

                # compute the approximated value of the integral
                lorenz_integral = quadrature(f = lorenz_parabola,
                                             grid = nw)

                gini = 1 - 2*lorenz_integral

                gini
        }

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        gini_result <- future_map_dfr(.x = data_split, .f = ~gini_parabolas(.x), .progress = T)

        gini_result <- tibble(ID   = rownames(t(gini_result)),
                              gini = t(gini_result)[,1])

        if(is.null(groups)){
                gini_result <- gini_result %>%
                        dplyr::select(gini)
        }else{
                gini_result <- gini_result %>%
                        dplyr::select(ID, gini) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        gini_result
}


