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
                        group_by(ID, faixas_renda) %>%
                        summarise(min_faixa = min(min_faixa),
                                  midpoint  = mean(midpoint, na.rm = T),
                                  max_faixa = max(max_faixa),
                                  n         = sum(n)) %>%
                        ungroup() %>%
                        arrange(ID, min_faixa)
        }

        data_split <- split(data_pnad, f = data_pnad$ID)

        #data_i = data_split[[9]]

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

                # Interpolando
                make_parabolafunctions <- function(i, x, beta_matrix){
                        env_tmp       <- new.env(parent = emptyenv())
                        env_tmp$linha <- i
                        env_tmp$beta  <- beta_matrix

                        f = function(x) env_tmp$beta[env_tmp$linha,1] + env_tmp$beta[env_tmp$linha,2]*x + env_tmp$beta[env_tmp$linha,3]*(x^2)
                        f
                }

                functions <- list()
                for(j in 1:nrow(betas)){
                        functions[[j]] <- make_parabolafunctions(i = j,
                                                                 x = x,
                                                                 beta_matrix = betas)
                }

                size = length(p)
                threasholds <- cbind(1:(size-2),3:size)
                p_threasholds <- p[threasholds]
                dim(p_threasholds) = dim(threasholds)

                is_between <-function(x){
                        which(x >= p_threasholds[,1] & x <= p_threasholds[,2])
                }

                lorenz_i <- function(y){

                        functionsNumbers <- is_between(y)

                        lorenz_value_temp = NULL

                        count <- 1
                        for(functionsNumber in functionsNumbers){
                                lorenz_value_temp[count] <- functions[[functionsNumber]](y)
                                count = count + 1
                        }

                        expression = paste(combine_function,"(",paste(lorenz_value_temp, collapse = ", "),")")

                        eval(expr = parse(text = expression))
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


