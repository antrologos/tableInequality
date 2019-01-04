#' @export

calc_gini_parabolicInterp <- function(data_pnad,
                                      groups = NULL,
                                      combinacao_de_interpolacoes = "mean"){

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


        bin_mids_unbound <- with(data_pnad, {
                getMids(ID = ID,
                        hb = n,
                        lb = min_faixa,
                        ub = max_faixa,
                        alpha_bound = numeric(0))
        })

        colMax <- function(...){
                apply(cbind(...), 1, max, na.rm=T)
        }

        data_pnad <- data_pnad %>% mutate(ID = as.character(ID))

        midPoints_highestBin <- data_pnad %>%
                ungroup() %>%
                filter(is.na(max_faixa)) %>%
                dplyr::select(ID, min_faixa) %>%
                rename(lowerbound = min_faixa) %>%
                left_join(y = tibble(ID              = data_pnad$ID %>% unique() %>% as.character(), #bin_mids_unbound$mids$ID %>% unique() %>% as.character(),
                                     alpha           = bin_mids_unbound$alpha,
                                     alpha_bounded_1 = colMax(1,alpha),
                                     alpha_bounded_2 = colMax(2,alpha)),
                          by = "ID") %>%
                mutate(aritMean = lowerbound * (alpha_bounded_2/(alpha_bounded_2 - 1)),
                       median   = lowerbound*(2^(1/alpha_bounded_1)),
                       geoMean  = lowerbound*(exp(1/alpha_bounded_1)),
                       HarMean  = lowerbound*(1 + (1/alpha_bounded_1)))

        pnads_midpoints <- data_pnad %>%
                ungroup() %>%
                left_join(y  = midPoints_highestBin %>% dplyr::select(ID, aritMean:HarMean),
                          by = "ID") %>%
                mutate(arit_midPoints = (min_faixa + max_faixa)/2,
                       geom_midPoints = (min_faixa*max_faixa)^(1/2),

                       arit_midPoints_aritMean = ifelse(is.na(arit_midPoints), aritMean, arit_midPoints),
                       arit_midPoints_median   = ifelse(is.na(arit_midPoints), median,   arit_midPoints),
                       arit_midPoints_geoMean  = ifelse(is.na(arit_midPoints), geoMean,  arit_midPoints),
                       arit_midPoints_HarMean  = ifelse(is.na(arit_midPoints), HarMean,  arit_midPoints),

                       geom_midPoints_aritMean = ifelse(is.na(geom_midPoints), aritMean, geom_midPoints),
                       geom_midPoints_median   = ifelse(is.na(geom_midPoints), median,   geom_midPoints),
                       geom_midPoints_geoMean  = ifelse(is.na(geom_midPoints), geoMean,  geom_midPoints),
                       geom_midPoints_HarMean  = ifelse(is.na(geom_midPoints), HarMean,  geom_midPoints)
                ) %>%
                dplyr::select(ID, n, arit_midPoints_aritMean:geom_midPoints_HarMean)

        pnads_midpoints_byGroups <- split(pnads_midpoints, pnads_midpoints$ID)

        var_names <- pnads_midpoints_byGroups[[1]] %>%
                dplyr::select(arit_midPoints_aritMean:geom_midPoints_HarMean) %>%
                names()



        #i=1
        #var = var_names[1]

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        gini_result <- future_map_dfr(1:length(pnads_midpoints_byGroups), .f = function(i){

                nome <- names(pnads_midpoints_byGroups)[i]
                data <- pnads_midpoints_byGroups[[i]] %>%
                        filter(n > 0)

                if(nrow(data) > 1){
                        map(var_names, .f = function(var){

                                n_classes = data$n
                                p_medios  = data[[var]]

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

                                # Interpolando
                                lista_interpolacoes <- list()
                                p_interpolacao = seq(from = 0, to = 1, by = 10e-4)
                                for(i in 1:(length(p)-2)){

                                        p_select = p_interpolacao[p_interpolacao >= p[i] & p_interpolacao <= p[i + 2]]

                                        lista_interpolacoes[[i]] = data.frame(segmento = i,
                                                                              p_select,
                                                                              betas[i, 1] + betas[i, 2]*p_select + betas[i, 3]*(p_select^2))
                                        names(lista_interpolacoes[[i]]) = c("segmento", "p" , paste0("r_",i))

                                }

                                df_interpolacoes <- lapply(lista_interpolacoes, function(x) setNames(x, c("segmento", "p" , "r"))) %>%
                                        do.call("rbind", .) %>%
                                        mutate(segmento = as.factor(segmento)) %>%
                                        filter(complete.cases(.)) %>%
                                        as_tibble() %>%
                                        mutate(r = ifelse(r < 0, 0, r))

                                #df_interpolacoes[r < 0, r := 0 ]

                                # Resultado parcial das interpolacoes
                                #graph_lorenz_severalInterpolations <- graph_lorenz_points +
                                #        geom_line(data = df_interpolacoes,
                                #                  aes(x = p, y = r, color = segmento), size = .8)


                                merge_interpolacoes = lista_interpolacoes[[1]] %>% dplyr::select(-segmento)
                                for(i in 2:length(lista_interpolacoes)){
                                        merge_interpolacoes = merge(merge_interpolacoes,
                                                                    lista_interpolacoes[[i]] %>% dplyr::select(-segmento),
                                                                    by="p",all.x = T, all.y = T)
                                }
                                merge_interpolacoes <- as.matrix(merge_interpolacoes)
                                merge_interpolacoes[which(merge_interpolacoes < 0)] <- 0
                                merge_interpolacoes <- as.data.frame(merge_interpolacoes)

                                # Interpolacoes finais da renda acumulada
                                r_interpolado = as.matrix(merge_interpolacoes[-1])

                                r_interpolado = apply(r_interpolado, 1, eval(parse(text = combinacao_de_interpolacoes)), na.rm = T) # mudar a função aqui de max para mean altera o resultado

                                # Lorenz interpolada
                                #graph_lorenz_finalInterpolation <- graph_lorenz_points +
                                #        geom_line(data = NULL, aes(x = p_interpolacao, y = r_interpolado), size = 1.2, alpha = .5)


                                ############################################################

                                # Gini
                                alpha = NULL
                                for(i in 2:length(r_interpolado)){
                                        # Integracao discreta trapezoidal
                                        alpha[i-1] = ((r_interpolado[i] + r_interpolado[i-1])*(p_interpolacao[i] - p_interpolacao[i-1]))/2
                                }

                                gini = 2*(0.5 - sum(alpha))


                                #estimativas_i <- list(graph_lorenz_points                = graph_lorenz_points,
                                #                      graph_lorenz_severalInterpolations = graph_lorenz_severalInterpolations,
                                #                      graph_lorenz_finalInterpolation    = graph_lorenz_finalInterpolation,
                                #                      data_interpolation                 = tibble(p_interpolacao, r_interpolado),
                                #                      gini                               = gini
                                #)
                                #estimativas[[grupo_i$ID]] <- estimativas_i

                                gini

                                }) %>%

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

        #gini_result <- bind_cols(tibble(ID = names(pnads_midpoints_byGroups)), gini_result)

        if(is.null(groups)){
                gini_result <- gini_result %>%
                        gather(key = `Tipo de Média`,
                               value = "gini",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)

        }else{
                gini_result <- gini_result %>%
                        separate(col = ID, into = groups, sep = "_") %>%
                        gather(key = `Tipo de Média`,
                               value = "gini",
                               arit_midPoints_aritMean:geom_midPoints_HarMean)
        }

        gini_result
}


