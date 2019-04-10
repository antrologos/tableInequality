#' @export

calc_theil_stepbins <- function(data_pnad, groups = NULL){

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

        data_split <- split(data_pnad, f = data_pnad$ID)

        interpolate <- function(data){
                interpolation <- binsmooth::stepbins(bEdges  = data$max_faixa,
                                                     bCounts = data$n)

                list(PDF = interpolation[[1]],
                     min = min(data$min_faixa),
                     max = interpolation$E)
        }

        calc_mean_from_pdf <- function(pdf_function, xmin, xmax){
                integrate(f = function(x) x*pdf_function(x), lower =  xmin, upper = xmax, subdivisions = 5000)$value
        }

        calc_theil_from_pdf <- function(pdf_function, xmean, xmin, xmax){
                theil_func <- function(x) (x/xmean)*log((x/xmean))

                theil <- integrate(f = function(x) theil_func(x)*pdf_function(x),
                                   lower =  0,
                                   upper = xmax,
                                   subdivisions = 5000)$value
                theil
        }

        calc_theil_from_data <- function(data){
                interpolation <- interpolate(data)

                xmean          <- calc_mean_from_pdf(pdf_function = interpolation$PDF,
                                                     xmin          = interpolation$min,
                                                     xmax          = interpolation$max)

                theil         <- calc_theil_from_pdf(pdf_function = interpolation$PDF,
                                                     xmean        = xmean,
                                                     xmin         = interpolation$min,
                                                     xmax         = interpolation$max)
                theil
        }

        theil_result <- tableInequality:::future_map_parallel(.x = data_split,
                                                              .f = ~calc_theil_from_data(.x),
                                                              .progress = T) %>%
                tibble(ID = names(.),  theil = unlist(.))

        if(is.null(groups)){
                theil_result <- theil_result %>%
                        dplyr::select(theil)
        }else{
                theil_result <- theil_result %>%
                        dplyr::select(ID, theil) %>%
                        separate(col = ID, into = groups, sep = "_")
        }

        theil_result

}

