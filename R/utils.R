# Numeric calculation of the Lorenz Curve
lorenz <- function(x, PDF_func, lowerbound = 0, upperbound) {

        # grand mean
        grid_mean = mvQuad::createNIGrid(dim=1, type="GLe", level=75)
        mvQuad::rescale(grid_mean, domain = matrix(c(lowerbound, upperbound), ncol=2))
        mu = mvQuad::quadrature(f = function(y) y*PDF_func(y), grid = grid_mean)

        # lorenz value for one observation
        lorenz_i = function(y){
                grid_lorenz = mvQuad::createNIGrid(dim=1, type="GLe", level=75)
                mvQuad::rescale(grid_lorenz, domain = matrix(c(lowerbound, y), ncol=2))

                mvQuad::quadrature(f = function(y) (1/mu)*y*PDF_func(y),
                                   grid = grid_lorenz)
        }

        # lorenz value for a vector
        lorenz_vector = Vectorize(lorenz_i)

        lorenz_vector(x)
}

# function for numeric calculation of a inverse function
# (the quantile function is the inverse of the CDF)
inverse = function (f, lower = -100, upper = 100, extendInt = "no") {
        function (y) uniroot((function(x) f(x) - y), lower = lower, upper = upper,  extendInt = extendInt)$root
}



# Numeric calculation of the quantile function
quantile = function(p, CDF_func, max_x){
        # quantile for one observation
        quantile_i = inverse(CDF_func, lower = 0, upper = max_x, extendInt = "yes")

        quantile_vector = Vectorize(quantile_i)

        # quantile for a vector
        #future_map_dbl(p, function(p) quantile_i(p)$root)
        quantile_vector(p)
}


# Gini index
gini_numerical_integration <- function(PDF_func,  CDF_func, max_x){

        # NUMERICAL INTEGRAL - QUADRATURE
        # create grid
        grid_quantiles_to_lorenz = mvQuad::createNIGrid(dim=1, type="nLe", level=25)

        # compute the approximated value of the integral
        lorenz_integral = mvQuad::quadrature(f = function(x) lorenz(quantile(x,
                                                                     CDF_func = CDF_func,
                                                                     max_x = max_x),
                                                            PDF_func = PDF_func,
                                                            lowerbound = 0,
                                                            upperbound = max_x),
                                     grid = grid_quantiles_to_lorenz)

        gini = 1 - 2*lorenz_integral

        gini
}


future_map_parallel <- function(.x, .f, ..., .progress = FALSE, .options = future_options()){

        if(!any(c("multiprocess", "multicore", "multisession", "cluster") %in% class(plan()))){
                plan(multiprocess)
        }

        future_map(.x, .f, ..., .progress = FALSE, .options = future_options())
}


slopes_MCIB = function(lower_i, upper_i, n_i){

        b_i = seq_along(n_i)
        freq_per_money <- n_i/(upper_i - lower_i) #freq_per_money

        # Estimating the slopes
        M_data <- tibble(lower  = lower_i,
                         upper  = upper_i,
                         f_prev   = c(NA,freq_per_money[-last(b_i)]),
                         f_actual = freq_per_money,
                         f_next   = c(freq_per_money[-1],NA),

                         # Slope from the previous bracket to the actual
                         m_1_0    = (f_actual - f_prev)/(upper - lower),

                         # Slope from the actual bracket to the next one
                         m_2_1    = (f_next   - f_actual)/(upper - lower),

                         # Final value = the average
                         m = (m_2_1+m_1_0)/2) %>%

                # For brackets with no neighbours, we do not use the average
                mutate(m = ifelse(is.na(m) & is.na(m_1_0) & !is.na(m_2_1), m_2_1, m),
                       m = ifelse(is.na(m) & is.na(m_2_1) & !is.na(m_1_0), m_1_0, m))

        # Slopes
        m = M_data$m

        # Minor correction
        higher_both_neighbours = with(M_data, {
                (f_actual > f_prev) & (f_actual > f_next)
        })

        lower_both_neighbours = with(M_data, {
                (f_actual < f_prev) & (f_actual < f_next)
        })

        m[higher_both_neighbours] = 0
        m[lower_both_neighbours]  = 0
        m[1] = 0

        m
}

intercepts_MCIB = function(lower_i, upper_i, n_i, slopes){
        # Intercepts
        c = n_i/(upper_i - lower_i) - slopes*((upper_i + lower_i)/2)
        c
}


estimate_slope_and_intercept_MCIB = function(lower_i = lower_i, upper_i = upper_i, n_i = n_i){

        N = sum(n_i)

        alpha_pareto <- getMids(ID = "1",
                                hb = n_i,
                                lb = lower_i,
                                ub = upper_i,
                                alpha_bound = 2)$alpha

        beta_pareto = last(lower_i)
        pareto_upper_bound = exp( log(beta_pareto) - log(1 - 0.9995)/alpha_pareto)

        # Initial values for slope and intercept
        m_initial = slopes_MCIB(lower_i = lower_i, upper_i = upper_i, n_i = n_i)

        # A function to check if the estimated slopes produce negative probability densities.
        pdf_MCIB_slopeTest = function(m){

                c_initial = intercepts_MCIB(lower_i = lower_i,
                                            upper_i = upper_i,
                                            n_i = n_i,
                                            slopes = m)

                # The values we will use to evaluate the probability densities are at the
                # edge of each bracket.
                y_lower = lower_i + .000000000001
                y_upper = upper_i[-length(upper_i)] - .000000000001

                y = c(y_lower, y_upper) %>% sort()

                n_brackets = length(lower_i)
                i = sapply(y, function(x){
                        which((x >= lower_i) & (x < upper_i))
                })

                i = as.numeric(i)
                i[y >= lower_i[n_brackets]] <- n_brackets

                m_i = m[i]
                c_i = c_initial[i]

                probDensity_closedBrackets = (m_i*y + c_i)/N

                f_pareto_lastBracket = function(y){
                        beta_pareto = lower_i[n_brackets]
                        log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                        exp(log_PDF)
                }

                probDensity_openBracket = (n_i[n_brackets]/N)*f_pareto_lastBracket(y)

                probDensity = ifelse(i < n_brackets,
                                     probDensity_closedBrackets,
                                     probDensity_openBracket)

                probDensity = ifelse(y < lower_i[1], 0, probDensity)

                probDensity_test = NA
                probDensity_test[probDensity <  0] <- -1
                probDensity_test[probDensity >= 0] <-  1

                slope_test = tibble(probDensity_test, i) %>%
                        group_by(i) %>%
                        summarise(test = min(probDensity_test)) %>%
                        arrange(i) %>%
                        .$test

                slope_test
        }

        pdf_MCIB_correctSlope = function(m_particular, m_order){

                m = m_initial

                m[m_order] = m_particular

                # The values we will use to evaluate the probability densities are at the
                # edge of each bracket.
                y_lower = lower_i + .000000000001
                y_upper = upper_i[-length(upper_i)] - .000000000001

                y = c(y_lower, y_upper) %>% sort()

                n_brackets = length(lower_i)
                i = sapply(y, function(x){
                        which((x >= lower_i) & (x < upper_i))
                })

                i = as.numeric(i)
                i[y >= lower_i[n_brackets]] <- n_brackets

                c = intercepts_MCIB(lower_i = lower_i,
                                    upper_i = upper_i,
                                    n_i     = n_i,
                                    slopes = m)

                m_i = m[i]
                c_i = c[i]

                probDensity_closedBrackets = (m_i*y + c_i)/N

                f_pareto_lastBracket = function(y){
                        beta_pareto = lower_i[n_brackets]
                        log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                        exp(log_PDF)
                }

                probDensity_openBracket = (n_i[n_brackets]/N)*f_pareto_lastBracket(y)

                probDensity = ifelse(i < n_brackets,
                                     probDensity_closedBrackets,
                                     probDensity_openBracket)

                probDensity = ifelse(y < lower_i[1], 0, probDensity)

                probDensity_test = tibble(probDensity, i, y) %>%
                        dplyr::filter(i == m_order) %>%
                        dplyr::filter(probDensity < 0)

                if(nrow(probDensity_test) > 1) stop("This segment only gives negative densities")

                y = probDensity_test$y
                p = probDensity_test$probDensity

                new_m = m[m_order]

                n  = n_i[m_order]
                up = upper_i[m_order]
                lw = lower_i[m_order]

                infinitesimal = .Machine$double.eps^0.75

                new_m = (infinitesimal -n/(up - lw))/(y - ((up + lw)/2))
                new_m
        }

        m_generate_positive = pdf_MCIB_slopeTest(m_initial)>0
        problematic_m <- which(!m_generate_positive)

        m_corrected <- m_initial

        if(length(problematic_m) >0 ){
                for(k in 1:length(problematic_m)){
                        m_order = problematic_m[k]
                        m_particular = m_corrected[m_order]
                        m_corrected[m_order] = pdf_MCIB_correctSlope(m_particular = m_particular, m_order = m_order)
                }
        }

        new_test <- pdf_MCIB_slopeTest(m_corrected)>0
        new_test <- new_test[-which(is.na(upper_i))]
        test_stillProducesNegative = any(!(new_test))

        #if(test_stillProducesNegative == T){
        #        stop("The slopes could could not be constrained to produce only positive probability densities")
        #}

        if(test_stillProducesNegative == T){
                m_generate_positive = pdf_MCIB_slopeTest(m_corrected) > 0
                problematic_m <- which(!m_generate_positive)

                if(length(problematic_m) >0 ){
                        for(k in 1:length(problematic_m)){
                                m_order = problematic_m[k]
                                m_particular = m_corrected[m_order]
                                m_corrected[m_order] = pdf_MCIB_correctSlope(m_particular = m_particular, m_order = m_order)
                        }
                }
        }

        if(test_stillProducesNegative == T){
                message("The slopes could could not be constrained to produce only positive probability densities")
        }

        # re-estimating the intercepts
        c = intercepts_MCIB(lower_i = lower_i, upper_i = upper_i, n_i = n_i,
                            slopes = m_corrected)

        list(m = m_corrected, c =c)
}



check_known_groupMeans_DF <- function(data_pnad, groups = NULL, known_groupMeans = NULL){

        if(!is.null(known_groupMeans)){

                checkMeansDF <- is.data.frame(known_groupMeans)
                if(checkMeansDF == FALSE){
                        stop("'known_groupMeans' is not a data.frame")
                }


                if(!is.null(groups)){
                        checkIDvars <- all(groups %in% names(known_groupMeans))
                        if(checkIDvars == FALSE){
                                stop("Groups variables are not contained in 'known_groupMeans'")
                        }

                        known_groupMeans <- known_groupMeans %>%
                                unite(col = ID, groups)

                        checkUniqueIDs = length(known_groupMeans$ID) == length(unique(known_groupMeans$ID))
                        if(checkUniqueIDs == FALSE){
                                stop("Group IDs defined by the group variables are not unique id the 'known_groupMeans' data.frame")
                        }

                }else{
                        known_groupMeans$ID = 1
                }

                testDataID = all(unique(data_pnad$ID) %in% known_groupMeans$ID)
                if(testDataID == FALSE){
                        stop("There are groups listed in the data, but not in the 'known_groupMeans' data.frame")
                }

                testMeanID = all(known_groupMeans$ID %in% unique(data_pnad$ID))
                if(testMeanID == FALSE){
                        stop("There are groups listed in the 'known_groupMeans' data.frame, but not in the data")
                }

                checkMeanVar <-"mean" %in% names(known_groupMeans)
                if(checkMeanVar == FALSE){
                        stop("'known_groupMeans' must have a column named 'mean', which contains the group means")
                }

                return(known_groupMeans)

        }else{
                NULL

        }
}



make_pdf_gpinter <- function(data_i, known_groupMeans_checked){

        ID_i = data_i$ID %>% unique()
        p    = data_i$n/sum(data_i$n)

        prob_quantis <- tibble(p_cum = c(0, cumsum(p[-length(p)])),
                               q     =  data_i$min_faixa)

        prob_quantis <- prob_quantis %>%
                dplyr::filter(p_cum < 1)
        min_p <- last(which(prob_quantis$p_cum == 0))
        prob_quantis <- prob_quantis[min_p:nrow(prob_quantis),]

        if(sum(prob_quantis$p_cum > 0) < 3){
                return(NA)
        }

        if(!is.null(known_groupMeans_checked)){
                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q,
                                       average = known_groupMean_i),
                        silent = T)
        }else{
                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q),
                        silent = T)
        }

        if("try-error" %in% class(pareto_threshold_fitted)){
                return(NA)
        }

        function(y) fitted_density(dist = pareto_threshold_fitted, x = y)


}



make_PDFpareto_lastBracket = function(data_i, topBracket_method_chosen, known_groupMeans_checked, m, c){

        ID_i    = data_i$ID %>% unique()
        lower_i = data_i$min_faixa
        upper_i = data_i$max_faixa
        n_i     = data_i$n

        N = sum(n_i)

        lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
        upper_i = c(lower_i[-1], NA)

        if(is.null(known_groupMeans_checked)){ # If there is no information about the grand mean

                test = is.character(topBracket_method_chosen) & topBracket_method_chosen %in% c("gpinter", "RPME")
                if(test == FALSE){
                        stop("'topBracket_method' must be 'gpinter' or 'RPME'")
                }


                if(topBracket_method_chosen == "RPME"){
                        beta_pareto = last(lower_i)

                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha

                        f_pareto_lastBracket = function(y){
                                # Probability density that uses alpha obtained from robust two-point estimation strategy
                                if(!is.na(alpha_pareto)){
                                        log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                                        PDF     = exp(log_PDF)

                                        (last(n_i)/N)*PDF
                                }else{
                                        0
                                }
                        }
                }

                if(topBracket_method_chosen == "gpinter"){
                        f_pareto_lastBracket = make_pdf_gpinter(data_i, known_groupMeans_checked)
                }

        }else{ # If there is information about the grand mean

                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                group_means_by_integral = {
                        ((m/(3*n_i))*(upper_i^3) + (c/(2*n_i))*(upper_i^2)) - ((m/(3*n_i))*(lower_i^3) + (c/(2*n_i))*(lower_i^2))
                }

                mean_last_group = (1/last(n_i))*(N*known_groupMean_i - sum(n_i*group_means_by_integral, na.rm = T))

                beta_pareto <- last(lower_i)

                if(mean_last_group > beta_pareto){
                        # Estimating alpha using the information about the grand mean
                        alpha_pareto = mean_last_group/(mean_last_group - beta_pareto)
                }else{
                        # If there is any problem or impossibility to estimate alpha using the mean-constrained strategy,
                        # we will use the robust two-point estimate
                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha
                }

                f_pareto_lastBracket = function(y){
                        # Probability density that uses the alpha estimated taking into account the grand mean
                        if(!is.na(alpha_pareto)){
                                log_PDF = log(alpha_pareto) + alpha_pareto*log(beta_pareto) - (alpha_pareto+1)*log(y)
                                PDF     = exp(log_PDF)

                                (last(n_i)/N)*PDF
                        }else{
                                0
                        }
                }

        }

        f_pareto_lastBracket
}


get_pareto_upper_bound = function(data_i,
                                  topBracket_method_chosen,
                                  known_groupMeans_checked,
                                  m, c){

        ID_i    = data_i$ID %>% unique()
        lower_i = data_i$min_faixa
        upper_i = data_i$max_faixa
        n_i     = data_i$n

        N = sum(n_i)

        lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
        upper_i = c(lower_i[-1], NA)

        beta_pareto = last(lower_i)

        if(is.null(known_groupMeans_checked)){ # If there is no information about the grand mean

                if(topBracket_method_chosen == "RPME"){

                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha
                }


        }else{
                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                group_means_by_integral = {
                        ((m/(3*n_i))*(upper_i^3) + (c/(2*n_i))*(upper_i^2)) - ((m/(3*n_i))*(lower_i^3) + (c/(2*n_i))*(lower_i^2))
                }

                mean_last_group = (1/last(n_i))*(N*known_groupMean_i - sum(n_i*group_means_by_integral, na.rm = T))

                if(mean_last_group > beta_pareto){
                        # Estimating alpha using the information about the grand mean
                        alpha_pareto = mean_last_group/(mean_last_group - beta_pareto)
                }else{
                        # If there is any problem or impossibility to estimate alpha using the mean-constrained strategy,
                        # we will use the robust two-point estimate
                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha
                }

        }

        if(topBracket_method_chosen == "gpinter"){
                alpha_pareto = NA
        }

        pareto_upper_bound = exp( log(beta_pareto) - log(1 - 0.995)/alpha_pareto)


        if(is.na(pareto_upper_bound) & last(n_i) > 0){
                pareto_upper_bound = Inf
        }

        if(is.na(pareto_upper_bound) & last(n_i) == 0){
                pareto_upper_bound = last(lower_i)
        }

        pareto_upper_bound
}


make_cdf_gpinter <- function(data_i, known_groupMeans_checked){

        ID_i = data_i$ID %>% unique()
        p    = data_i$n/sum(data_i$n)

        prob_quantis <- tibble(p_cum = c(0, cumsum(p[-length(p)])),
                               q     =  data_i$min_faixa)

        prob_quantis <- prob_quantis %>%
                dplyr::filter(p_cum < 1)
        min_p <- last(which(prob_quantis$p_cum == 0))
        prob_quantis <- prob_quantis[min_p:nrow(prob_quantis),]

        if(sum(prob_quantis$p_cum > 0) < 3){
                return(NA)
        }

        if(!is.null(known_groupMeans_checked)){

                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q,
                                       average = known_groupMean_i),
                        silent = T)
        }else{
                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q),
                        silent = T)
        }

        if("try-error" %in% class(pareto_threshold_fitted)){
                return(NA)
        }

        function(y) fitted_cdf(dist = pareto_threshold_fitted, x = y)
}


make_CDFpareto_lastBracket = function(data_i, topBracket_method_chosen, known_groupMeans_checked, m, c){

        ID_i    = data_i$ID %>% unique()
        lower_i = data_i$min_faixa
        upper_i = data_i$max_faixa
        n_i     = data_i$n

        N = sum(n_i)

        lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
        upper_i = c(lower_i[-1], NA)

        cum_p = cumsum(n_i)/N
        cum_p = c(0,cum_p[-length(cum_p)])

        if(is.null(known_groupMeans_checked)){ # If there is no information about the grand mean

                test = is.character(topBracket_method_chosen) & topBracket_method_chosen %in% c("gpinter", "RPME")
                if(test == FALSE){
                        stop("'topBracket_method' must be 'gpinter' or 'RPME'")
                }

                if(topBracket_method_chosen == "RPME"){
                        beta_pareto = last(lower_i)

                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha

                        f_pareto_lastBracket = function(y){
                                # Cumulative density function that uses alpha obtained from robust two-point
                                # estimation strategy
                                if(!is.na(alpha_pareto)){
                                        CDF_pareto = 1 - (beta_pareto/y)^alpha_pareto
                                        last(cum_p) + (last(n_i)/N)*CDF_pareto
                                }else{
                                        last(cum_p)
                                }
                        }
                }

                if(topBracket_method_chosen == "gpinter"){
                        f_pareto_lastBracket = make_cdf_gpinter(data_i, known_groupMeans_checked)
                }


        }else{ # If there is information about the grand mean

                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                group_means_by_integral = {
                        ((m/(3*n_i))*(upper_i^3) + (c/(2*n_i))*(upper_i^2)) - ((m/(3*n_i))*(lower_i^3) + (c/(2*n_i))*(lower_i^2))
                }

                mean_last_group = (1/last(n_i))*(N*known_groupMean_i - sum(n_i*group_means_by_integral, na.rm = T))

                beta_pareto <- last(lower_i)
                if(mean_last_group > beta_pareto){
                        # Estimating alpha using the information about the grand mean
                        alpha_pareto = mean_last_group/(mean_last_group - beta_pareto)
                }else{
                        # If there is any problem or impossibility to estimate alpha using the mean-constrained strategy,
                        # we will use the robust two-point estimate
                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha
                }

                f_pareto_lastBracket = function(y){
                        # Cumulative density function that uses the alpha estimated taking into account the grand mean
                        if(!is.na(alpha_pareto)){
                                CDF_pareto = 1 - (beta_pareto/y)^alpha_pareto
                                last(cum_p) + (last(n_i)/N)*CDF_pareto
                        }else{
                                last(cum_p)
                        }
                }
        }
        f_pareto_lastBracket
}


make_quantileFunction_gpinter <- function(data_i, known_groupMeans_checked){

        ID_i = data_i$ID %>% unique()
        p    = data_i$n/sum(data_i$n)

        prob_quantis <- tibble(p_cum = c(0, cumsum(p[-length(p)])),
                               q     =  data_i$min_faixa)

        prob_quantis <- prob_quantis %>%
                dplyr::filter(p_cum < 1)
        min_p <- last(which(prob_quantis$p_cum == 0))
        prob_quantis <- prob_quantis[min_p:nrow(prob_quantis),]

        if(sum(prob_quantis$p_cum > 0) < 3){
                return(NA)
        }

        if(!is.null(known_groupMeans_checked)){

                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q,
                                       average = known_groupMean_i),
                        silent = T)
        }else{
                pareto_threshold_fitted <- try(
                        thresholds_fit(p = prob_quantis$p_cum, threshold = prob_quantis$q),
                        silent = T)
        }

        if("try-error" %in% class(pareto_threshold_fitted)){
                return(NA)
        }

        function(p) fitted_quantile(dist = pareto_threshold_fitted, probs = p)
}


make_quantileFunctionPareto_lastBracket = function(data_i, topBracket_method_chosen, known_groupMeans_checked, m, c){

        ID_i    = data_i$ID %>% unique()
        lower_i = data_i$min_faixa
        upper_i = data_i$max_faixa
        n_i     = data_i$n

        N = sum(n_i)

        lower_i = rowMeans(cbind(lower_i,c(NA, upper_i[-length(upper_i)])),na.rm = T)
        upper_i = c(lower_i[-1], NA)

        cum_p = cumsum(n_i)/N
        cum_p = c(0,cum_p[-length(cum_p)])

        if(is.null(known_groupMeans_checked)){ # If there is no information about the grand mean

                test = is.character(topBracket_method_chosen) & topBracket_method_chosen %in% c("gpinter", "RPME")
                if(test == FALSE){
                        stop("'topBracket_method' must be 'gpinter' or 'RPME'")
                }

                if(topBracket_method_chosen == "RPME"){
                        beta_pareto = last(lower_i)

                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha

                        f_pareto_lastBracket = function(p){
                                # quantile function that uses alpha obtained from robust two-point estimation strategy
                                y =  beta_pareto/((1 - (p - last(cum_p))/(last(n_i)/N))^(1/alpha_pareto)) #quantile function for pareto
                                y
                        }

                }

                if(topBracket_method_chosen == "gpinter"){
                        f_pareto_lastBracket = make_quantileFunction_gpinter(data_i, known_groupMeans_checked)
                }


        }else{ # If there is information about the grand mean

                known_groupMean_i = known_groupMeans_checked[known_groupMeans_checked$ID == ID_i,]$mean

                group_means_by_integral = {
                        ((m/(3*n_i))*(upper_i^3) + (c/(2*n_i))*(upper_i^2)) - ((m/(3*n_i))*(lower_i^3) + (c/(2*n_i))*(lower_i^2))
                }

                mean_last_group = (1/last(n_i))*(N*known_groupMean_i - sum(n_i*group_means_by_integral, na.rm = T))

                beta_pareto <- last(lower_i)
                if(mean_last_group > beta_pareto){
                        # Estimating alpha using the information about the grand mean
                        alpha_pareto = mean_last_group/(mean_last_group - beta_pareto)
                }else{
                        # If there is any problem or impossibility to estimate alpha using the mean-constrained strategy,
                        # we will use the robust two-point estimate
                        alpha_pareto <- getMids(ID = ID_i,
                                                hb = n_i,
                                                lb = lower_i,
                                                ub = upper_i,
                                                alpha_bound = 2)$alpha
                }

                f_pareto_lastBracket = function(p){
                        # quantile function function that uses the alpha estimated taking into account the grand mean
                        y =  beta_pareto/((1 - (p - last(cum_p))/(last(n_i)/N))^(1/alpha_pareto)) #quantile function for pareto
                        y
                }

        }

        f_pareto_lastBracket
}
