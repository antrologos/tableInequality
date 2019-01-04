
# Numeric calculation of the Lorenz Curve
lorenz <- function(x, PDF_func, lowerbound = 0, upperbound) {

        # grand mean
        #mu = pracma::integral(function(z){z*PDF_func(z)}, lowerbound, upperbound)
        mu = integrate(f = function(z){z*PDF_func(z)},
                       lower = lowerbound,
                       upper = upperbound,
                       stop.on.error = FALSE,
                       subdivisions = 2000)$value

        # lorenz value for one observation
        #lorenz_i = function(y) { (1/mu)*pracma::integral(function(y){y*PDF_func(y)}, lowerbound, y) }
        lorenz_i = function(y) { (1/mu)*integrate(f = function(y){y*PDF_func(y)},
                                                  lower = lowerbound,
                                                  upper = y,
                                                  stop.on.error = FALSE,
                                                  subdivisions = 2000)$value }

        # lorenz value for a vector
        future_map_dbl(x, lorenz_i)
}

# function for numeric calculation of a inverse function
# (the quantile function is the inverse of the CDF)
inverse = function (f, lower = -100, upper = 100, extendInt = "no") {
        function (y) uniroot((function(x) f(x) - y), lower = lower, upper = upper,  extendInt = extendInt)
}

# Numeric calculation of the quantile function
quantile = function(p, CDF_func, max_x){
        # quantile for one observation
        quantile_i = inverse(CDF_func, lower = 0, upper = max_x, extendInt = "yes")

        # quantile for a vector
        future_map_dbl(p, function(p) quantile_i(p)$root)
}


# Gini index
gini_numerical_integration <- function(PDF_func,  CDF_func, max_x){
        #1 - 2*pracma::integral(function(x){
        #        lorenz(
        #                quantile(x, CDF_func = CDF_func, max_x = max_x),
        #
        #                      PDF_func   = PDF_func,
        #               upperbound = max_x)}, xmin = 0, xmax = 1)


        1 - 2*integrate(f = function(x){
                lorenz(
                        quantile(x, CDF_func = CDF_func, max_x = max_x),
                        PDF_func   = PDF_func,
                        upperbound = max_x)},
                lower = 0,
                upper = 1,
                stop.on.error = FALSE,
                subdivisions = 2000)$value


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
        m = slopes_MCIB(lower_i = lower_i, upper_i = upper_i, n_i = n_i)
        c = intercepts_MCIB(lower_i = lower_i, upper_i = upper_i, n_i = n_i, slopes = m)

        # Now we have to check if the estimated slopes produce negative probability densities.
        # If so, the slopes will be constrained to produce only non-negative densities

        # But the new slopes must have same sign as those initially estimated
        # Here we keep the sign information:
        m_positive = m[-length(m)] > 0

        # A function that is optimized only if all the probabily densities are positive and
        # all the constrained slopes have the same sign as the initial ones
        pdf_MCIB_test = function(m){

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

                probDensity_test = NA
                probDensity_test[probDensity< 0]   <- -1
                probDensity_test[probDensity >= 0] <- 1

                probDensity_test[which(! ( (m > 0) ==  m_positive))] <- -1

                -sum(probDensity_test)
        }

        # Optimizing the test function
        m_corrected = optim(fn = pdf_MCIB_test, par = m[-length(m)])$par
        m[1:(length(m)-1)] <- m_corrected

        # re-estimating the intercepts
        c = intercepts_MCIB(lower_i = lower_i, upper_i = upper_i, n_i = n_i, slopes = m)

        list(m = m, c =c)
}


