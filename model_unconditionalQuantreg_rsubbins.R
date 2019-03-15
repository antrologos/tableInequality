#' @export

#c1970 %>%
#        group_by(municipalityCurrent) %>%
#        summarise(N = sum(wgtperson)) %>%
#        top_n(10, N)

#data = c1970[municipalityCurrent == 2738835]

#write_csv(data, "c:/dropbox/dados_teste_1970.csv")

#formula = log(totalIncome2010Values) ~ male + as.factor(educationAttainment_aggreg)
#weight = "wgtperson"
#tau = seq(.1,.9,.1)
#kernel = "gaussian"
#bandwidth = "nrd0"
#returnData = F

wtd.urq <- function (formula, data, tau = .5, weight = NULL, kernel = "gaussian",
                     bandwidth = "nrd0",
                     returnData = F){

        data <- data %>%
                dplyr::select(c(all.vars(formula), weight)) %>%
                filter(complete.cases(.)) %>%
                as.data.frame()

        dep = all.vars(formula)[1]

        if(is.null(weight)){
                weight <- "weight"
                data[[weight]] <- 1
        }

        y   = model.matrix(formula[1:2], data = data)[,2]
        X   = model.matrix(formula[c(1,3)], data = data)
        wgt = data[[weight]]

        q = Hmisc::wtd.quantile(y, wgt, probs = tau)

        density_values <- density(x = y,
                                  kernel = kernel,
                                  bw =  bandwidth,
                                  weights = wgt/sum(wgt))

        pdf <- function(x) approx(density_values$x, density_values$y, xout = x, rule = 2)$y
        fq <- pdf(q)
        w = sqrt(wgt)

        indicator <- function(condition) ifelse(condition, 1, 0)
        c1 = NULL
        for (i in 1:length(tau)) {
                RIF = q[i] + ((tau[i] - indicator(y < q[i]))/fq[i])
                c = crossprod(solve(crossprod(w*X)), crossprod(w*X, w*RIF)) %>% t()
                c1 = rbind(c1, c)
        }

        RIF = t(c1)
        colnames(RIF) <- paste("tau=", tau)

        if(returnData == T){
                fit = list(coefficients = RIF, tau = tau, formula = formula, data = data)
        }else{
                fit = list(coefficients = RIF, tau = tau, formula = formula)
        }

        fit$call <- sys.call()
        class(fit) <- "wtd.urq"
        return(fit)
}


tau = seq(.01,.99,.025)
teste = quantreg::rq(formula,
                     data = data,
                     tau = tau,
                     weights = wgtperson)
plot(teste)

teste2 = urq(formula, data, tau = tau)
plot(tau, teste2$coefficients[1,], type = "o")
plot(tau, teste2$coefficients[2,], type = "o")
plot(tau, teste2$coefficients[3,], type = "o")
plot(tau, teste2$coefficients[4,], type = "o")
plot(tau, teste2$coefficients[5,], type = "o")

tau = seq(.05, .95, .05)
teste3 = wtd.urq(formula, data, tau = tau, weight = "wgtperson", bandwidth = "nrd")
plot(tau, teste3$coefficients[1,], type = "o")
plot(tau, teste3$coefficients[2,], type = "o")
plot(tau, teste3$coefficients[3,], type = "o")
plot(tau, teste3$coefficients[4,], type = "o")
plot(tau, teste3$coefficients[5,], type = "o")




#

plot(seq(.01,.95,.05), teste$coefficients[2, ], type = "l")
