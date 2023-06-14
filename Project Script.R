calc_AIC <- function(model){
  log_likelihood <- model@fit$LLH
  
  # Calculate the number of parameters
  num_parameters <- length(model@fit$coef)
  
  # Calculate the AIC
  aic <- -2 * log_likelihood + 2 * num_parameters
  return(aic)
}
CRPS <- function(x,mean,sigma){
  standart_value <- as.numeric((x-mean)/sigma)
  score <- -sigma*(1/sqrt(pi)-2*dnorm(standart_value)-standart_value*(2*pnorm(standart_value)-1))
  return(score)
}

scatterhist <- function(x, y, xlab = '', ylab = '', xl = NULL, yl = NULL, m = NULL){
  plot(x, y, xlab = xl, ylab = yl, pch = ' .', cex = 1.5, main = m)
  hist(x, freq = FALSE, xlab = xl, main = xl)
  hist(y, freq = FALSE, xlab = yl, main = yl)
}

energy_score <- function(X, y){
  error1 <- 0
  error2 <- 0
  for(i in 1:length(X)){
    error1 <- error1 + abs(X[[i]] - y)
    if(i < length(X)){
      error2 <- error2 + abs(X[[i]] - X[[i+1]])
    }
  }
  error_sum <- error1/length(X) + error2/(2*(length(X)-1))
  return(error_sum)
}


#commodities <- read.csv('C:\\Users\\Acer\\OneDrive - ADA University\\Documents\\GitHub\\cs2\\commodities.csv')
commodities <- read.csv('C:/Users/alkim/OneDrive/Documents/GitHub/cs2/commodities.csv')
commodities$date <- as.Date(commodities$date)
NGas_ts <- xts(commodities$NGas, order.by = commodities$date)
oil_ts <- xts(commodities$oil, order.by = commodities$date)
coal_ts <- xts(commodities$coal, order.by = commodities$date)
####Plotting the time series basic####
plot.xts(NGas_ts, type = 'l')
plot.xts(oil_ts, type = 'l')
plot.xts(coal_ts, type = 'l')
hist(coal_ts, breaks = 50)
hist(NGas_ts, breaks = 50)
hist(oil_ts, breaks = 50)
############AUTOCORRELATION###########
###First Moments###
acf(commodities$coal)
acf(commodities$oil)
acf(commodities$NGas)
###Second Moments###
acf(commodities$coal**2)
acf(commodities$oil**2)
acf(commodities$NGas**2)
############PARTIAL AUTOCORRELATION###########
###First Moments###
pacf(commodities$coal)
pacf(commodities$oil)
pacf(commodities$NGas)
###Second Moments###
pacf(commodities$coal**2)
pacf(commodities$oil**2)
pacf(commodities$NGas**2)

#PART B (you know from where lol)
# Install and load the 'rugarch' package
library(rugarch)
library(ensembleBMA)
library(scoringutils)
library(ggplot2)
library(PerformanceAnalytics)
library(copula)
library(scoringRules)
library(ggridges)
#part b
### Train datalarini cikardik

NGas_ts_train <- NGas_ts[1:2500]
oil_ts_train <- oil_ts[1:2500]
coal_ts_train <- coal_ts[1:2500]
NGas_ts_test <- NGas_ts[2501:2700]
oil_ts_test <- oil_ts[2501:2700]
coal_ts_test <- coal_ts[2501:2700]

# Specify the ARMA(p,q) and GARCH(r,s) orders
p <- c(0:2)  # Autoregressive order
q <- c(0:2)  # Moving average order
r <- c(0:2)  # GARCH order
s <- c(0:2)  # GARCH order for lagged conditional variances

####Models##########
# Specify the specification for the ARMA-GARCH model
AIC_values_NGas <- list()
AIC_values_oil <- list()
AIC_values_coal <- list()
index_list <- list()
counter <- 1
for(i in p){
  for(j in q){
    for(k in r){
      for(l in s){
        if(k == 0 & l == 0){
          next
        }
        else{
          spec <- ugarchspec(mean.model = list(armaOrder = c(i, j)),
                             variance.model = list(garchOrder = c(k, l)))
          fit_NGas <- ugarchfit(spec, NGas_ts_train, solver = 'hybrid')
          fit_oil <- ugarchfit(spec, oil_ts_train, solver = 'hybrid')
          fit_coal <- ugarchfit(spec, coal_ts_train, solver = 'hybrid')
          AIC_values_NGas[[counter]] <- calc_AIC(fit_NGas)
          AIC_values_oil[[counter]] <- calc_AIC(fit_oil)
          AIC_values_coal[[counter]] <- calc_AIC(fit_coal)
          index_list[[counter]] <- paste(i,j,k,l)
          counter <- counter + 1
          print(paste(i,j,k,l))
        }
      }
    }
  }
}
AIC_values_index <- as.data.frame(cbind(unlist(index_list),unlist(AIC_values_NGas),unlist(AIC_values_oil),unlist(AIC_values_coal)))
AIC_values_index$V2 <- as.numeric(AIC_values_index$V2) 
AIC_values_index$V3<- as.numeric(AIC_values_index$V3) 
AIC_values_index$V4<- as.numeric(AIC_values_index$V4) 
##Best model with parameters p=2, q=1, r=1, s=1 for NGas
best_model_NGas <- ugarchspec(mean.model = list(armaOrder = c(2, 1)),
                                           variance.model = list(garchOrder = c(1, 1)))
best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts_train, solver = 'hybrid')

##Best model with parameters p=0, q=1, r=1, s=1 for oil
best_model_oil <- ugarchspec(mean.model = list(armaOrder = c(0, 1)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts_train, solver = 'hybrid')

##Best model with parameters p=2, q=2, r=1, s=1 for coal
best_model_coal <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts_train, solver = 'hybrid')


# Reporting the estimated parameters
print('NGas (2,1,1,1):')
coef(best_model_NGas_fit)
print('oil (0,1,1,1):')
coef(best_model_oil_fit)
print('coal (2,2,1,1):')
coef(best_model_coal_fit)
#part c
best_model_NGas_fit_sigma <- xts(sigma(best_model_NGas_fit)**2, order.by = commodities$date[1:2500])
plot.xts(best_model_NGas_fit_sigma)
best_model_oil_fit_sigma <- xts(sigma(best_model_oil_fit)**2, order.by = commodities$date[1:2500])
plot.xts(best_model_oil_fit_sigma)
best_model_coal_fit_sigma <- xts(sigma(best_model_coal_fit)**2, order.by = commodities$date[1:2500])
plot.xts(best_model_coal_fit_sigma)
#part d-e
#####ITERATIVE FORECASTING WITHOUT REESTIMATION######
fitted_forecast_NGas <- list()
residuals_NGas <- list()
mu_forecast_NGas <- list()
sigma_forecast_NGas <- list()
ar_forecast_NGas <- list()
fitted_forecast_oil <- list()
residuals_oil <- list()
mu_forecast_oil <- list()
sigma_forecast_oil <- list()
ar_forecast_oil <- list()
fitted_forecast_coal <- list()
residuals_coal <- list()
mu_forecast_coal <- list()
sigma_forecast_coal <- list()
ar_forecast_coal <- list()
CRPS_NGas <- list()
CRPS_NGas_ar <- list()
CRPS_oil <- list()
CRPS_oil_ar <- list()
CRPS_coal <- list()
CRPS_coal_ar <- list()
#for (i in 1:(length(NGas_ts) - 2500)){
for (i in 1:200){
  best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts[1:(2500+i-1)], solver = 'hybrid')
  residuals_NGas[[i]] <- residuals(best_model_NGas_fit)
  forecast_NGas <- ugarchforecast(best_model_NGas_fit, n.ahead = 1)
  fitted_forecast_NGas[[i]] <- fitted(forecast_NGas)
  sigma_forecast_NGas[[i]] <- sigma(forecast_NGas)
  ar_model_NGas <- ar.ols(NGas_ts[1:(2500+i-1)],order = 1)
  ar_forecast_NGas[[i]] <- as.numeric(predict(ar_model_NGas, n.ahead = 1)$pred)
  
  
  best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts[1:(2500+i-1)], solver = 'hybrid')
  residuals_oil[[i]] <- residuals(best_model_oil_fit)
  forecast_oil <- ugarchforecast(best_model_oil_fit, n.ahead = 1)
  fitted_forecast_oil[[i]] <- fitted(forecast_oil)
  sigma_forecast_oil[[i]] <- sigma(forecast_oil)
  ar_model_oil <- ar.ols(oil_ts[1:(2500+i-1)],order = 1)
  ar_forecast_oil[[i]] <- as.numeric(predict(ar_model_oil, n.ahead = 1)$pred)
  
  best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts[1:(2500+i-1)], solver = 'hybrid')
  residuals_coal[[i]] <- residuals(best_model_coal_fit)
  forecast_coal <- ugarchforecast(best_model_coal_fit, n.ahead = 1)
  fitted_forecast_coal[[i]] <- fitted(forecast_coal)
  sigma_forecast_coal[[i]] <- sigma(forecast_coal)
  ar_model_coal <- ar.ols(coal_ts[1:(2500+i-1)],order = 1)
  ar_forecast_coal[[i]] <- as.numeric(predict(ar_model_coal, n.ahead = 1)$pred)
  
  CRPS_NGas[[i]] <- GaussCrps(fitted_forecast_NGas[[i]], sigma_forecast_NGas[[i]], NGas_ts[2500+i])
  CRPS_NGas_ar[[i]] <- GaussCrps(ar_forecast_NGas[[i]], mean(NGas_ts[1:(2500+i-1)]), NGas_ts[2500+i])
  CRPS_oil[[i]] <- GaussCrps(fitted_forecast_oil[[i]], sigma_forecast_oil[[i]], oil_ts[2500+i])
  CRPS_oil_ar[[i]] <- GaussCrps(ar_forecast_oil[[i]], mean(oil_ts[1:(2500+i-1)]), sd(oil_ts[1:(2500+i-1)]))
  CRPS_coal[[i]] <- GaussCrps(fitted_forecast_coal[[i]], sigma_forecast_coal[[i]], coal_ts[2500+i])
  CRPS_coal_ar[[i]] <- GaussCrps(ar_forecast_coal[[i]], mean(coal_ts[1:(2500+i-1)]), sd(coal_ts[1:(2500+i-1)]))
  #print(paste0('%',(i/(length(NGas_ts) - 2500))*100))
  print(paste0('%',i/2))
}
dates <- rep(commodities$date[2501:2700],each = 1000)
sim_forecasts_NGas <- as.matrix(rnorm(1000,fitted_forecast_NGas[[1]], sigma_forecast_NGas[[1]]), ncol = 1)
sim_forecasts_oil <- as.matrix(rnorm(1000,fitted_forecast_oil[[1]], sigma_forecast_oil[[1]]), ncol = 1)
sim_forecasts_coal <- as.matrix(rnorm(1000,fitted_forecast_coal[[1]], sigma_forecast_coal[[1]]), ncol = 1)
for(i in 2:200){
  sim_forecasts_NGas <- rbind(sim_forecasts_NGas, as.matrix(rnorm(1000,fitted_forecast_NGas[[i]], sigma_forecast_NGas[[i]]), ncol = 1))
  sim_forecasts_oil <- rbind(sim_forecasts_oil, as.matrix(rnorm(1000,fitted_forecast_oil[[i]], sigma_forecast_oil[[i]]), ncol = 1))
  sim_forecasts_coal <- rbind(sim_forecasts_coal, as.matrix(rnorm(1000,fitted_forecast_coal[[i]], sigma_forecast_coal[[i]]), ncol = 1))
}
sim_forecasts_NGas_df <- as.data.frame(cbind(dates, sim_forecasts_NGas))
sim_forecasts_NGas_df$dates <- as.Date(sim_forecasts_NGas_df$dates)
sim_forecasts_oil_df <- as.data.frame(cbind(dates, sim_forecasts_oil))
sim_forecasts_oil_df$dates <- as.Date(sim_forecasts_oil_df$dates)
sim_forecasts_coal_df <- as.data.frame(cbind(dates, sim_forecasts_coal))
sim_forecasts_coal_df$dates <- as.Date(sim_forecasts_coal_df$dates)
ggplot(sim_forecasts_NGas_df[95000:130000,], aes(x = V2, y = as.factor(dates))) +
  geom_density_ridges(rel_min_height = 0.005)
ggplot(sim_forecasts_NGas_df[1:30000,], aes(x = V2, y = as.factor(dates))) +
  geom_density_ridges(rel_min_height = 0.005)



#############SECOND PART###################
# part h
best_model_NGas <- ugarchspec(mean.model = list(armaOrder = c(2, 1)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts_train, solver = 'hybrid')
best_model_NGas_fit_std_residuals <- residuals(best_model_NGas_fit)/sigma(best_model_NGas_fit)
#pit_NGas <- pdist("norm",best_model_NGas_fit_std_residuals, mu = 0, sigma = 1, skew = coef(best_model_NGas_fit)["skew"], shape=coef(best_model_NGas_fit)["shape"])
pit_NGas <- pobs(residuals(best_model_NGas_fit))
plot_pit(pit_NGas)

best_model_oil <- ugarchspec(mean.model = list(armaOrder = c(0, 1)),
                             variance.model = list(garchOrder = c(1, 1)))
best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts_train, solver = 'hybrid')
best_model_oil_fit_std_residuals <- residuals(best_model_oil_fit)/sigma(best_model_oil_fit)
#pit_oil <- pdist("norm",best_model_oil_fit_std_residuals, mu = 0, sigma = 1, skew = coef(best_model_oil_fit)["skew"], shape=coef(best_model_oil_fit)["shape"])
pit_oil <- pobs(residuals(best_model_oil_fit))
plot_pit(pit_oil)

best_model_coal <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts_train, solver = 'hybrid')
best_model_coal_fit_std_residuals <- residuals(best_model_coal_fit)/sigma(best_model_coal_fit)
#pit_coal <- pdist("norm",best_model_coal_fit_std_residuals, mu = 0, sigma = 1, skew = coef(best_model_coal_fit)["skew"], shape=coef(best_model_coal_fit)["shape"])
pit_coal <- pobs(residuals(best_model_coal_fit))
plot_pit(pit_coal)

residuals_df <- as.data.frame(cbind(pit_NGas, pit_oil, pit_coal))
colnames(residuals_df) <- c('NGas', 'oil', 'coal')

cor(residuals_df, method = 'pearson')
cor(residuals_df, method = 'kendall')
cor(residuals_df, method = 'spearman')
#NGas vs. oil
ggplot(residuals_df)+geom_point(aes(NGas,oil)) + theme_minimal() + labs(x ='Natural Gas', y = 'Oil')# + geom_smooth(method = 'lm')
ggplot()+geom_point(aes(best_model_NGas_fit_std_residuals, best_model_oil_fit_std_residuals)) + theme_minimal() + labs(x ='Natural Gas', y = 'Oil')# + geom_smooth(method = 'lm')

ggplot(residuals_df)+geom_point(aes(NGas,coal)) + theme_minimal() + labs(x ='Natural Gas', y = 'Coal')# + geom_smooth(method = 'lm')
ggplot()+geom_point(aes(best_model_NGas_fit_std_residuals, best_model_coal_fit_std_residuals)) + theme_minimal() + labs(x ='Natural Gas', y = 'Coal')# + geom_smooth(method = 'lm')

ggplot(residuals_df)+geom_point(aes(coal,oil)) + theme_minimal() + labs(x ='Coal', y = 'Oil')# + geom_smooth(method = 'lm')
ggplot()+geom_point(aes(best_model_coal_fit_std_residuals, best_model_oil_fit_std_residuals)) + theme_minimal() + labs(x ='Coal', y = 'Oil')# + geom_smooth(method = 'lm')

hist(best_model_NGas_fit_std_residuals, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)
hist(best_model_coal_fit_std_residuals, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)
hist(best_model_oil_fit_std_residuals, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)

hist(residuals_df$NGas, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)
hist(residuals_df$coal, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)
hist(residuals_df$oil, breaks = 50, xlab = NULL, ylab = NULL, main = NULL, axes = FALSE)

ggplot(residuals_df,aes(NGas,coal))+geom_point()# + geom_smooth(method = 'lm')
ggplot(residuals_df,aes(coal,oil))+geom_point()# + geom_smooth(method = 'lm')

scatterhist(best_model_NGas_fit_std_residuals, best_model_oil_fit_std_residuals, xlab = 'NGas', ylab = 'Oil')
#part i
normal_fit <- fitCopula(normalCopula(dim = 3, dispstr = 'un'), cbind(as.numeric(pit_coal), as.numeric(pit_NGas), as.numeric(pit_oil)), method = 'ml')
normal_fit_aic <- -2 * normal_fit@loglik + 2 * 3 
t_fit <- fitCopula(tCopula(dim = 3, dispstr = 'un'), cbind(as.numeric(pit_coal), as.numeric(pit_NGas), as.numeric(pit_oil)), method = 'ml')
t_fit_aic <- -2 * t_fit@loglik + 2 * 3
#part j
prob_forecast_NGas_t <- list()
prob_forecast_coal_t <- list()
prob_forecast_oil_t <- list()

prob_forecast_NGas_normal <- list()
prob_forecast_coal_normal <- list()
prob_forecast_oil_normal <- list()

prob_forecast_NGas_indep <- list()
prob_forecast_coal_indep <- list()
prob_forecast_oil_indep <- list()
for(i in 1:200){
  pit_NGas <- pobs(unlist(residuals_NGas[[i]]))
  pit_coal <- pobs(unlist(residuals_coal[[i]]))
  pit_oil <- pobs(unlist(residuals_oil[[i]]))
  
  #independent copula
  indep_copula <- indepCopula(dim = 3)
  pse_obs_indep <- rCopula(1000, indep_copula)
  prob_forecast_NGas_indep[[i]] <- qnorm(pse_obs_indep[,1], mean = fitted_forecast_NGas[[i]], sd = sigma_forecast_NGas[[i]])
  prob_forecast_coal_indep[[i]] <- qnorm(pse_obs_indep[,2], mean = fitted_forecast_coal[[i]], sd = sigma_forecast_coal[[i]])
  prob_forecast_oil_indep[[i]] <- qnorm(pse_obs_indep[,3], mean = fitted_forecast_oil[[i]], sd = sigma_forecast_oil[[i]])
  
  #gaussian copula
  normal_fit <- fitCopula(normalCopula(dim = 3, dispstr = 'un'), cbind(as.numeric(pit_NGas), as.numeric(pit_coal), as.numeric(pit_oil)), method = 'ml')
  normal_fit_est <- normal_fit@estimate
  pse_obs_normal <- rCopula(1000, tCopula(param = normal_fit_est,dim = 3, dispstr = 'un'))
  prob_forecast_NGas_normal[[i]] <- qnorm(pse_obs_normal[,1], mean = fitted_forecast_NGas[[i]], sd = sigma_forecast_NGas[[i]])
  prob_forecast_coal_normal[[i]] <- qnorm(pse_obs_normal[,2], mean = fitted_forecast_coal[[i]], sd = sigma_forecast_coal[[i]])
  prob_forecast_oil_normal[[i]] <- qnorm(pse_obs_normal[,3], mean = fitted_forecast_oil[[i]], sd = sigma_forecast_oil[[i]])
  
  #t-copula
  t_fit <- fitCopula(tCopula(dim = 3, dispstr = 'un'), cbind(as.numeric(pit_NGas), as.numeric(pit_coal), as.numeric(pit_oil)), method = 'ml')
  t_fit_est <- t_fit@estimate[1:3]
  pse_obs_t <- rCopula(1000, tCopula(param = t_fit_est,dim = 3, dispstr = 'un'))
  prob_forecast_NGas_t[[i]] <- qnorm(pse_obs_t[,1], mean = fitted_forecast_NGas[[i]], sd = sigma_forecast_NGas[[i]])
  prob_forecast_coal_t[[i]] <- qnorm(pse_obs_t[,2], mean = fitted_forecast_coal[[i]], sd = sigma_forecast_coal[[i]])
  prob_forecast_oil_t[[i]] <- qnorm(pse_obs_t[,3], mean = fitted_forecast_oil[[i]], sd = sigma_forecast_oil[[i]])
  
  print(paste0('%',i/2))
  }

fitCopula(gumbelCopula(dim = 3), cbind(as.numeric(pit_coal), as.numeric(pit_NGas), as.numeric(pit_oil)), method = 'ml')
rCopula(1000,fitCopula(tCopula(dim = 3), cbind(as.numeric(pit_coal), as.numeric(pit_NGas), as.numeric(pit_oil)), method = 'ml'))


# Create the joint distribution using the independent copula
pse_obs <- qnorm(rCopula(1000, tCopula(param = t_fit_est,dim = 3, dispstr = 'un')))

qnorm(rCopula(1, indepCopula(dim = 3)))

es_sample(as.vector(cbind(as.numeric(NGas_ts[[2501]]),as.numeric(coal_ts[[2501]]),as.numeric(oil_ts[[2501]]))),t(as.matrix(cbind(prob_forecast_NGas_t[[1]], prob_forecast_coal_t[[1]], prob_forecast_oil_t[[1]]))))