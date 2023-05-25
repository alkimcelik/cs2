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
  score <- sigma*(1/sqrt(pi)-2*dnorm(standart_value)-standart_value*(2*pnorm(standart_value)-1))
  return(score)
}
commodities <- read.csv('C:/Users/alkim/OneDrive/Documents/GitHub/cs2/commodities.csv')
commodities$date <- as.Date(commodities$date)
####Plotting the time series basic####
plot(commodities$coal, type = 'l')
plot(commodities$oil, type = 'l')
plot(commodities$NGas, type = 'l')
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

#part b

# Create the xts object
xts_obj <- xts(commodities$NGas, order.by = commodities$date)
### Train datalarini cikardik
NGas_ts <- xts(commodities$NGas, order.by = commodities$date)
oil_ts <- xts(commodities$oil, order.by = commodities$date)
coal_ts <- xts(commodities$coal, order.by = commodities$date)
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
##Best model with parameters p=2, q=1, r=1, s=1
best_model_NGas <- ugarchspec(mean.model = list(armaOrder = c(2, 1)),
                                           variance.model = list(garchOrder = c(1, 1)))
best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts_train, solver = 'hybrid')

best_model_oil <- ugarchspec(mean.model = list(armaOrder = c(0, 1)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts_train, solver = 'hybrid')

best_model_coal <- ugarchspec(mean.model = list(armaOrder = c(2, 2)),
                              variance.model = list(garchOrder = c(1, 1)))
best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts_train, solver = 'hybrid')


# Access the estimated parameters
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
#part d

######ONE SHOT FORECASTING################
forecast_NGas <- ugarchforecast(best_model_NGas_fit, n.ahead = 200)
fitted_forecast_NGas <- fitted(forecast_NGas)
sigma_forecast_NGas <- sigma(forecast_NGas)

forecast_oil <- ugarchforecast(best_model_oil_fit, n.ahead = 200)
fitted_forecast_oil <- fitted(forecast_oil)
sigma_forecast_oil <- sigma(forecast_oil)


forecast_coal <- ugarchforecast(best_model_coal_fit, n.ahead = 200)
fitted_forecast_coal <- fitted(forecast_coal)
sigma_forecast_coal <- sigma(forecast_coal)
CRPS_NGas <- list()
CRPS_oil <- list()
CRPS_coal <- list()
for(i in 1:200){
  CRPS_NGas[[i]] <- CRPS(forecast_NGas@forecast$seriesFor[i], coef(best_model_NGas_fit)[1], sigma(best_model_NGas_fit)[1])
  CRPS_oil[[i]] <- CRPS(forecast_oil@forecast$seriesFor[i], coef(best_model_oil_fit)[1], sigma(best_model_oil_fit)[1])
  CRPS_coal[[i]] <- CRPS(forecast_coal@forecast$seriesFor[i], coef(best_model_coal_fit)[1], sigma(best_model_coal_fit)[1])
}

#####ITERATIVE FORECASTING WITH REESTIMATION######
fitted_forecast_NGas <- list()
mu_forecast_NGas <- list()
sigma_forecast_NGas <- list()
fitted_forecast_oil <- list()
mu_forecast_oil <- list()
sigma_forecast_oil <- list()
fitted_forecast_coal <- list()
mu_forecast_coal <- list()
sigma_forecast_coal <- list()

#for (i in 1:(length(NGas_ts) - 2500)){
for (i in 1:200){
  best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts[1:(2500+i-1)], solver = 'hybrid')
  forecast_NGas <- ugarchforecast(best_model_NGas_fit, n.ahead = 1)
  fitted_forecast_NGas[[i]] <- fitted(forecast_NGas)
  sigma_forecast_NGas[[i]] <- sigma(forecast_NGas)
  
  best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts[1:(2500+i-1)], solver = 'hybrid')
  forecast_oil <- ugarchforecast(best_model_oil_fit, n.ahead = 1)
  fitted_forecast_oil[[i]] <- fitted(forecast_oil)
  sigma_forecast_oil[[i]] <- sigma(forecast_oil)
  
  best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts[1:(2500+i-1)], solver = 'hybrid')
  forecast_coal <- ugarchforecast(best_model_coal_fit, n.ahead = 1)
  fitted_forecast_coal[[i]] <- fitted(forecast_coal)
  sigma_forecast_coal[[i]] <- sigma(forecast_coal)
  
  CRPS_NGas[[i]] <- CRPS(forecast_NGas@forecast$seriesFor, coef(best_model_NGas_fit)[1], sd(NGas_ts[1:(2500+i-1)]))
  CRPS_oil[[i]] <- CRPS(forecast_oil@forecast$seriesFor, coef(best_model_oil_fit)[1], sd(oil_ts[1:(2500+i-1)]))
  CRPS_coal[[i]] <- CRPS(forecast_coal@forecast$seriesFor, coef(best_model_coal_fit)[1], sd(coal_ts[1:(2500+i-1)]))
  #print(paste0('%',(i/(length(NGas_ts) - 2500))*100))
  print(paste0('%',i/2))
}



#CRPS
# Create a vector of observed values
observed <- c(1.2, 2.3, 0.8, 1.5, 3.1)

# Create a matrix of forecast probabilities
# Each row represents a different forecast, and each column represents a probability value
# Replace 'forecast_probs' with your actual forecast probabilities
forecast_probs <- matrix(c(0.1, 0.3, 0.4, 0.2,
                           0.2, 0.2, 0.3, 0.3,
                           0.3, 0.1, 0.2, 0.4,
                           0.4, 0.3, 0.2, 0.1,
                           0.1, 0.1, 0.4, 0.4), nrow = 5, byrow = TRUE)

# Calculate the CRPS for each forecast
crps <- vector("numeric", length = nrow(forecast_probs))
for (i in 1:nrow(forecast_probs)) {
  cum_probs <- cumsum(forecast_probs[i, ])
  crps[i] <- sum((cum_probs - (observed[i] <= cum_probs))^2)
}

# Calculate the average CRPS
average_crps <- mean(crps)

# Print the CRPS values
print(crps)
print(average_crps)