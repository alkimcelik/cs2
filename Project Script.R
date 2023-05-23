calc_AIC <- function(model){
  log_likelihood <- model@fit$LLH
  
  # Calculate the number of parameters
  num_parameters <- length(model@fit$coef)
  
  # Calculate the AIC
  aic <- -2 * log_likelihood + 2 * num_parameters
  return(aic)
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
install.packages("rugarch")
library(rugarch)

### Train datalarini cikardik
NGas_ts <- ts(commodities$NGas, start = c(2010, 03, 16), end = c(2020, 10, 27), freq = 365)
oil_ts <- ts(commodities$oil, start = c(2010, 03, 16), end = c(2020, 10, 27), freq = 365)
coal_ts <- ts(commodities$coal, start = c(2010, 03, 16), end = c(2020, 10, 27), freq = 365)
NGas_ts_train <- as.ts(NGas_ts[1:2500], start = c(2010, 03, 16), end = c(2020, 01, 17),freq = 365)
oil_ts_train <- as.ts(oil_ts[1:2500], start = c(2010, 03, 16), end = c(2020, 01, 17),freq = 365)
coal_ts_train <- as.ts(coal_ts[1:2500], start = c(2010, 03, 16), end = c(2020, 01, 17),freq = 365)
NGas_ts_test <- as.ts(NGas_ts[1:2500], start = c(2020, 01, 18), end = c(2020, 10, 27),freq = 365)
oil_ts_test <- as.ts(oil_ts[1:2500], start = c(2020, 01, 18), end = c(2020, 10, 27),freq = 365)
coal_ts_test <- as.ts(coal_ts[1:2500], start = c(2020, 01, 18), end = c(2020, 10, 27),freq = 365)

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
##Best model with parameters p=2, q=1, r=1, s=1
best_model_NGas <- ugarchspec(mean.model = list(armaOrder = c(2, 1)),
                                           variance.model = list(garchOrder = c(1, 1)))
best_model_NGas_fit <- ugarchfit(best_model_NGas, NGas_ts_train, solver = 'hybrid')

####Coal Model##########
# Specify the specification for the ARMA-GARCH model
AIC_values <- list()
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
          fit <- ugarchfit(spec, coal_ts_train, solver = 'hybrid')
          AIC_values[[counter]] <- calc_AIC(fit)#Akaike
          index_list[[counter]] <- paste(i,j,k,l)
          counter <- counter + 1
          print(paste(i,j,k,l))
        }
      }
    }
  }
}
AIC_values_index <- as.data.frame(cbind(unlist(AIC_values), unlist(index_list)))
best_model_coal <- ugarchspec(mean.model = list(armaOrder = c(i, j)),
                              variance.model = list(garchOrder = c(k, l)))
best_model_coal_fit <- ugarchfit(best_model_coal, coal_ts_train, solver = 'hybrid')

####oil Model##########
# Specify the specification for the ARMA-GARCH model
AIC_values <- list()
index_list <- list()
counter <- 1
for(i in p){
  for(j in q){
    for(k in r){
      for(l in s){
          spec <- ugarchspec(mean.model = list(armaOrder = c(i, j)),
                             variance.model = list(garchOrder = c(k, l)))
          fit <- ugarchfit(spec, oil_ts_train, solver = 'hybrid')
          AIC_values[[counter]] <- calc_AIC(fit)#Akaike
          index_list[[counter]] <- paste(i,j,k,l)
          counter <- counter + 1
          print(paste(i,j,k,l))
      }
    }
  }
}

best_model_oil <- ugarchspec(mean.model = list(armaOrder = c(i, j)),
                              variance.model = list(garchOrder = c(k, l)))
best_model_oil_fit <- ugarchfit(best_model_oil, oil_ts_train, solver = 'hybrid')

# Access the estimated parameters
estimated_params <- coef(fit)
print(estimated_params)

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