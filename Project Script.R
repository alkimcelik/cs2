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

# Create a time series object from your data
# Replace 'your_data' with your actual time series data
ts_data <- ts(your_data)

# Specify the ARMA(p,q) and GARCH(r,s) orders
p <- 1  # Autoregressive order
q <- 1  # Moving average order
r <- 1  # GARCH order
s <- 1  # GARCH order for lagged conditional variances

# Specify the specification for the ARMA-GARCH model
spec <- ugarchspec(mean.model = list(armaOrder = c(p, q)),
                   variance.model = list(garchOrder = c(r, s)))

# Fit the ARMA-GARCH model to the data
fit <- ugarchfit(spec, ts_data)

# Print the model summary
show(fit)

# Access the estimated parameters
estimated_params <- coef(fit)
print(estimated_params)