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