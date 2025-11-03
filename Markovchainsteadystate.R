# Fit SPI or SPEI for historical Data
# Perform Markov Chain Analysis (stationary and aperiodic)
# Original Code and newspi implementation by Dr. Venki Uddammeri
# Adapted by Janie Davis for downscaled Projections
# Runs by Taha

library(SPEI)
library(markovchain)
library(lubridate)
library(zoo)
library(fitdistrplus)
library(actuar)

path = '/home/taha/Documents/Fluxes'
fname = 'fluxes_downscaled_Chist.csv'           # 'fluxes_downscaled_Chist.csv' # fluxes_downscaled_Uhist.csv # fluxoptimfinal_VC.csv
acc = 3 # monthly accumulation for standardized index
setwd(path)
a <- read.csv(fname)
head(a)

fllogis <- function(x)
{
  x <- as.vector(na.omit(x))
  mu <- min(x) - 200 # arbitrary centring
  xz <- x - mu
  fit <- fitdist(xz,distr='llogis',start=list(shape=1,scale=1))
  xx <- list('loc'=mu,'fit'=fit)
  return(xx)
}

# For spi of P - PET and avoiding infinity/ nan values through log logistic, 
# primarily for future scenarios, but to keep methodology consistent with copula analysis
spinew <- function(x, scale, ref.start = NULL, ref.end = NULL, na.rm = FALSE)
{ 
  xs <- rollsum(x,scale,align='right',fill=NA)
  xts <- ts(xs,start= ref.start,frequency=12)
  xfit <- window(xts,ref.start,ref.end)
  xmat <- matrix(xfit,ncol=12,byrow=TRUE)
  zza <- apply(xmat,2,fllogis)  # Fit functions for each month
  locx <- c()
  shapex <- c()
  scalex <- c()
  for(i in seq(1,12,1)){
    locx <- c(locx,zza[[i]]$loc)
    shapex <- c(shapex,zza[[i]]$fit$estimate[1])
    scalex <- c(scalex,zza[[i]]$fit$estimate[2])
  }
  parms <- data.frame(locx,shapex,scalex)
  xm <- matrix(xts,ncol=12,byrow=TRUE)
  ny <- floor(length(xts) / 12) # number of years
  spix <- mat.or.vec(ny,12) # 205 130
  for(i in seq(1,12,1)){
    xma <- xm[,i] - parms[i,1]
    qut <- pllogis(xma,shape=parms[i,2],scale=parms[i,3])
    qut <- ifelse(qut==0,1E-10,qut)
    qut <- ifelse(qut==1,(1-1E-10),qut)
    spix[,i] <- qnorm(qut)
  }
  spixts <- ts(as.vector(t(spix)),start=ref.start,frequency=12)
  return(spixts)
}



# Get start and end times
a$Date <- as.Date(paste(a$Year, a$Month, "01", sep="-"))
bgn <- a$Date[1] # first date in the input CSV
endx <- as.Date("2014-12-31") # last date in the input CSV

# Subset 'a' to include only rows up to 'endx'
a <- a[a$Date <= endx, ]

spidate <- bgn %m+% months(acc-1) # this will make the SPI start after enough months have passed
strtyr <- year(spidate) # the start year
strtmo <- month(spidate) # the start month
# the sample CSV has one column for all dates
# mine has a year column and a month column so will need to adjust this part

pbcsd <- a$W # - a$PET # pull the column of values you want to use # change this


histts <- ts(pbcsd,start=c(1895,1),frequency=12)  # monthly Time series
SPI <- spi(histts,acc,ref.start=c(1895,1),ref.end=c(2014,12)) # spi for historical through SPEI package
# spinew for P - PET
summary(SPI) 
#spicoef <- SPI$coefficients  # store parameters of the gamma fit

# Write gamma parameters to CSV
#coeff_2d <- spicoef[, 1, ]
#coeff_t <- t(coeff_2d)
#spi_df <- as.data.frame(coeff_t)
#spi_df$Month <- month.abb
#spi_df <- spi_df[, c("Month", "alpha", "beta")]
#write.csv(spi_df, "gammacoeff_Uhist_W.csv", row.names = FALSE) # change name with scenario


spifit <- na.omit(SPI$fitted) # remove NA values change to just SPI for spi
spits <- ts(spifit,start=c(strtyr,strtmo),frequency=12) # create monthly time series
plot(spits,xlab='Date',ylab='SPI3',col='blue',main='3 Month SPI BMT Downscaled') # plot monthly time series
grid()
abline(h=0,lty=2,lwd=2) # put reference line at 0

# Create a 3 state series and create stationary Markov Chain
cutoff <- c(-Inf,-0.5,0.5,Inf) # normal is from -0.5 to 0.5; dry is less, wet is more
labels <- c('Dry','Normal','Wet')
spistate <- cut(spits,breaks=cutoff,labels=labels) # separate the indexed values into the three categories
spistate

# Create a Transition Matrix
spimc <- markovchainFit(data = spistate)
spitm <- spimc$estimate

# Compute Markovian Statistics
# steady state
spiss <- steadyStates(spitm) # probability
spiss
# Mean recurrance time
spimrt <- 1/spiss # the units are in months
spimrt
# Mean First Passage Time
mfpt <- meanFirstPassageTime(spitm)  # in months
mfpt

# Absorbing states - Check
absorbingStates(spitm)  # should not be because there is no 0 in transition matrix

# Create time-series and write to a dataframe
spistatets <- ts(spistate,start=c(strtyr,strtmo),frequency=12) # time series of state values
datex <- seq.Date(spidate,endx,by='month')

zz <- data.frame(datex,spits,spistatets)

out = '/home/taha/Documents/MarkovChains/Steady'
setwd(out)
write.csv(zz,'markstates_CNRM_SSMI.csv',row.names=FALSE) # change name for scenario and index