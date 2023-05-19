# CLEAR SCREEN/WORKSPACE
cat("\014"); rm(list = ls()); gc()
# SET DEFAULTS: display options, font and y axis label rotation
options(digits = 12); options(scipen = 999);  options(max.print=10000)
windowsFonts("Palatino" = windowsFont("Palatino Linotype")); par(las = 1, family = "Palatino")
# INSTALL PACMAN PACKAGE MANAGER IF NOT INSTALLED (Note: may need to disable windows firewall for packages to install)
if (!"pacman" %in% installed.packages()){install.packages("pacman"); cat("pacman installed\n")}
# LOAD REQUIRED PACKAGES
# source("./R_help_functions.R")
# source("https://raw.github.com/tonybreyal/Blog-Reference-Functions/master/R/bingSearchXScraper/bingSearchXScraper.R")
# ls()
#> [1] "bingSearchXScraper"
#> 
pacman::p_load(forecast, matlab, tidyverse, readxl, arrow, curl, mFilter) 
cat("\014")
set.seed(123)

# A. Computer Exercises 1) ------
N = 1e4
T = 170 + 1
# generate (T+1)*(N) matrix of standard normal random numbers
matNorm = matrix(rnorm(T*N,0,1),nrow = T,ncol = N)
# make pure RW process
matRW   = apply(matNorm, 2, cumsum)
# make time trend variable
tt = 1:(T-1)

# make storage space for coefficient estimates as well as tstats
rho0   = zeros(N,1); rhoC   = zeros(N,1); rhoCT   = zeros(N,1)
tstat0 = zeros(N,1); tstatC = zeros(N,1); tstatCT = zeros(N,1)

tic()
# DF regressions
for (ii in 1:N) {
  y1 = matRW[1:T-1,ii]
  y  = matRW[2:T,  ii]
  # run the regressions and store rho_hat as well as tstat = (rho-1)/se(rho)
  df0 = lm(y ~ 0 + y1)   
  rho0[ii]   = df0$coefficients[1]
  tstat0[ii] = (rho0[ii] - 1)/sqrt(diag(vcov(df0)))
  
  dfC = lm(y ~ y1)       
  rhoC[ii] = dfC$coefficients[2]
  tstatC[ii] = (rhoC[ii] - 1)/sqrt(diag(vcov(dfC)))[2]
  
  dfCT = lm(y ~ y1 + tt)  
  rhoCT[ii] = dfCT$coefficients[2]
  tstatCT[ii] = (rhoCT[ii] - 1)/sqrt(diag(vcov(dfCT)))[2]
}
toc()

# density plots ----
ug = seq(-4, 4, length=1e3)
KS0 <- density(rho0)
plot(ug,dnorm(ug), type = "l", lwd = 2, las = 1, ylab="", xlab="Distribution of the t-statistics",
     xlim=c(-6, 4), ylim=c(0, .6), col='black', 
     cex.lab=1.3, cex.axis=1.3
     )
lines(density(tstat0),  col='dodgerblue3', lwd=2)
lines(density(tstatC),  col='brown2'     , lwd=2)
lines(density(tstatCT), col='orange'     , lwd=2)
abline(h=0)
# Add a legend
legend(-6.1, .6, legend = c("No-constant and no trend", "with constant", "with constant and trend", "N(0,1)"), 
       col = c("dodgerblue3", "brown2", "orange", "black" ),
       text.font=1, lwd=2, 
       lty=1, cex=1)
cat("\014")

# %% COMPUTE PERCENTILS OF CRITICAL VALUES
pctls = c(1, 2.5, 5, 7.5, 10, 50, 90, 92.5, 95, 97.5, 99)/100 
dgts  = 4
DF = rbind( qnorm(pctls), 
            quantile(tstat0,  pctls),
            quantile(tstatC,  pctls),
            quantile(tstatCT, pctls)
)
round(DF,digits = 4)

# A. Computer Exercises 2) ----
# get the US data from GitHub # download.file(.) seems to create problems reading the file
# temp_file_p = tempfile(); curl_download(url = 'https://github.com/4db83/TSEF-code/raw/main/data/real_gdp_US_2022Q4.parquet',destfile = temp_file_p)
temp_file_  = tempfile(); curl_download(url = 'https://github.com/4db83/TSEF-code/raw/main/data/real_gdp_US_2022Q4.xlsx', destfile = temp_file_)
cat("\014")
# read US data
# this below uses the parquet file system as an alternative
# usdata      = read_parquet(temp_file_) # this uses the parquet file system as an alternative
# usdata$Time = as.Date(usdata$Time); colnames(usdata)[1] = "Date"

usdata = read_excel(temp_file_)
# make y and dy from US Data
usdata$dates = as.Date(usdata$dates, "%d-%b-%Y") 
# rename to dates to Date
colnames(usdata)[1] = "Date"
usdata$y  = log(usdata$gdpc1)
usdata$dy = 100*(usdata$y - lag(usdata$y,1))
## TRUNCATE THE SAMPLE TO 'Q1-1947', 'Q4-2019' ----
cat("\014")
usdata = usdata[usdata$Date >= "1947-01-01" & usdata$Date <= "2019-12-31", ]
tail(usdata)

# ESTIMATE THE VARIOUS ARMA MODELS
# NOTE: in R, the unconditional mean instead of the intercept is reported in the estimates
arma11 = Arima(usdata$dy, order=c(1,0,1))
arma01 = Arima(usdata$dy, order=c(0,0,1))
arma10 = Arima(usdata$dy, order=c(1,0,0))

# BN-DECOMPOSITIONS BASED ON THE LECTURE NOTES EXAMPLES ------

# BN Based on an AR(1)
mu_AR1      = arma10$coef[2] 
trnd_BN_AR1 = usdata$y + arma10$coef[1]/(1 - arma10$coef[1])*(usdata$dy - mu_AR1)
cycl_BN_AR1 = usdata$y - trnd_BN_AR1

# BN Based on an MA(1)
cycl_BN_MA1 =  -arma01$coef[1]*arma01$residuals

# BN Based on an ARMA(1)
mu_ARMA11     = arma11$coef[3];
kk = arma11$coef[1]*(usdata$dy - mu_ARMA11) + arma11$coef[2]*arma11$residuals
trnd_BN_ARMA11 = usdata$y + kk/(1 - arma11$coef[1])
cycl_BN_ARMA11 = usdata$y - trnd_BN_ARMA11

# HP-FILTER
hp = hpfilter(usdata$y,1600)

# join in a  tibble and save to parquet file 
cycles = as_tibble(usdata$Date)
cycles$BN_AR1     = cycl_BN_AR1
cycles$BN_MA1     = cycl_BN_MA1
cycles$BN_ARMA11  = cycl_BN_ARMA11
cycles$HP         = 100*hp$cycle
# head(cycles)
# write_parquet(cycles,"R_cycles.parquet")

# PLOTTING ---------
plot(usdata$Date,  100*hp$cycle, type = "l", lwd = 2, las = 1, ylab="", xlab="",
     ylim=c(-6, 4.5), col='dodgerblue3', tck=-0.015, 
     cex.lab=1.3, cex.axis=1.3
)
lines(usdata$Date, cycl_BN_AR1,    col='brown2' , lwd=2)
lines(usdata$Date, cycl_BN_MA1,    col='orange' , lwd=2)
lines(usdata$Date, cycl_BN_ARMA11, col='green4' , lwd=2)
abline(h=0)
# Add a legend
legend(10700, 4.8, 
       legend = c("HP-Filter Cycle", 
                  "BN-AR(1) Cycle", 
                  "BN-MA(1) Cycle",
                  "BN-ARMA(1,1) Cycle"
                  ), 
       col = c("dodgerblue3", "brown2", "orange","green4" ),
       text.font=1, lwd=2, 
       lty=1, cex=1)

cat("\014")

# print the DF critical values table here 
round(DF,digits = 4)






























