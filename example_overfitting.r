# Script: example_overfitting.R ----
# Simulate WN series and try to fit up to ARMA(Pmax,Qmax) using auto.arima function from forecast package
# NOTE: To run this script, you also need "R_help_functions.R" available at: https://github.com/4db83/code-TSEF.
# clear screen/workspace
cat("\014"); rm(list = ls()); gc()
# SET WORKING DIRECTORY PATH IF NEED
# this is mine, your need to set your path
# setwd("D:/_teaching/_current.teaching/_SU.TSEF/lectures/code-TSEF")

# INSTALL PACMAN PACKAGE MANAGER IF NOT INSTALLED. 
# (Note: you may need to disable windows firewall to allow installation)
if (!"pacman" %in% installed.packages()){install.packages("pacman"); cat("pacman installed\n")}
# LOAD REQUIRED PACKAGES
pacman::p_load(polynom, matlab, forecast)
source("./R_help_functions.R")

Pmax  = 5;
Qmax  = 5;
T     = 3e2;
a0    = rep(0,Pmax);
b0    = rep(0,Qmax);
T0    = proc.time();
Nsim  = 10;

# space allocation;
arma.terms = matrix(0,Nsim,2)    
# set simulation seed
set.seed(123)
dx = matrix( rnorm(T*Nsim), T, Nsim )

for (i in 1:Nsim){
    aout = auto.arima( dx[,i],
        d=0, D=0, max.Q=0, max.P=0, 
        max.p=Pmax,
        max.q=Qmax,
        start.p=a0,
        start.q=b0,
        seasonal=FALSE,
        trace=FALSE,
        stepwise=FALSE,
        ic=c("aic"), max.order=20 )
    # store the ARMA terms
    arma.terms[i,] = aout$arma[1:2]
    # print iterations to screen
    cat("Iteration number =", (i), " ARMA terms: ", aout$arma[1:2], "\n")
    t.tmp = Nsim*(proc.time() - T0)
    if (i == 1)
    cat("Total time is approximately: ", t.tmp[3], "Seconds\n")
}
# tT = proc.time() - T0
cat("Elapased time is:", (proc.time() - T0)[3], "\n" )
print(arma.terms)

