# Script: example_acf_pacf_ar2.R
# NOTE: To be able to run this code, you need the R_help_functions.R available from my Github page.
# clear screen and workspace
cat("\014"); rm(list = ls()); gc()
# SET WORKING DIRECTORY PATH IF NEED
# this is mine, your need to set your path
# setwd("D:/_teaching/_current.teaching/_SU.TSEF/lectures/code-TSEF")

# INSTALL PACMAN PACKAGE MANAGER IF NOT INSTALLED. 
# (Note: you may need to disable windows firewall to allow installation)
if (!"pacman" %in% installed.packages()){install.packages("pacman"); cat("pacman installed\n")}
# LOAD REQUIRED PACKAGES
pacman::p_load(polynom, matlab)
source("./R_help_functions.R")

# AR Lag polynomial
a1 =  1.50 ; a2 = -0.56 
aL = c(1, -a1, -a2)
# companion matrix Phi
Phi = matrix( c(a1,a2,1,0), nrow = 2,ncol = 2)

# lag-polynomial and roots
cat("Lag Polynomial: ", gsub("x","L",as.polynomial(aL)), "\n")
lag.roots = round(polyroot(aL),4)
if (sum(Im(lag.roots)==0)) {
  cat("Roots of Lag polynomial are:", Re(lag.roots), "\n") } else {
	cat("Roots of Lag polynomial are:", lag.roots, "\n" ) 
}

# factored-polynomial and roots
cat("Factored Polynomial: ", gsub("x","\u03BB",as.polynomial(fliplr(aL))), '\n')
fact.roots = round(polyroot(fliplr(aL)),4)
if (sum(Im(fact.roots)==0)) {
  cat("Roots of Lag polynomial are:", Re(fact.roots), "\n") } else {
  cat("Roots of Lag polynomial are:", fact.roots, "\n" ) 
}

# companion matrix Phi and its roots/eigenvalues
char.roots = eigen(Phi)[1] # returns a list
cat("Eigenvalues of Phi are:", unlist(char.roots), "\n")

# plot theoretical ACF/PACF
plot.acf0(aL,1,50)

