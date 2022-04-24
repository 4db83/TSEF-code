# clear screen and workspace
cat("\014"); rm(list = ls()); gc()
# set working directory if need

# INSTALL PACMAN PACKAGE MANAGER IF NOT INSTALLED (Note: may need to disable windows firewall)
if (!"pacman" %in% installed.packages()){install.packages("pacman"); cat("pacman installed\n")}
# MAKE SOURCE DATA DIRECTORY 
# if (!dir.exists(raw_data)){cat("Making director"); dir.create(raw_data)}
# LOAD REQUIRED PACKAGES 
pacman::p_load(polynom, matlab)
source("./R_help_functions.R")

# Lecture Example AR(2): phi_1 = 1.5, phi_2 = -0.56
phi_1 =  1.5 
phi_2 = -0.56

#phi_1 =  1.4
#phi_2 = -0.85
# AR Lag polynomial
phi.L = c(1, -phi_1, -phi_2)
cat("AR Lag Polynomial is: "); print(as.polynomial(phi.L))
# MA Lag polynomial
thet.L = c(1, 0.3, 0.1)
cat("MA Lag Polynomial is: "); print(as.polynomial(thet.L))

# lag polynomial roots of the process are (need to specify 1 - phi_L -phi_L^2 here)
lag.roots = round(polyroot(phi.L),4)
if (Im(lag.roots[1])>1e-15)    {
    cat("Roots of Lag polynomial are:", lag.roots, "\n")
} else {
    cat("Roots of Lag polynomial are:", Re(lag.roots), "\n" )
}

# plot phi.L -------------------------------------------------------------------
par(mfrow=c(1,2), mar = c(3,3,2,0.8),oma = c(1,1.2,1,1), mgp = c(1.7,0.5,0))
plot(polynomial(phi.L), xlim=c(1, 1.5), ylim=c(-0.01,0.01))
abline(h=c(0), lty=c(1), col=c(1))
abline(v=c(lag.roots), lty=c(2), col=c(2))

# compute the factored roots as 
fact.roots = round(polyroot(fliplr(phi.L)),4)
if (Im(fact.roots[1])>1e-15)    {
    cat("Factored Roots are:", fact.roots, "\n")
} else {
    cat("Factored Roots are:", Re(fact.roots), "\n")
}

plot(polynomial(fliplr(phi.L)), xlim=c(0.5, 0.9), ylim=c(-0.01,0.01))
abline(h=c(0), lty=c(1), col=c(1))
abline(v=c(fact.roots), lty=c(2), col=c(2))

# theoretical ACF
plot.acf0(phi.L,thet.L)

# plot of simulated series of sample size 200.
#set.seed(1234);
B = 200; T = 200; # B = is the burn in period
y.ar2 = arima.sim(n = T+B, list(ar = -phi.L[-1], ma = thet.L[-1]),sd = 1); 
y.ar2 = y.ar2[-(1:B)]; # kill the burn-in period
plot(y.ar2,type='l')

# plot the sample acf and pacf
plot.acf(y.ar2)

# Psi inf and Pi inf from Brockwell and Davies page 87.
N.terms = 10
phi.L   = cbind(1,-.5);
thet.L  = cbind(1,.4);
ma.out = arma2ma(phi.L, thet.L, N.terms)
ar.out = arma2ar(phi.L, thet.L, N.terms)
round(cbind(ma.out,ar.out),4)




