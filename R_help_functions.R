# R helper / add on functions.
# add these to the top of the script.
######################################################################################

# ARMA to AR(infinity)
arma2ar = function(ar.L,ma.L ,max.lag){
	#if (is.null(max.lag)) max.lag=10
	ar.inf = ARMAtoMA(ma=ar.L[-1],ar=-ma.L[-1],max.lag)
	ar.inf.out = rbind(1,as.matrix(ar.inf))
	return(ar.inf.out) 
}

# ARMA to MA(infinity)
arma2ma = function(ar.L,ma.L ,max.lag){
	#if (is.null(max.lag)) max.lag=10
	ma.inf = ARMAtoMA(ar=-ar.L[-1],ma=ma.L[-1],max.lag)
	ma.inf.out = rbind(1,as.matrix(ma.inf))
	return(ma.inf.out) 
}

# plot theoretical ACF/PACF
plot.acf0 = function(ARpolynomial, MApolynomial, max.lag = 50){
	# if (is.null(max.lag)) max.lag = 51
	ARterms = -ARpolynomial[-1]
	MAterms =  MApolynomial[-1]
	# add one to it to get the required lag structure because of the zapsmall function
	max.lag = max.lag+1
	ACF 	= ARMAacf(ARterms,MAterms,(max.lag),pacf=FALSE)[2:max.lag]
	PACF 	= zapsmall(ARMAacf(ARterms,MAterms,max.lag,pacf=TRUE))[1:(max.lag-1)]
	LAG 	= 1:(max.lag-1)
	minA 	= min(ACF)
	minP	= min(PACF)
	LW 		= 6.6
	colr 	= "skyblue2"
	minu 	= min(minA,minP)-.01
	old.par <- par(no.readonly = TRUE)
	par(mfrow=c(1,2), mar = c(3,3,2,0.8),oma = c(1,1.2,1,1), mgp = c(1.7,0.5,0))
	plot(LAG, ACF, type="h",ylim=c(minu,1) ,lwd=LW, las=1, lend = 1, col=colr, xlim=c(1, max.lag),
	     xlab="Lag")
	#main=paste("Series: ",deparse(substitute(series))))
	abline(h=c(0), lty=c(1), col=c(1))
	plot(LAG, PACF, type="h",ylim=c(minu,1) ,lwd=LW, las=1, lend = 1, col=colr, xlim=c(1, max.lag),
	     xlab="Lag")
	abline(h=c(0), lty=c(1), col=c(1))
	on.exit(par(old.par))
	ACF<-round(ACF,2); PACF<-round(PACF,2)    
	return(invisible(cbind(ACF, PACF))) 
}

# plot sample ACF/PACF
plot.acf=function(data_in, max.lag = 50){
	num=length(data_in)
	if (max.lag > (num-1)) stop("Number of lags exceeds number of observations")
	ACF 	= acf(data_in, max.lag, plot=FALSE)$acf[-1]
	PACF 	= pacf(data_in, max.lag, plot=FALSE)$acf
	#LAG=1:max.lag/frequency(data_in)
	LAG 	= 1:max.lag
	minA 	= min(ACF)
	minP 	= min(PACF)
	U 		= 2/sqrt(num)
	L 		=-U
	LW 		= 6.6
	colr 	= "skyblue2"
	minu 	= min(minA,minP,L)-.01
	old.par <- par(no.readonly = TRUE)
	par(mfrow=c(1,2), mar = c(3,3,2,0.8),oma = c(1,1.2,1,1), mgp = c(1.7,0.5,0))
	plot(LAG, ACF, type="h",ylim=c(minu,1) ,lwd=LW, las=1, lend = 1, col=colr, xlim=c(1, max.lag))
	  #main=paste("Series: ",deparse(substitute(series))))
	  abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,2,2))
	plot(LAG, PACF, type="h",ylim=c(minu,1) ,lwd=LW, lend = 1, col=colr, xlim=c(1, max.lag))
	  abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,2,2))
	on.exit(par(old.par))  
	ACF<-round(ACF,2); PACF<-round(PACF,2)    
	return(invisible(cbind(ACF, PACF))) 
}
