
library(splines)
##periodic-spline:
get.pbas<- function( n , S=24, dK = S/4, ord=4){
	## ord=4 --> cubic splines
	## dK = equidistance distance
	## S= season
	## support will be 1:n
	stp<- dK
	x<- 1:S ## must be sorted!
	lb<- x[1]
	ub<- x[length(x)]
	knots<- seq(lb-(0)*stp,ub+(1)*stp, by=stp)
	derivs<- numeric(length(x))
	## some stuff adjusted from pbc-package to be faster
	degree<- ord-1
	nKnots = length(knots)
	Aknots = c( knots[1] - knots[nKnots] + knots[nKnots - degree:1] , knots,  knots[nKnots] + knots[1:degree + 1] - knots[1] )
	basisInterior <- splineDesign(Aknots, x, ord, derivs) 
	basisInteriorLeft <- basisInterior[, 1:(ord-1), drop = FALSE]
	basisInteriorRight <- basisInterior[, (ncol(basisInterior) - ord+2):ncol(basisInterior), drop = FALSE]
	basis <- cbind(basisInterior[, -c(1:(ord-1), (ncol(basisInterior) - ord+2):ncol(basisInterior)), drop = FALSE], basisInteriorLeft + basisInteriorRight)
	t(array(t(basis), dim= c(dim(basis)[2] ,n)))
}

# ## EXAMPLE hourly data
# S<- 24
# A<- 365.24
# n.total<- 3*S*A
# BAS.daily<- get.pbas(n.total, S=S, dK= S/8 ) 
# BAS.annual<- get.pbas(n.total, S=S*A, dK= S*A/6 )  
# 
# ## example plot...
# ts.plot(BAS.annual, col=rainbow(dim(BAS.annual)[2]))

