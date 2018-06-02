`EWScore` <- function(d,doNorm=0) {

	w = 1/(d[,5]^2); # Weights based on #reads and regional MAD
	
	if(doNorm == 0) { # Don't normalize log2 ratio
		EWS = d[,4]*w/sqrt(w);
		d = data.frame(d,d[,4],w,EWS);
	} else { # Normalize
		dd = (d[,4]-mean(d[,4]))/var(d[,4]); # Normalize
		EWS = dd*w/sqrt(w);
		d = data.frame(d,dd,w,EWS);
	}	
	return(d);
}