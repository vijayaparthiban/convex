# TODO: Get centromere boundaries for each chromosome
# 
# Author: pv1
###############################################################################

`GetCentroBoundaries` <- function(d,centromeres) {
	
	chr = unique(d[,1]); SE = NULL;
	for(c in chr) {	
		uSE = centromeres[as.character(centromeres[,1])==c,];
		dd = d[d[,1]==c & d[,3] < uSE[1,2],]; 
		if(nrow(dd) > 0) { # before the centromere
			SEtmp = data.frame(c,min(dd[,2]),max(dd[,3]),nrow(dd)); SE = rbind(SE,SEtmp);
		}
		dd = d[d[,1]==c & d[,2] > uSE[1,3],]; 
		if(nrow(dd) > 0) { # after the centromere
			SEtmp = data.frame(c,min(dd[,2]),max(dd[,3]),nrow(dd)); SE = rbind(SE,SEtmp);
		}
	}	
	names(SE) = c("chr","start","end","num_probes"); return(SE);
}



