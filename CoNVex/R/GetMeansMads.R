# TODO: getMADMeanRatio
# 
# Author: pv1
###############################################################################

`GetMeansMads` <- function(x) {
	sid = as.character(x[7]);
	fi = as.numeric(x[length(x)-1]); ti = as.numeric(x[length(x)]);
	aXp = allX[c(fi:ti),];
	sample_column = which(names(allX)==sid); ## Column of the CNV sample
	mean_logratio_all = as.numeric(colMeans(aXp[,5:length(aXp)])); mean_ratio_sample = mean_logratio_all[sample_column-4];
	mad_mean_ratio = mad(mean_logratio_all); 
	return(c(round(mean_ratio_sample,digits=2),round(mad_mean_ratio,digits=2)));
}