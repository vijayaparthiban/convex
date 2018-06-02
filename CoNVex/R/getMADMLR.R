# TODO: get regional Mean Log Ratio, Mean EW score
# 
# Author: pv1
###############################################################################

getMADMLR <- function(Fd_rep,aX,aXscore) {
	
	mr_mads = NULL;
	for(i in 1:nrow(Fd_rep)) {
		chr = as.character(Fd_rep[i,1]); start = Fd_rep[i,2]; end = Fd_rep[i,3]; 
		aXp = aX[as.character(aX[,2])==as.character(chr) & (aX[,3]>=start) & (aX[,4]<=end),];
		aXsp = aXscore[aXp[,1],];
		sample_column = which(names(aX) == as.character(Fd_rep[i,7]))-7; ## Column of the CNV sample
		mean_logratio_all = as.numeric(colMeans(aXp[,8:length(aXp)])); mean_ratio_sample = mean_logratio_all[sample_column];
		mean_adm_all = as.numeric(colMeans(aXsp[,8:length(aXsp)]));
		mad_mean_ratio = mad(mean_logratio_all); mad_mean_admscore = mad(mean_adm_all);
		mr_mads = rbind(mr_mads,c(mean_ratio_sample,mad_mean_ratio,mad_mean_admscore));
	}
	return(mr_mads);
}

