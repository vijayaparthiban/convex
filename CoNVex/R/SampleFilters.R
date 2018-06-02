# TODO: CNV filters
# 
# Author: pv1
###############################################################################

SampleFilters <- function(SamplesAll, sd_dist=3, EO_MAD_cutoff=0.2, type=1) {
	
	# N = sd_dist
	# Expected Columns for SamplesAll
	# Sample IDs with col name sample_id
	# SampleMeans
	# SampleMads
	## Type 1 -- SampleMeans vs MADs filter - Remove samples that are N*sd away from the Lowess and return the rest
	if(type==1) { 
		SamplesAll = SamplesAll[SamplesAll$SampleMads<=EO_MAD_cutoff,]
		## Type 1 -- Return samples that are mean +/- N*sd times from the Lowess
		SamplesAllSorted = SamplesAll[order(SamplesAll$SampleMeans),];
		mads = SamplesAllSorted$SampleMads; means = SamplesAllSorted$SampleMeans;
		modlow <- lowess(mads ~ means);	modlow_resid <- mads-modlow$y;
		SamplesAllResid <- data.frame(SamplesAllSorted,modlow_resid); names(SamplesAllResid)[ncol(SamplesAllResid)] = "LowessResidual";	
		SamplesAll2Return <- SamplesAllResid[abs(SamplesAllResid$LowessResidual)<(sd(modlow_resid)*sd_dist),];
	}
	## Type 2 -- Remove samples with calls that are N*sd away from mean Dels/Dups ratio
	if(type==2) {
		ddratio = round((SamplesAll$num_dels/SamplesAll$num_dups),digits=2);
		higherdd = mean(ddratio)+sd_dist*sd(ddratio); lowerdd = mean(ddratio)-sd_dist*sd(ddratio);
		ddratio_outliers = which(ddratio>=higherdd | ddratio<=lowerdd);
		SamplesAll2Return = SamplesAll[-c(ddratio_outliers),]
	}
	## Type 3 -- Remove samples with calls that are N*sd from Known % and return the rest
	if(type==3) {
		higherdd = mean(SamplesAll$known_pc)+sd_dist*sd(SamplesAll$known_pc); lowerdd = mean(SamplesAll$known_pc)-sd_dist*sd(SamplesAll$known_pc);
		kpc_outliers = which(SamplesAll$known_pc>=higherdd | SamplesAll$known_pc<=lowerdd);
		SamplesAll2Return = SamplesAll[-c(kpc_outliers),]
	}
	## Type 4 -- Remove samples with calls that are N*sd from #calls and return the rest
	if(type==4) {
		higherdd = mean(SamplesAll$num_calls)+sd_dist*sd(SamplesAll$num_calls); lowerdd = mean(SamplesAll$num_calls)-sd_dist*sd(SamplesAll$num_calls);
		nc_outliers = which(SamplesAll$num_calls>=higherdd | SamplesAll$num_calls<=lowerdd);
		SamplesAll2Return = SamplesAll[-c(nc_outliers),]
	}
	return(SamplesAll2Return);
}


