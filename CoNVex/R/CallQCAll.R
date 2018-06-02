# TODO: Apply call QC filters to CNV calls
# 
# Author: pv1
###############################################################################

# Input CNV calls data frame with the following columns (with 'column names' in brackets):
# - CoNVex score (convex_score)
# - Mean log2 ratio (mean_l2r)
# - MAD of log2 ratio
# - Internal frequency proportion at 50% reciprocal overlap (rc50_internal_freq)
# - common_forward and common_backward overlap proportion (output of recipeB package OR CoNVex's ChrOverlap() function)

`CallQCAll` <- function(cnv_calls = data.frame()) {
	
	cnv_calls_names = c("convex_score","mean_l2r","mad_l2r","internal_freq","rc50_internal_freq","common_forward","common_backward");
	present = cnv_calls_names %in% names(cnv_calls);
	rownames(cnv_calls) = c(1:nrow(cnv_calls))
	
	if( (nrow(cnv_calls)>0 && length(unique(present)) == 1) && (unique(present)) ) {
		# All required columns are present
		# RULES:
		# ALL CNVs - COMMON + RARE #	
		cnv_size = cnv_calls$end - cnv_calls$start
		cnv_calls_qc = data.frame(cnv_calls,cnv_size,"PASS"); names(cnv_calls_qc) = c(names(cnv_calls),"cnv_size","call_qc");	
		
		dx1 = which(cnv_calls_qc$cnv_type=='DUP' & cnv_calls_qc$convex_score<7)  # [you can set this in Step 6: swt_dup=7 in the last step]
		dx2 = which(cnv_calls_qc$cnv_size>500000 & cnv_calls_qc$convex_score<10)
		dx3 = which(cnv_calls_qc$cnv_size>200000 & cnv_calls_qc$convex_score<10 & cnv_calls$num_probes >= 10)
		dx4 = which(cnv_calls_qc$mean_l2r >= 1.5)
		dx5 = which(cnv_calls_qc$common_forward == 0 & cnv_calls_qc$rc50_internal_freq > 0.05)
		dx5a = which(cnv_calls_qc$mean_l2r > 0 & cnv_calls_qc$cnv_type == 'DEL')
		dx5b = which(cnv_calls_qc$mean_l2r < 0 & cnv_calls_qc$cnv_type == 'DUP')
		cnv_calls_qc$call_qc = as.character(cnv_calls_qc$call_qc); cnv_calls_qc$call_qc[unique(c(dx1,dx2,dx3,dx4,dx5,dx5a,dx5b))] = "FAIL";
		
		# RARE CNVs (common_forward<0.5) #
		cnv_calls_rare = cnv_calls_qc[cnv_calls_qc$common_forward < 0.5,]
		dx6 = which(cnv_calls_rare$cnv_type=='DUP' & cnv_calls_rare$num_probes == 1 & cnv_calls_rare$convex_score < 10)
		dx7 = which(cnv_calls_rare$convex_score < 10 & cnv_calls_rare$num_probes == 1 & cnv_calls_rare$cnv_size < 500)
		dx8 = which(cnv_calls_rare$mean_l2r > -0.5 & cnv_calls_rare$cnv_type == 'DEL' & cnv_calls_rare$convex_score < 10)
		dx9 = which(cnv_calls_rare$mean_l2r < 0.29 & cnv_calls_rare$cnv_type == 'DUP' & cnv_calls_rare$convex_score < 10)
		ccr = rownames(cnv_calls_rare[unique(c(dx6,dx7,dx7,dx8,dx9)),])
		cnv_calls_qc$call_qc = as.character(cnv_calls_qc$call_qc); cnv_calls_qc$call_qc[as.numeric(ccr)] = "FAIL";
			
		# ALL - RARE (common_forward<0.5) - this is done separately below because this may go under different FAIL category in future 
		dx10 = which(cnv_calls_rare$cnv_size<15000 & cnv_calls_rare$convex_score<10 & cnv_calls_rare$internal_freq > 0.01)
		dx11 = which(cnv_calls_rare$cnv_type=='DUP' & cnv_calls_rare$num_probes > 5 & cnv_calls_rare$convex_score < 8)
		ccr2 = rownames(cnv_calls_rare[unique(c(dx10,dx11)),])
		cnv_calls_qc$call_qc = as.character(cnv_calls_qc$call_qc); cnv_calls_qc$call_qc[as.numeric(ccr2)] = "FAIL";
		
		# Return the final cnv_calls_qc
		return(cnv_calls_qc);		
	} else {
		write("CNV calls data frame with the following name(s) required (column names must match):", stderr());
		write(cnv_calls_names[!present], stderr());	
	}	
}