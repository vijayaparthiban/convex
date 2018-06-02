# TODO: Correction - Multivariate splines
# 
# Author: pv1
###############################################################################


`GAMCorrection` <- function(L2Rfile,features_file,RDfile,BPfile,doNorm=0,gc_only=0) {
	sampl2r = read.delim(file=L2Rfile, header=FALSE);
	allXSfeat = read.delim(file=features_file, header=TRUE);
	Reads = read.delim(file=RDfile, header=FALSE);
	Breakpoints = read.delim(file=BPfile, header=TRUE);
	
	# No. of Reads per Target and Cut breaks - break points vector
	ReadsallXS = Reads[rownames(Reads) %in% allXSfeat[,1],] ## MedianReadsInBins, MedianMADInBins, ADMweights, BreakPoints
	nnmax = c(Breakpoints[,2],max(ReadsallXS[,4])+1)
	
	# Cut the Number of Reads into bins
	xc = cut(ReadsallXS[,4]+0.01,breaks=nnmax)
	xxc_cut = data.frame(seq(1:length(ReadsallXS[,1])),ReadsallXS,xc)
	names(xxc_cut)[length(xxc_cut)] = "CutRanges"
	xcu = sort(unique(xc))
		
	# Correction
	require(mgcv);
	if(gc_only==0) {
		model = gam(sampl2r$V1 ~ s(allXSfeat[,5]) + s(allXSfeat[,6],allXSfeat[,7])) # linear correction
	} else {
		model = gam(sampl2r$V1 ~ s(allXSfeat[,5])) # linear correction
	}
	CR_MADs = data.frame(xxc_cut[,1:6],model$residual,xxc_cut$CutRanges)
	names(CR_MADs) = c("RowCount","Chr","Start","End","Reads","Depths","ModelResid","CutRanges")
	
	xclmad = rep(0.0001,length(xcu))
	for(i in 1:length(xcu)) {
		tmpc = CR_MADs[CR_MADs[,8]==xcu[i],]	
		xclmad[i] = mad(tmpc[,7],constant=1)	
	}	
	# Cut Ranges and their ADM weights in one DF
	BPranges = data.frame(xclmad,xcu)
	names(BPranges) = c("SampleMADsInBins","CutRanges")
	
	# Merge and sort it in Chr ascending order
	xtm = merge(x=CR_MADs,y=BPranges,by="CutRanges")
	xtm = xtm[with(xtm,order(xtm[,2])),]
	
	ADMRep = EWScore(xtm[,c(3:5,8:9)], doNorm=0);	
	return(ADMRep);
}