# TODO: This script generates log2 ratio for all samples
# Task outline:
# -- Median depth of exons across samples - between Males/Females separately for Chr X 
# Author: pv1
###############################################################################

`SampleMedian` <- function(RDfiles, sample_ids, sample_gender, regions_file, chrX=0) {
	
	sortx = read.table(file=regions_file, header=TRUE); 
	allX = data.frame(seq(1:length(sortx[,1])), sortx); names(allX) = c("rowcount", "chr","start","end","GC","dG","Tm");
	
	# Read read depth from each sample's files
	allXD = lapply(RDfiles, function(x) read.table(x, header=FALSE, check.names=FALSE, colClasses=c("character",rep("integer",3),"numeric"), nrows=nrow(allX), comment.char="", sep="\t")[,5])
	allX = data.frame(allX,do.call('cbind',allXD)); names(allX)[8:ncol(allX)] = sample_ids; 
	rm(allXD); gc();
	
	# Separate Autosomes and X,Y
	allX[,8:length(allX)] = allX[,8:length(allX)]+0.01 # Add 0.01 to avoid '0' values
	allXS = allX[allX[,2] != 'Y',]; allXS = allXS[allXS[,2] != 'X',]; # Autosomes only - remove X and Y ######
	allXSX = allX[allX[,2] == 'X',]; allXSY = allX[allX[,2] == 'Y',]; # X only and Y only
	
	rm(allX); rm(sortx); gc() # Cleanup
	
	allXSmedian = apply(allXS[,8:length(allXS)],1,median) # Median depth - Autosomes only
	allXSreturn = data.frame(allXS[,1:7],allXSmedian,allXSmedian) # Two columns are returned - For autosomes these are same. For X and Y, these will be different.
	names(allXSreturn) = c("rowcount", "chr","start","end","GC","dG","Tm","Median_Male","Median_Female");
			
	if(chrX == 0) {
		return(allXSreturn);
	} else {
		Males = sample_ids[which(sample_gender=='M')]; 	Females = sample_ids[which(sample_gender=='F')]; 
		allXSX_Males = allXSX[,names(allXSX) %in% Males]; allXSX_Females = allXSX[,names(allXSX) %in% Females];	
		allXSXmedian_Males = apply(allXSX_Males[,1:length(allXSX_Males)],1,median) # Median depth - Males Chr X only
		allXSXmedian_Females = apply(allXSX_Females[,1:length(allXSX_Females)],1,median) # Median depth - Females Chr X only
		for(s in 8:length(allXSX)) { # for each sample's Chr X - males and females separately
			if( any(Males == names(allXSX)[s]) ) { # If the sample is male
				allXSX[,s] = log2(allXSX[,s]/allXSXmedian_Males)
			} else {
				allXSX[,s] = log2(allXSX[,s]/allXSXmedian_Females)
			}
		}
		allXSXret1 = data.frame(allXSX[,1:7],allXSXmedian_Males,allXSXmedian_Females);
		names(allXSXret1) = c("rowcount", "chr","start","end","GC","dG","Tm","Median_Male","Median_Female");
		allXSreturn = rbind(allXSreturn,allXSXret1);
		return(allXSreturn);
	}	
}