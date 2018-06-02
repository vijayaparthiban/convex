# TODO: This script calculates log2 ratio for all samples using a subset of highly correlated samples within the analysis batch
# Task outline:
# -- Median depth of exons across samples - between Males/Females separately for Chr X 
# Author: pv1
###############################################################################

`SampleLogRatioV2` <- function(RDfiles, sample_ids, sample_gender, regions_file, chrX=0, min_samples=25, RPKM=0, sample_means_file="") {
	
	if(length(sample_ids) < min_samples) {
		print(paste("There are less than",min_samples,"samples totally; Please use SampleLogRatio() instead!"));
		q("no", status=2, runLast=FALSE);
	}
	
	sortx = read.table(file=regions_file, header=TRUE); 
	allX = data.frame(seq(1:length(sortx[,1])), sortx); names(allX) = c("rowcount", "chr","start","end","GC","dG","Tm");
	
	if(RPKM!=0) {
		RPKM=1;
	}
		
	# Read read depth from each sample's files
	allXD = lapply(RDfiles, function(x) read.table(x, header=FALSE, check.names=FALSE, colClasses=c("character",rep("integer",3),rep("numeric",2)), nrows=nrow(allX), comment.char="", sep="\t")[,(5+RPKM)])
	allX = data.frame(allX,do.call('cbind',allXD)); names(allX)[8:ncol(allX)] = sample_ids; 
	rm(allXD); gc();
	
	# Separate Autosomes and X,Y
	allX[,8:length(allX)] = allX[,8:length(allX)]+0.01 # Add 0.01 to avoid '0' values
	allXS = allX[allX[,2] != 'Y',]; allXS = allXS[allXS[,2] != 'X',]; # Autosomes only - remove X and Y ######
	allXSX = allX[allX[,2] == 'X',]; allXSY = allX[allX[,2] == 'Y',]; # X only and Y only
	
	rm(allX); rm(sortx); gc() # Cleanup
	
	if(nchar(sample_means_file)>0) {
		sample_means = colMeans(allXS[,8:length(allXS)]);
		smmdf = data.frame(sample_ids,round(as.numeric(sample_means), digits=2)); names(smmdf) = c("sample_id","SampleMeans")
		write.table(smmdf,file=sample_means_file,row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)		
	}
	
	cor_matrix = cor(allXS[,8:ncol(allXS)]);
	allXS2 = allXS;
	
	for(i in 1:length(sample_ids)) {
		cor_sample = sort(cor_matrix[,i], decreasing=TRUE); highly_correlated = names(cor_sample)[2:24]; # top 20 excluding the first (same) sample 
		allXSmed = apply(allXS2[,names(allXS2) %in% highly_correlated],1,median);
		allXS[,i+7] = log2(allXS[,i+7]/allXSmed);
	}
	rm(allXS2); gc();
	
	if(chrX == 0) {
		return(allXS);
	} else {
		Males = sample_ids[which(sample_gender=='M')]; 	Females_Unknown = sample_ids[which(sample_gender!='M')]; 
		
		if(length(Males) < min_samples) {
			print(paste("There are less than",min_samples,"samples totally; Please use SampleLogRatio() instead!"));
			q("no", status=2, runLast=FALSE);
		}
		if(length(Females_Unknown) < min_samples) {
			print(paste("There are less than",min_samples,"samples totally; Please use SampleLogRatio() instead!"));
			q("no", status=2, runLast=FALSE);
		}
		allXSX_Males = allXSX[,names(allXSX) %in% Males]; 
		allXSX_Females = allXSX[,names(allXSX) %in% Females_Unknown];	
		cor_matrix_males = cor(allXSX_Males); cor_matrix_females = cor(allXSX_Females);
		
		#allXSXmedian_Males = apply(allXSX_Males[,1:length(allXSX_Males)],1,median) # Median depth - Males Chr X only
		#allXSXmedian_Females = apply(allXSX_Females[,1:length(allXSX_Females)],1,median) # Median depth - Females Chr X only
		
		for(i in 1:length(sample_ids)) { # for each sample's Chr X - males and females separately
			sid = names(allXSX)[i+7];
			if( any(Males == sid) ) { # If the sample is male
				si = which(colnames(cor_matrix_males)==sid);
				cor_sample = sort(cor_matrix_males[,si], decreasing=TRUE); highly_correlated = names(cor_sample)[2:24]; # top 20 excluding the first (same) sample 
				allXSXmedian_Males = apply(allXSX_Males[,names(allXSX_Males) %in% highly_correlated],1,median);
				allXSX[,i+7] = log2(allXSX[,i+7]/allXSXmedian_Males);
			} else {
				si = which(colnames(cor_matrix_females)==sid);
				cor_sample = sort(cor_matrix_females[,si], decreasing=TRUE); highly_correlated = names(cor_sample)[2:24]; # top 20 excluding the first (same) sample 
				allXSXmedian_Females = apply(allXSX_Females[,names(allXSX_Females) %in% highly_correlated],1,median);
				allXSX[,i+7] = log2(allXSX[,i+7]/allXSXmedian_Females);
			}
		}
		return(rbind(allXS,allXSX));
	}	
}