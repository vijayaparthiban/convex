# TODO: Add comment
# 
# Author: pv1
###############################################################################

`getLibPath` <- function() {

	convex_path="NA"
	lp = .libPaths(); 
	for(i in 1:length(lp)) {
		lf = list.files(path=lp[i],pattern="CoNVex");
		if(length(lf) > 0) {
			convex_path=lp[i];
		}
	}
	return(convex_path);	
}
