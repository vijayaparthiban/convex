# TODO: Get/Set CLASSPATH to run CoNVex's Java programs
# 
# Author: pv1
###############################################################################

getClassPath <- function(convex_path="") {
	if(convex_path=="") {
		convex_path = getLibPath();
		java_lib_path = paste(convex_path,"/CoNVex/java/lib/",sep="");
		jars = list.files(path=java_lib_path,pattern=".jar"); 
		if(convex_path=="NA") { # Note: getLibPath() may return 'NA' if CoNVex is not found
			print("Can't find CoNVex installation path; please specify explicitly for this function!");
		}
	} else {
		java_lib_path = paste(convex_path,"/CoNVex/java/lib/",sep="");
		jars = list.files(path=java_lib_path,pattern=".jar"); 
	}
	if(length(jars)>0) {
		jar_path = paste(java_lib_path,jars,collapse=":",sep="");
		jar_path_bash = paste("export CLASSPATH=$CLASSPATH:",jar_path,":",sep="");
		jar_path_csh = paste("setenv CLASSPATH $CLASSPATH':'",jar_path,":",sep="");
		
		print("If your default shell is bash, add this line to .bash_profile:"); print(jar_path_bash);
		print("***");
		print("If your default shell is csh, add this line to .cshrc:"); print(jar_path_csh);
		print("NOTE: If you don't have a defined CLASSPATH for other Java programs, remove $CLASSPATH from commands above")
		print("***")
		print("Which shell are you using? Use this command in Unix terminal/shell to find out:"); print("finger -m <unix_user_id>");	
	}
}