\name{getLibPath}
\alias{getLibPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get CoNVex installation Path }
\description{
Search and return CoNVex installation path (R library path) using this function. 
}
\usage{
getLibPath()
}
\details{
getLibPath() iteratively searches through all folders in .libPaths() for a functional CoNVex installation, assuming you installed CoNVex in one of these locations.
 If it fails, it returns 'NA'. Please use convex_path = "/path/to/ConvexRpackage" explicitly for getClassPath() or other functions.
}
\references{ R Library path }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This path can be used as CoNVex's installation folder for locating bundled external data and binaries. 
If folders and files are used from a custom location, this function is not required. 
}
\seealso{ 
.libPaths(), getClassPath()
}
\examples{
> convex_path = getLibPath()
> print(convex_path)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ getLibPath, sw_lin64, sw_osx }
