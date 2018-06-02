\name{getClassPath}
\alias{getClassPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Helps you setup the $CLASSPATH environment variable }
\description{
Automatically searches for jar files in the installation path and shows $CLASSPATH environment variable settings.
}
\usage{
getClassPath(convex_path="")
}
\arguments{
  \item{convex_path}{ CoNVex R package installation path (default: convex_path = "" - this tries to automatically retrieve it from your default R lib path) .  }
}
\details{
 getClassPath() helps you run Java programs bundled with CoNVex. Please use convex_path = "/path/to/ConvexRpackage" explicitly, if
 getLibPath() does not find the CoNVex installation path! 
}
\references{ R Library path }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{ getClassPath() internally uses getLibPath() which iteratively searches through all folders in .libPaths() for a functional CoNVex installation, assuming you installed CoNVex in one of these locations.
}
\seealso{ 
.libPaths(), getLibPath()
}
\examples{
> getClassPath()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ getLibPath, sw_lin64, sw_osx }
