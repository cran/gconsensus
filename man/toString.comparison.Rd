\name{toString.comparison}
\alias{toString.comparison}
\title{
  Builds a detailed description string of the comparison object.
}
\description{
  This function converts a gconsensus comparison object into a string. The returned value contains a verbatim description of the object. This function is used to display the result of comparing several gconsensus values.
}
\usage{
\method{toString}{comparison}(x, ...)
}
\arguments{
  \item{x}{comparison, object to be described.}
  \item{\dots}{any additional parameters.}
}
\value{
  string, it contains a verbatim description of the comparison object.
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{comparison.gconsensus}}, \code{\link{print.comparison}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ toString }
\keyword{ comparison } % __ONLY ONE__ keyword per line
