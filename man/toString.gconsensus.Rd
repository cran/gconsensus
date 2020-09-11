\name{toString.gconsensus}
\alias{toString.gconsensus}
\title{
  Builds a detailed description string of the gconsenus object.
}
\description{
  This function converts a gconsensus object into a string. The returned value contains a verbatim description of the object. This function is used to display the result of estimating the gconsensus value.
}
\usage{
\method{toString}{gconsensus}(x, ...)
}
\arguments{
  \item{x}{gconsensus, object to be described.}
  \item{\dots}{any additional parameters.}
}

\value{
  string, it contains a verbatim description of the consensus object.
}

\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}

\seealso{
  \code{\link{gconsensus}}, \code{\link{print.gconsensus}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ toString }
\keyword{ gconsensus } % __ONLY ONE__ keyword per line
