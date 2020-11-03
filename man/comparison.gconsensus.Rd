\name{comparison.gconsensus}
\alias{comparison.gconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtains a set of generalized consensus estimates by different statistical methods for comparison purposes.
}
\description{
Obtains a set of generalized consensus estimates.
}
\usage{
comparison.gconsensus(x, methods, build.model = NULL, get.samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{consensus object, containing the definition of the participants' information.}
  \item{methods}{vector of strings, containing the statistical methods to be used in obtaining the consensus values.}
  \item{build.model}{function, NULL by default, the current function to build a JAGS measurand model.}
  \item{get.samples}{function, NULL by default, the current function to obtain CORA samples out of the JAGS model.}
}
\details{
The default function requires the default parameters (x, \dots), the next version will replace the formal parameters with the default signature.
}
\value{
  \item{fit}{a data.frame object, consensus value estimates by the specified methods}
  \item{gconsensus}{the current estimated gconsensus as provided in the x parameter}
  \item{total.included.participants}{a numeric, the current number of included participating sources}
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{gconsensus}}, \code{\link{print.comparison}}, \code{\link{plot.comparison}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ comparison.gconsensus }% use one of  RShowDoc("KEYWORDS")
