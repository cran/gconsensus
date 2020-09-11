\name{plot.comparison}
\alias{plot.comparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot a gconsensus comparison object
}
\description{
Shows the confidence intervals associated with the different gconsensus values specified in the comparison.
}
\usage{
\method{plot}{comparison}(x, ...)
}
\arguments{
  \item{x}{gconsensus comparison object, containing the estimates for a set of gconsensus value.}
  \item{\dots}{any additional graphical parameter.}
}
\details{
  The plot contains the confidence intervals associated with the set of gconsensus values under comparison.
}
\value{
  a plot object.
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}
\seealso{
  \code{\link{comparison.gconsensus}}, \code{\link{print.comparison}}, \code{\link{plot}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ comparison.gconsensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ gconsensus }% __ONLY ONE__ keyword per line
\keyword{ plot }% __ONLY ONE__ keyword per line
