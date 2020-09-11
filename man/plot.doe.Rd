\name{plot.doe}
\alias{plot.doe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a "degrees of equivalence" object
}
\description{
Plots a chart showing the confidence interval for the unilateral degrees of equivalance associated to each participating laboratory.
}
\usage{
\method{plot}{doe}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A doe object.}
  \item{\dots}{Any additional graphical parameter.}
}
\details{
Plot a doe object using the options \code{display.order} and \code{display.shownames}.
}
\value{
A plot of the doe object using the options displor.order and display.shownames, showing the confidence interval for the unilateral degree of equivalence attached to each participant.
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{doe.gconsensus}}, \code{\link{plot}}, \code{\link{print.doe}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ doe.gconsensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ gconsensus }% __ONLY ONE__ keyword per line
\keyword{ plot }% __ONLY ONE__ keyword per line
