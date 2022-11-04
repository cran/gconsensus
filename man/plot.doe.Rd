\name{plot.doe}
\alias{plot.doe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a "degrees of equivalence" object
}
\description{
Plots a char showing the confidence interval for the unilateral degrees of equivalance associated to each participating laboratory.
}
\usage{
\method{plot}{doe}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A doe object}
  \item{\dots}{ Any additional graphical parameter}
}
\details{
Plot a unilateral degree of equivalence (doe) object.
  Formating options of the plot:
  \code{display.order} a string, this is used to order the data source, valid values are: "location" for sorting based on the repoted means, "dispersion" for sorting based on reported variances, any other value for sorting based on typing order. Default value is "location".
  \code{display.shownames} a logical, if true then the participants names are used to identify the ploted values, otherwise the  code value is used instead. Default value is FALSE.
  \code{display.orientation} a string, this is used to build the plot vertically or horizontally, valid values are: "horizontal" otherwise vertical is assumed, default value is "horizontal".
  \code{display.expandedUncertainty} a logical, if true then expanded uncertainty is used for plots, otherwise standard uncertainty is used. Default value is FALSE.
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
 \code{\link{doe.gconsensus}}, \code{\link{plot}}, \code{\link{print.doe}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ doe.gconsensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ plot }% __ONLY ONE__ keyword per line
\keyword{ print.doe }% __ONLY ONE__ keyword per line
