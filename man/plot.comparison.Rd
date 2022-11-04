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
  Formating options of the plot:
  \code{display.order} a string, this is used to order the data source, valid values are: "location" for sorting based on the repoted means, "dispersion" for sorting based on reported variances, any other value for sorting based on typing order. Default value is "location".
  \code{display.shownames} a logical, if true then the participants names are used to identify the ploted values, otherwise the  code value is used instead. Default value is FALSE.
  \code{display.orientation} a string, this is used to build the plot vertically or horizontally, valid values are: "horizontal" otherwise vertical is assumed, default value is "horizontal".
  \code{display.expandedUncertainty} a logical, if true then expanded uncertainty is used for plots, otherwise standard uncertainty is used. Default value is
}
\value{
  a plot object
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
