\name{plot.gconsensus}
\alias{plot.gconsensus}
\title{
  Plot a gconsensus object
}
\description{
Shows the probability distribution of the estimated gconsensus value.
}
\usage{
\method{plot}{gconsensus}(x, ...)
}
\arguments{
  \item{x}{gconsensus object, containing the estimated gconsensus value by the specified statistical method.}
  \item{\dots}{any additional graphical parameter}
}
\details{
  The plot shows the confidence interval for each of the participants' values, the consensus value and its confidence interval, on the right side the probability distribution associated with the consensus value is showed and a normal distribution with the same expected value and variace.
  formating options of the plot:
  \code{display.order} a string, this is used to order the data source, valid values are: "location" for sorting based on the repoted means, "dispersion" for sorting based on reported variances, any other value for sorting based on typing order. Default value is "location".
  \code{display.shownames} a logical, if true then the participants names are used to identify the ploted values, otherwise the  code value is used instead. Default value is FALSE.
  \code{display.orientation} a string, this is used to build the plot vertically or horizontally, valid values are: "horizontal" otherwise vertical is assumed, default value is "horizontal".
  \code{display.length.out} a numeric, this represents the number of points to build the mixture distribution, default value = 101.
}
\value{
  a plot object.
}
\author{
\packageAuthor{gconsensus}
Maintainer: \packageMaintainer{gconsensus}
}
\seealso{
  \code{\link{gconsensus}}, \code{\link{print.gconsensus}}, \code{\link{plot}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gconsensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ plot }% __ONLY ONE__ keyword per line
