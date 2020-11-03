\name{vr.mle}
\alias{vr.mle}
\title{
  Finds the maximum likelihood estimate solution described by Vangel-Rukhin for the one way random effects model.
}
\description{
  This function finds the mle solution to the one way random effects model.
}
\usage{
vr.mle(xi, si2, ni, labi=c(1:length(xi)), 
   max.iter=1000, tol=.Machine$double.eps^0.5, 
   init.mu=mean(xi), init.sigma2=var(xi), 
   trace=FALSE, alpha=0.05)
}
\arguments{
  \item{xi}{numeric vector, represents the mean values.}
  \item{si2}{numeric vector, represents the variances associated with a single measurement.}
  \item{ni}{integer vector, represents the number of observations associated with the reported mean values.}
  \item{labi}{vector, containing the associated labels of the participanting laboratories, source of the reported values (mean, variances, number of observations)}
  \item{max.iter}{integer, maximum number of iterations allowed.}
  \item{tol}{numeric, relative tolerance.}
  \item{init.mu}{numeric, initial consensus value.}
  \item{init.sigma2}{numeric, initial between variance.}
  \item{trace}{logic, indicates if traceable information must be shown during the execution.}
  \item{alpha}{numeric, significance level.}
}

\value{
  \item{mu}{estimated consensus value by the method of maximum likelihood}
  \item{u.mu}{standard uncertainty estimation attached to the consensus value}
  \item{kp}{estimated expansion factor for the specified configuration options}
}
\author{
\packageAuthor{gconsensus}
Maintainer: \packageMaintainer{gconsensus}
}
\seealso{
  See also \code{\link{gconsensus}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gconsensus } % __ONLY ONE__ keyword per line
