\name{doe.gconsensus}
\alias{doe.gconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtains the unilateral degrees of equivalence out of a generalized consensus object.
}
\description{
Obtains the unilateral degrees of equivalence from a generalized consensus object.
}
\usage{
doe.gconsensus(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A generalized consensus object.}
}
\details{
The unilateral degrees of equivalence are obtained from the participants data and the estimated generalized consensus value described in the gconsensus value.
}
\value{
A doe object containing the unilateral degrees of equivalence.
  \item{fit}{a data.frame with the source code, source lab name, value, expanded uncertainty, unit, expansion factor, probability coverage and tau.}
  \item{gconsensus}{the generalized consensus object the evaluation is based on.}

}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{gconsensus}}, \code{\link{print.doe}}, \code{\link{plot.doe}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gconsensus }% use one of  RShowDoc("KEYWORDS")

