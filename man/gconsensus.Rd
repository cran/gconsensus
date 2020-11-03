\name{gconsensus}
\alias{gconsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Consensus Contructor
}
\description{
Compute a set of estimates for the generalized consensus value from a set of sources.
}
\usage{
gconsensus(ilab, method = "mean", 
  build.model = NULL, get.samples = NULL,
  config = list(alpha = 0.05, 
    expansion.factor.type = "naive", 
    unreliable.uncertainties = FALSE, 
    MC_samples = 1e5, 
    MC_seed = NA, 
    MC_use.HKSJ.adjustment = FALSE,
    filename = "hb_consensus_model.txt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ilab}{ilab object, containing the definition of the participating source/experts/laboratories.}
  \item{method}{string, the current statistical method for estimating the consensus value. The list of current valid methods is:
 \code{mean} [15], \code{grand.mean} [15], \code{median} [15], \code{h15} for Huber estimator [15], \code{MCM.mean} for Monte Carlo mean, \code{MCM.median} for Monte Carlo median, \code{GD1} for the original Graybill-Deal estimator [1], \code{GD2} for the modified Graybill-Deal estimator as described by Sinha [13], \code{GD3} for the modified Graybill-Deal estimator as described by Zhang [14], \code{GD4} for the modified Graybill-Deal estimator as described by Zhang [14], \code{DL1} for original DerSimonian-Laid estimator [2], \code{DL2} for the modified DerSimonian-Laird estimator as described by Horn-Horn-Duncan [3], \code{PM} for original Paul-Mandel estimator [5], \code{MPM} for the modified Paul-Mandel estimator [6], \code{VRMLE} for Vangel-Rukhin MLE [7, 8, 9, 10, 11], \code{BOB} for the type B on Bias estimator [4], \code{SE} for the Schieller-Eberhardt [12], MCM.LP for linear pool, HB for Hierarchical Bayesian}

  \item{build.model}{a function (NULL by default) to build a JAGS model}
  \item{get.samples}{a function (NULL by dafault) to obtain the CORA samples based on the JAGS model}

  \item{config}{list, contains the configuration options used in the computational process. The complete list of options is:
		\code{alpha} for the significance level, default value = 0.05,
		\code{expansion.factor.type} for the type of coverage factor to be used: "naive" which use a coverage factor=2, "large sample" which uses a normal quantile, "small sample" which uses a t quantile,
		\code{unreliable.uncertainties} an indicator for unreliables uncertainties, if TRUE then adjustment is used [15],
		\code{MC_samples} the number of cycles used in the simple Monte Carlo simulation,
		\code{MC_seed} the seed used to initializing the random number generator in the Monte Carlo simulation,
		\code{MC_use.HKSJ.adjustment} a flag indicating if small sample adjustments are required [16],
		\code{filename} a string containing the name of the file for communication with JAGS package.
%		\code{MC_use.parallel.computing} a flag indicating if parallel computing should be used in the Monte Carlo simulation, when available,
%		\code{MC_use.adaptative.method} a flag indicating if simple or adaptative algorithm is used in the Monte Carlo simulation,
%		\code{MC_sample.size.increment} numeric, the size of the sample to be incremented in search of the convergence of the adaptative Monte Carlo simulation,
%		\code{MC_sample.size.upper.bound} largest sample size allowed to commit in the adaptative Monte Carlo simulation.
}
}
\details{
The consensus estimation will vary with the specified method.
Display can be customized by using the options mechanism. Available options are:
{\code{display.signif.digits} an integer representing the number of significant digits to be obtained in the adaptative approach of the Monte Carlo simulation,
		\code{display.order} a string ("name", "code", "location", "dispersion") representing the order to be used for displaying the uncertainty contribution,
		\code{display.shownames} a boolean, if true source names are used otherwise source codes are used,
		\code{display.tab.size} an integer representing the column width used in formating data tables,
		\code{display.orientation} a string ("horizontal", "vertical"), plot orientation.
}
}
\value{
  \item{fit}{a data.frame with the value, expanded uncertainty, unit, expansion factor, probability coverage and tau, the between effect uncertainty.}
  \item{method}{statistical method used to estimate the consensus value.}
  \item{subset}{the subset of included values in the consensus value estimation.}
  \item{ilab}{the current participants' information.}
  \item{config}{the current configuration options provided for the consensus value estimation process.}
  \item{study}{the study the consensus belongs to.}
  \item{measurand}{the measurand the consensus belongs to.}
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{doe.gconsensus}}, \code{\link{print.gconsensus}}, \code{\link{plot.gconsensus}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gconsensus }% use one of  RShowDoc("KEYWORDS")

