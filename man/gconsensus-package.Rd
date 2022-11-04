\name{gconsensus-package}
\alias{gconsensus-package}
\docType{package}
\title{
\packageTitle{gconsensus}
}
\description{
\packageDescription{gconsensus}
Partially based on the results published in Control clinical trials (1985). <https://doi.org/10.1016/0197-2456(86)90046-2>
For a complete list of references see the reference section.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{gconsensus}
\packageIndices{gconsensus}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
\packageAuthor{gconsensus}

Maintainer: \packageMaintainer{gconsensus}
}
\references{
1. Graybill and Deal (1959), \emph{Combining Unbiased Estimators}, Biometrics, 15, pp. 543-550.

2. DerSimonian and Laird (1986), \emph{Meta-analysis in Clinical Trials}, Controlled Clinical Trials, 7, pp. 177-188.

3. R. A. Horn, S. A. Horn and D. B. Duncan (1975), \emph{Estimating heteroscedastic variance in linear models}. Journ. Amer. Statist. Assoc. 70, 380

4. M. S. Levenson, D. L. Banks, K. R. Eberhardt, L. M. Gill, W. F. Guthrie, H. K. Liu, M. G. Vangel, J. H. Yen, and N. F. Zhang (2000), \emph{An ISO GUM Approach to Combining Results from Multiple Methods}, Journal of Research of the National Institute of Standards and Technology, Volume 105, Number 4.

5. John Mandel and Robert Paule (1970), \emph{Interlaboratory Evaluation of a Material with Unequal Number of Replicates}, Analytical Chemistry, 42, pp. 1194-1197.

6. Robert Paule and John Mandel (1982), \emph{Consensus Values and Weighting Factors}, Journal of Research of the National Bureau of Standards, 87, pp. 377-385.

7. Andrew Rukhin (2009), \emph{Weighted Means Statistics in Interlaboratory Studies}, Metrologia, Vol. 46, pp. 323-331.

8. Andrew Ruhkin (2003), \emph{Two Procedures of Meta-analysis in Clinical Trials and Interlaboratory Studies}, Tatra Mountains Mathematical Publications, 26, pp. 155-168.

9. Andrew Ruhkin and Mark Vangel (1998), \emph{Estimation of a Common Mean and Weighted Means Statistics}, Journal of the American Statistical Association, Vol. 93, No. 441.

10. Andrew Ruhkin, B. Biggerstaff, and Mark Vangel (2000), \emph{Restricted Maximum Likelihood Estimation of a Common Mean and Mandel-Paule Algorithm}, Journal of Statistical Planning and Inference, 83, pp. 319-330.

11. Mark Vangel and Andrew Ruhkin (1999), \emph{Maximum Likelihood Analysis for Heteroscedastic One-Way Random Effects ANOVA in Interlaboratory Studies}, Biometrics 55, 129-136.

12. Susannah Schiller and Keith Eberhardt (1991), \emph{Combining Data from Independent Analysis Methods}, Spectrochimica, ACTA 46 (12).

13. Bimal Kumar Sinha (1985), \emph{Unbiased Estimation of the Variance of the Graybill-Deal Estimator of the Common Mean of Several Normal Populations}, The Canadian Journal of Statistics, Vol. 13, No. 3, pp. 243-247.

14. Nien-Fan Zhang (2006), \emph{The Uncertainty Associated with The Weighted Mean of Measurement Data}, Metrologia, 43, PP. 195-204.

15. CCQM (2013), \emph{CCQM Guidance note: Estimation of a Consensus KCRV and associated Degrees of Equivalence}, Version 10.

16. Knapp G. and Hartung J. (2003), \emph{Improved tests for a random effects meta-regression with a single covariate}, Statistics in Medicine, Vol 22, Issue 17, pp 2693-2710  
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
\keyword{ package }

\seealso{
 \code{\link{gconsensus}}
}
\examples{
ilab <- list(
  data = data.frame(
    participant = c("BAM", "IRMM", "LGC", "NARL", "NIST", "NMIJ", "NRC"),
    code = paste0("p", c(1:7)),
    method = rep(1, 7), 
    value = c(10.21, 10.9, 10.94, 10.58, 10.81, 9.62, 10.8), 
    unit = rep("\\\\micro \\\\kilogram/\\\\kilogram", 7),
    expandedUnc = c(0.381, 0.250, 0.130, 0.410, 0.445, 0.196, 0.093), 
    n = rep(1, 7), 
    coverageFactor = rep(2, 7),
    coverageProbability = rep(0.95, 7), 
	method.type = rep("p", 7),
    included = rep(1, 7)),
	symbol = rep(1, 7),
	symbol.Fillcolor = rep(1, 7),
	comments = rep("", 7),
	info = list(variable = c("Software", "Institute", "Exercise", "Author", "Date", 
      "Measurand", "Units"), 
    value = c("consensus", "BIPM", "CCQM-K25", "Michele Schantz & Stephen Wise",
      "12-12-2003", "PCB 105", "\\\\micro \\\\kilogram/\\\\kilogram"))
)

aConsensus <- gconsensus(ilab,
  method = "DL1",
  config = list(
               signif.digits = 2,
			   display.order = "location",
			   display.orientation = "horizontal",
			   display.shownames = TRUE,
			   tab.size = 12,
			   alpha = 0.05,
               expansion.factor.type = "small.sample",
			   tau = mad(ilab$data$value),
               unreliable.uncertainties = FALSE,
               MC_samples = 1e5,
               MC_seed = 12345,
               MC_use.HKSJ.adjustment = TRUE,
			   MC_interval.type = "symmetric",
			   MC_use.parallel.computing = FALSE,
		       filename = "hb_consensus_model.txt"
          )
)

print(aConsensus)

}
