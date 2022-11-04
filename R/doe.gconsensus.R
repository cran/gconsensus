doe.gconsensus <- function(x) {
  mss <- rep(TRUE, length(x$ilab$data$value)) 
  
  subset <- (x$ilab$data$included[mss] == 1)
  mm <- sum(subset)
  if (x$config$expansion.factor.type == "naive") {
    coverageFactor <- 2
  } else if (x$config$expansion.factor.type == "large sample") {
    coverageFactor <- qnorm(1 - x$config$alpha/2)
  } else {
    coverageFactor <- qt(1 - x$config$alpha/2, mm - 1)
  }
  doe <- (x$ilab$data$value[mss] - x$fit$value)
  
  # if reliable estimates
  # for results not included in calculation of RV
  expandedUnc.doe <- sqrt(x$ilab$data$expandedUnc[mss]^2 + x$fit$expandedUnc^2)
  
  # for results included in calculation of RV
  expandedUnc.doe[subset] <- sqrt(1 - 2/mm)*
	sqrt(x$ilab$data$expandedUnc[mss][subset]^2 + x$fit$expandedUnc^2)
  
  # correction for results included in calculation of RV
  if (x$config$unreliable.uncertainties) {
    # if unreliable uncertainty estimates then
    expandedUnc.doe[subset] <- sqrt(1 + 2*(mm - 2)/pi)*x$fit$expandedUnc
  }
  
  labs <- x$ilab$data$participant[mss]
  codes <- x$ilab$data$code[mss]

  standardUnc.doe <- sqrt(x$ilab$data$expandedUnc[mss]^2 + x$fit$expandedUnc^2)/x$fit$coverageFactor
  standardUnc.doe[subset] <- sqrt(1 - 2/mm)*standardUnc.doe[subset]
  if (x$config$unreliable.uncertainties) {
    # if unreliable uncertainty estimates then
    standardUnc.doe[subset] <- sqrt(1 + 2*(mm - 2)/pi)*x$fit$expandedUnc/x$fit$coverageFactor
  }
  coverageFactor <- expandedUnc.doe/standardUnc.doe

  res <- list(fit = data.frame(code = codes, lab = labs, value = doe, 
				expandedUnc = expandedUnc.doe, 
				unit = rep(x$fit$unit, length(doe)), 
				coverageFactor = coverageFactor, 
				coverageProbability = rep(1 - x$config$alpha, length(doe)), 
				En = doe/expandedUnc.doe),
				gconsensus = x)
  class(res) <- "doe"
  
  return(res)
}
 