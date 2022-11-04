comparison.gconsensus <- function(x, methods, build.model = NULL,
                                 get.samples = NULL) {
  estimates <- data.frame(method = rep(NA, length(methods)), 
    value = NA, expandedUnc = NA, unit = NA, coverageFactor = NA, coverageProbability = NA, tau = NA)
  for (i in 1:length(methods)) {
    res <- gconsensus(x$ilab, method = methods[i],
                     build.model = x$build.model, get.samples = x$get.samples,
		         config = x$config)
    estimates$method[i] <- methods[i]
    estimates$value[i] <- res$fit$value
    estimates$expandedUnc[i] <- res$fit$expandedUnc
    estimates$unit[i] <- res$fit$unit # levels(res$fit$unit)[res$fit$unit]
    estimates$coverageFactor[i] <- res$fit$coverageFactor
    estimates$coverageProbability[i] <- res$fit$coverageProbability
    estimates$tau[i] <- res$fit$tau
  }

  res <- list(fit = estimates,
              gconsensus = x,
              total.included.participants = sum(x$ilab$data$included)
              )
  class(res) <- "comparison"
  return(res)
}
