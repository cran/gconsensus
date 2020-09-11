comparison.gconsensus <- function(x, methods, build.model = NULL,
                                 get.samples = NULL) {
  estimates <- data.frame(method = rep(NA, length(methods)), 
    value = NA, U = NA, unit = NA, k = NA, p = NA, tau = NA)
  for (i in 1:length(methods)) {
    res <- gconsensus(x$ilab, method = methods[i],
                     build.model = build.model, get.samples = get.samples,
		         config = x$config)
    estimates$method[i] <- methods[i]
    estimates$value[i] <- res$fit$value
    estimates$U[i] <- res$fit$U
    estimates$unit[i] <- levels(res$fit$unit)[res$fit$unit]
    estimates$k[i] <- res$fit$k
    estimates$p[i] <- res$fit$p
    estimates$tau[i] <- res$fit$tau
  }

  res <- list(fit = estimates,
              gconsensus = x,
              total.included.participants = sum(x$ilab$data$included)
              )
  class(res) <- "comparison"
  return(res)
}
