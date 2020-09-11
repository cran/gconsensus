doe.gconsensus <- function(x) {
  mss <- rep(TRUE, length(x$ilab$data$mean)) 
  
  subset <- (x$ilab$data$included[mss] == 1)
  mm <- sum(subset)
  if (x$config$expansion.factor.type == "naive") {
    k <- 2
  } else if (x$config$expansion.factor.type == "large sample") {
    k <- qnorm(1 - x$config$alpha/2)
  } else {
    k <- qt(1 - x$config$alpha/2, mm - 1)
  }
  doe <- (x$ilab$data$mean[mss] - x$fit$value)
  
  # if reliable estimates
  # for results not included in calculation of RV
  U.doe <- sqrt(x$ilab$data$k[mss]^2*x$ilab$data$sd[mss]^2 + x$fit$U^2)
  
  # for results included in calculation of RV
  U.doe[subset] <- sqrt(1 - 2/mm)*
	sqrt((x$ilab$data$k*x$ilab$data$sd)[mss][subset]^2 + x$fit$U^2)
  
  # correction for results included in calculation of RV
  if (x$config$unreliable.uncertainties) {
    # if unreliable uncertainty estimates then
    U.doe[subset] <- sqrt(1 + 2*(mm - 2)/pi)*x$fit$U
  }
  
  labs <- x$ilab$data$participant[mss]
  codes <- x$ilab$data$code[mss]

  u.doe <- sqrt(x$ilab$data$sd[mss]^2 + (x$fit$U/x$fit$k)^2)
  u.doe[subset] <- sqrt(1 - 2/mm)*u.doe[subset]
  if (x$config$unreliable.uncertainties) {
    # if unreliable uncertainty estimates then
    u.doe[subset] <- sqrt(1 + 2*(mm - 2)/pi)*x$fit$U/x$fit$k
  }
  k <- U.doe/u.doe

  res <- list(fit = data.frame(code = codes, lab = labs, value = doe, U = U.doe, 
				unit = rep("1", length(doe)), k = k, 
				p = rep(1 - x$config$alpha, length(doe)), En = doe/U.doe),
				gconsensus = x)
  class(res) <- "doe"
  
  return(res)
}
