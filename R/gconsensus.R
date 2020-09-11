## 2017-03-14: expansion type was included as a parameter
## 2020-01-21: hierarchical bayesian method was included
gconsensus <- function(ilab, method = "mean",
		      build.model = NULL, get.samples = NULL,
                      config = list(alpha = 0.05,
                                    expansion.factor.type = "naive",
                                    unreliable.uncertainties = FALSE,
                                    MC_samples = 1e5,
                                    MC_seed = NA,
                                    MC_use.HKSJ.adjustment = FALSE
                                    )
                      ) {
  mss <- rep(TRUE, length(ilab$data$mean))
  
  n <- length(ilab$data$mean[mss])
  subset <- ilab$data$included[mss] == 1
  
	mm.res <- .internal.mm.pdf(ilab$data[mss, ][subset, ], 
	  B = config$MC_samples, new.seed = config$MC_seed)
	
	x <- ilab$data$mean[mss][subset]
	u2 <- (ilab$data$sd[mss][subset])^2
	n <- ilab$data$n[mss][subset]
  
## 2017-03-14: this code was included to handle the selected type of expansion
  mm <- sum(subset)
  k <- 2
  if (config$expansion.factor.type == "naive") {
    k <- 2
  } else {
    if (config$expansion.factor.type == "large sample") {
      k <- qnorm(1 - config$alpha/2)
    } else {
      k <- qt(1 - config$alpha/2, mm - 1)
    }
  }
## end
 
	qq <- 1
  tau <- 0
	if (all(is.null(n))) { n <- rep(2, length(x)) }
	if (method == "mean") {
	  res <- .internal.mom(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	} else if (method == "grand.mean") {
	  res <- .internal.gmm(x, sqrt(u2), n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	} else if (method == "h15") {
	  res <- .internal.hh(x, sqrt(u2), n, alpha = config$alpha)
	} else if (method == "median") {
	  res <- .internal.my.median(x, u2 * n, n, alpha = config$alpha)
	} else if (method == "MCM.mean") {
	  res <- .internal.mm.mean(mm.res, length(x), alpha = config$alpha)
	  tau <- sqrt(res$var.b)
	} else if (method == "MCM.LP") {
	  res <- .internal.mm.linearpool(mm.res, length(x), alpha = config$alpha)
	  tau <- sqrt(res$var.b)
	} else if (method == "MCM.median") {
	  res <- .internal.mm.median(mm.res, length(x), alpha = config$alpha)
	  tau <- sqrt(res$var.b)
	} else if (method == "GD1") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "naive", 
	    alpha = config$alpha)
	} else if (method == "GD2") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "sinha", 
	    alpha = config$alpha)
	} else if (method == "GD3") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "zhang1", 
	    alpha = config$alpha)
	} else if (method == "GD4") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "zhang2", 
	    alpha = config$alpha)
	} else if (method == "DL1") {
	  res <- .internal.dersimonian.laird(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	  if (config$MC_use.HKSJ.adjustment) {
	    if (!is.null(res$w.i)) {
	      mm <- sum(subset)
		  wsse<- sum(res$w.i*(ilab$data$mean[mss][subset] - res$mu)^2)
	      qq <- sqrt(1/(mm - 1)*wsse)
	    } else {
	      qq <- 1
	    }
	  } else {
	    qq <- 1
	  }
	  res$u.mu <- qq*res$u.mu
	} else if (method == "DL2") {
	  res <- .internal.ss.dersimonian.laird(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	  if (config$MC_use.HKSJ.adjustment) {
	    if (!is.null(res$w.i)) {
	      mm <- sum(subset)
		  wsse<- sum(res$w.i*(ilab$data$mean[mss][subset] - res$mu)^2)
	      qq <- sqrt(1/(mm - 1)*wsse)
	    } else {
	      qq <- 1
	    }
	  } else {
	    qq <- 1
	  }
	  res$u.mu <- qq*res$u.mu
	} else if (method == "PM") {
	  res <- .internal.pmm(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	} else if (method == "MPM") {
	  res <- .internal.mpmm(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
#	} else if (method == "HB") {
#	  res <- .internal.hb(x, u2, n, build.model, get.samples, 
#			      alpha = config$alpha, seed = config$MC_seed)
#	  tau <- sqrt(res$var.b)
	} else if (method == "VRMLE") {
	  res <- vr.mle(x, u2 * n, n, alpha = config$alpha, tol = 1e-12, 
	    max.iter = 3000)
	  tau <- sqrt(res$var.b)
	} else if (method == "BOB") {
	  res <- .internal.bob(x, u2 * n, n, alpha = config$alpha, 
	    expansion = config$expansion.factor.type)
	} else if (method == "SE") {
	  res <- .internal.Schiller.Eberhardt(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	}
  
  res2 <- list(method = method)
  res2$subset <- subset
  res2$ilab <- ilab
  res2$ilab$data <- ilab$data[mss, ]
  res2$measurand <- ilab$info$value[ilab$info$variable == "Measurand"]
  res2$study <- ilab$info$value[ilab$info$variable == "Study"]
  res2$config <- config
  res.unit <- ilab$info$value[ilab$info$variable == "Units"]

  res2$fit <- data.frame(value = res$mu, U = qq*k*res$u.mu, 
		unit = res.unit, k = k, p = 1-config$alpha, tau = tau)

  class(res2) <- "gconsensus"
	return(res2)
}
