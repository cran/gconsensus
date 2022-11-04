## 2017-03-14: expansion type was included as a parameter
## 2020-01-21: hierarchical bayesian method was included
## 2020-11-02: hierarchical bayesian comments from Rpackages reviewers about allowing model_file to be user defined
## 2022-10-07: DSI terminology is implemented
## 2022-11-02: validation of parameters is implemented
gconsensus <- function(ilab, method = "mean",
		      build.model = NULL, get.samples = NULL,
                      config = list(alpha = 0.05,
                                    expansion.factor.type = "naive",
						tau = mad(ilab$data$value),
                                    unreliable.uncertainties = FALSE,
                                    MC_samples = 1e5,
						MC_burn_in = 1000,
                                    MC_seed = NA,
                                    MC_use.HKSJ.adjustment = FALSE,
						filename = "hb_consensus_model.txt"
                                   )
                      ) {
  
  mss <- rep(TRUE, length(ilab$data$value))
  
  n <- length(ilab$data$value[mss])
  subset <- ilab$data$included[mss] == 1
  
	mm.res <- .internal.mm.pdf(ilab$data[mss, ][subset, ], 
	  B = config$MC_samples, new.seed = config$MC_seed)
	
	x <- ilab$data$value[mss][subset]
	u2 <- (ilab$data$expandedUnc[mss][subset]/ilab$data$coverageFactor[mss][subset])^2
	n <- ilab$data$n[mss][subset]
  
## 2017-03-14: this code was included to handle the selected type of expansion
  mm <- sum(subset)
  coverageFactor <- 2
  if (config$expansion.factor.type == "naive") {
    coverageFactor <- 2
  } else {
    if (config$expansion.factor.type == "large sample") {
      coverageFactor <- qnorm(1 - config$alpha/2)
    } else {
      coverageFactor <- qt(1 - config$alpha/2, mm - 1)
    }
  }
## end
  
#  qq <- 1
  tau <- config$tau
	if (all(is.null(n))) { n <- rep(2, length(x)) }
	if (method == "mean") {
	  res <- .internal.mom(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	} else if (method == "grand.mean") {
	  res <- .internal.gmm(x, sqrt(u2), n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	} else if (method == "h15") {
	  res <- .internal.hh(x, sqrt(u2), n, alpha = config$alpha)
	  tau <- 0
	} else if (method == "median") {
	  res <- .internal.my.median(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
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
	  tau <- 0
	} else if (method == "GD2") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "sinha", 
	    alpha = config$alpha)
	  tau <- 0
	} else if (method == "GD3") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "zhang1", 
	    alpha = config$alpha)
	  tau <- 0
	} else if (method == "GD4") {
	  res <- .internal.graybill.deal(x, u2 * n, n, method = "zhang2", 
	    alpha = config$alpha)
	  tau <- 0
	} else if (method == "DL1") {
	  res <- .internal.dersimonian.laird(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
		if (tau > 0) {
	  if (config$MC_use.HKSJ.adjustment) {
	    if (!is.null(res$w.i)) {
	      mm <- sum(subset)
		  wsse<- sum(res$w.i*(ilab$data$value[mss][subset] - res$mu)^2)
	      qq <- sqrt(1/(mm - 1)*wsse)
	    } else {
	      qq <- 1
	    }
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
		if (tau > 0) {
	  if (config$MC_use.HKSJ.adjustment) {
	    if (!is.null(res$w.i)) {
	      mm <- sum(subset)
		  wsse<- sum(res$w.i*(ilab$data$value[mss][subset] - res$mu)^2)
	      qq <- sqrt(1/(mm - 1)*wsse)
	    } else {
	      qq <- 1
	    }
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
	} else if (method == "HB") {
	  res <- .internal.hb(x, u2, n, build.model, get.samples, 
			      alpha = config$alpha, tau = tau,
				seed = config$MC_seed,
			      MC_samples = config$MC_samples,
				MC_burn_in = config$MC_burn_in,
				file = file.path(tempdir(), config$filename))
	  tau <- sqrt(res$var.b)
	} else if (method == "VRMLE") {
	  res <- vr.mle(x, u2 * n, n, alpha = config$alpha, tol = 1e-12, 
	    max.iter = 3000)
	  tau <- sqrt(res$var.b)
	} else if (method == "BOB") {
	  res <- .internal.bob(x, u2 * n, n, alpha = config$alpha, 
	    expansion = config$expansion.factor.type)
	  tau <- sqrt(res$sb2)
	} else if (method == "SE") {
	  res <- .internal.Schiller.Eberhardt(x, u2 * n, n, alpha = config$alpha)
	  tau <- sqrt(res$sb2)
	}
  
  res2 <- list(method = method)
  res2$subset <- subset
  res2$ilab <- ilab
  res2$ilab$data <- ilab$data[mss, ]
  res2$measurand <- ilab$info$value[ilab$info$variable == "Measurand"]
  res2$exercise <- ilab$info$value[ilab$info$variable == "Exercise"]
  res2$config <- config
  # 2021-11-19: we need to keep the used build.model and get.samples methods 
  # to be passed to the report document method
  res2$build.model <- build.model
  res2$get.samples <- get.samples
  res.unit <- ilab$info$value[ilab$info$variable == "Units"]

  res2$fit <- data.frame(value = res$mu, expandedUnc = coverageFactor * res$u.mu, 
		unit = res.unit, coverageFactor = coverageFactor, coverageProbability = 1-config$alpha, tau = tau)

  class(res2) <- "gconsensus"
	return(res2)
}
