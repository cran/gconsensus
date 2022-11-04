##
## PDF Maker
## R implementation
## 2018-04-30: include the corrected uncerainty including the between labs 
## variance component
## 2021-11-19: include the handle of standard uncertainty and expended uncertainty
## as a display option. CCQM guidelines requires to use standard uncertainty.
##

plot.gconsensus <- function(x, ...) {
  display.order <- getOption("display.order")
  display.shownames <- getOption("display.shownames")
  display.orientation <- getOption("display.orientation")
  display.length.out <- getOption("display.length.out")
  display.tab.size <- getOption("display.tab.size")
  display.signif.digits <- getOption("display.signif.digits")
  display.expandedUncertainty <- getOption("display.expandedUncertainty")
  
  if (is.null(display.order)) display.order <- "location"
  if (is.null(display.shownames)) display.shownames <- FALSE
  if (is.null(display.orientation)) display.orientation <- "horizontal"
  if (is.null(display.length.out)) display.length.out <- 101
  if (is.null(display.tab.size)) display.tab.size <- 12
  if (is.null(display.signif.digits)) display.signif.digits <- 2
  if (is.null(display.expandedUncertainty)) display.expandedUncertainty <- FALSE
  
  mss <- rep(TRUE, length(x$ilab$data$value))
  
  n <- length(x$ilab$data$value[mss])
  subset <- x$ilab$data$included[mss] == 1

  if (display.order == "location") {
    os <- order(x$ilab$data$value[mss])
  } else if (display.order == "dispersion") {
    os <- order(x$ilab$data$expandedUnc[mss])
  } else {
    os <- order(x$ilab$data$participant[mss])
  }
  if (display.shownames) {
    xlab <- x$ilab$data$participant[mss]
  } else {
    xlab <- x$ilab$data$code[mss]
  }
  
  my.pch <- x$ilab$data$symbol[mss]
  my.color <- x$ilab$data$symbol.fillcolor[mss]

  p <- length(x$ilab$data$value[mss])
  ss <- x$ilab$data$included[mss] == 1
	mm <- sum(ss)
	coverageFactor <- 2
	if (x$config$expansion.factor.type == "naive") {
	  coverageFactor <- 2
	} else {
	  if (x$config$expansion.factor.type == "large sample") {
	    coverageFactor <- qnorm(1 - x$config$alpha / 2)
	  } else {
	    coverageFactor <- qt(1 - x$config$alpha / 2, mm - 1)
	  }
	}

	tau2 <- 0
	if (x$method == "DL1" | x$method == "DL2" | x$method == "PM" | 
		x$method == "MPM") {
	  tau2 <- x$sb2
	} else if (x$method == "VRMLE" | x$method == "MCM.median" |
		x$method == "HB" | x$method == "MPM.LP") {
		  tau2 <- x$var.b
	  #    print(paste0("tau2=", tau2))
	} else {
	  tau2 <- 0
	}
	
	sxlab <- c(xlab)
	for (i in 1:p) if (!subset[i]) {
	  sxlab[i] <- paste0("<", xlab[i], ">")
	} else sxlab[i] <- paste(xlab[i])

	wlim <- c(0, p+3)
	if (display.expandedUncertainty) {
	  zlim <- range(c(x$ilab$data$value - sqrt(x$ilab$data$expandedUnc^2 + (x$ilab$data$coverageFactor * x$fit$tau)^2), 
	                  x$fit$value - x$fit$expandedUnc),
	                c(x$ilab$data$value + sqrt(x$ilab$data$expandedUnc^2 + (x$ilab$data$coverageFactor * x$fit$tau)^2), 
	                  x$fit$value + x$fit$expandedUnc),
	                na.rm = TRUE)
	} else {
	  zlim <- range(c(x$ilab$data$value - sqrt((x$ilab$data$expandedUnc/x$ilab$data$coverageFactor)^2 + (x$fit$tau)^2), 
	                  x$fit$value - x$fit$expandedUnc/x$fit$coverageFactor),
	                c(x$ilab$data$value + sqrt((x$ilab$data$expandedUnc/x$ilab$data$coverageFactor)^2 + (x$fit$tau)^2), 
	                  x$fit$value + x$fit$expandedUnc/x$fit$coverageFactor),
	                na.rm = TRUE)
	  
	}

  wlab <- ""
#  if (display.shownames) wlab <- "source name" else wlab <- "source code"

  zlab <- x$ilab$info$value[x$ilab$info$variable == "Units"]
  
	if (display.orientation == "horizontal") {
	  xx <- c(1:p)
	  yy <- x$ilab$data$value[mss][os]

	  xlim <- wlim
	  ylim <- zlim
	  xlab <- wlab
	  ylab <- zlab
	  xaxis <- 1
	  yaxis <- 2
	  x.separator <- c(p+1, p+1)
	  y.separator <- c(-1, 1)*10*max(abs(x$ilab$data$value), na.rm = TRUE)
	} else {
	  xx <- x$ilab$data$value[mss][os]
	  yy <- c(1:p)
	  xlim <- zlim
	  ylim <- wlim
	  xlab <- zlab
	  ylab <- wlab
	  xaxis <- 2
	  yaxis <- 1
	  x.separator <- c(-1, 1)*10*max(abs(x$ilab$data$value), na.rm = TRUE)
	  y.separator <- c(p+1, p+1)
	}
 
  plot(xx, yy,
		xlim = xlim,
		ylim = ylim,
		axes = FALSE,
		xlab = xlab,
		ylab = ylab,
		main = paste0(x$exercise, " - ", x$measurand),
    cex = 1.5,
    pch = 20, #my.pch,
    ...
	)

  
	axis(yaxis)
	axis(xaxis, at = 1:p, labels = sxlab[os], las = 2)
  for (i in 1:p) if (!subset[os][i]) {
    axis(xaxis, at = i, labels = sxlab[os][i], las = 2, col.axis = 2)
  }
  box()
  lines(x.separator, y.separator)

  xx.pdf <- seq(min(zlim), max(zlim), length.out = display.length.out)
  yy.mat <- matrix(NA, display.length.out, p)
  for (i in (c(1:p)[ss])) 
    yy.mat[,i] <- dnorm(xx.pdf, x$ilab$data$value[i], x$ilab$data$expandedUnc[i])
  yy.pdf <- apply(yy.mat, 1, mean, na.rm = TRUE)
  mu.pdf <- sum(xx.pdf * yy.pdf) / sum(yy.pdf)
  u.pdf <- sqrt(sum((xx.pdf - mu.pdf) ^ 2 * yy.pdf)/sum(yy.pdf))
  yy.eq <- dnorm(xx.pdf, x$fit$value, x$fit$expandedUnc*sqrt(p)/x$fit$coverageFactor)
  
  tcol<- rgb(118, 238, 0, maxColorValue = 255, alpha=127)

  if (display.orientation == "horizontal") {

    if (display.expandedUncertainty) {
      polygon(c(1, p, p, 1, 1), x$fit$value+c(-1, -1, 1, 1, -1)*x$fit$expandedUnc, 
              border = NA, col = tcol)
    } else {
      polygon(c(1, p, p, 1, 1), x$fit$value+c(-1, -1, 1, 1, -1)*x$fit$expandedUnc/x$fit$coverageFactor, 
              border = NA, col = tcol)
    }
    lines(c(0, p+1), rep(x$fit$value, 2), col = "darkgreen")
    
    lines(p + 1 + 2*yy.pdf/max(c(yy.eq,yy.pdf)), xx.pdf, lwd = 2)
    lines(p + 1 + 2*yy.eq/max(c(yy.eq,yy.pdf)), xx.pdf, lwd = 2, col = 4)

    if (display.expandedUncertainty) {
      for (i in 1:p) {
        rect(i-0.25, 
             x$ilab$data$value[mss][os][i] - 
                x$ilab$data$expandedUnc[mss][os][i], 
             i+0.25,
             x$ilab$data$value[mss][os][i] + 
                x$ilab$data$expandedUnc[mss][os][i], 
            col = "pink", border = NA)
        lines(rep(i, 2),
              x$ilab$data$value[mss][os][i] +
                c(-1,1) * 
                sqrt(x$ilab$data$expandedUnc[mss][os][i] ^ 2 + (x$ilab$data$coverageFactor[mss][os][i] * x$fit$tau)^2),
              lwd = 2, col = "red")
      }
    } else {
      for (i in 1:p) {
        rect(i-0.25, 
             x$ilab$data$value[mss][os][i] - 
               x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i], 
             i+0.25,
             x$ilab$data$value[mss][os][i] + 
               x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i], 
             col = "pink", border = NA)
        lines(rep(i, 2),
              x$ilab$data$value[mss][os][i] +
                c(-1,1) * sqrt((x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i])^ 2 + x$fit$tau^2),
              lwd = 2, col = "red")
      }
    }
  } else { ## display.orientation == "vertical"
    if (display.expandedUncertainty) {
      polygon(x$fit$value+c(-1, -1, 1, 1, -1)*x$fit$expandedUnc, c(1, p, p, 1, 1),
              border = NA, col = tcol)
    } else {
      polygon(x$fit$value+c(-1, -1, 1, 1, -1)*x$fit$expandedUnc/x$fit$coverageFactor, c(1, p, p, 1, 1), 
              border = NA, col = tcol)
    }
    lines(rep(x$fit$value, 2), c(0, p+1), col = "darkgreen")
    
    lines(xx.pdf, p + 1 + 2*yy.pdf/max(c(yy.eq,yy.pdf)), lwd = 2)
    lines(xx.pdf, p + 1 + 2*yy.eq/max(c(yy.eq,yy.pdf)), lwd = 2, col = 4)
    
    if (display.expandedUncertainty) {
      for (i in 1:p) {
        rect(x$ilab$data$value[mss][os][i] - 
               x$ilab$data$expandedUnc[mss][os][i], i-0.25,
             x$ilab$data$value[mss][os][i] + 
               x$ilab$data$expandedUnc[mss][os][i], i+0.25,
             col = "pink", border = NA)
        lines(x$ilab$data$value[mss][os][i] +
                c(-1,1) * 
                sqrt(x$ilab$data$expandedUnc[mss][os][i] ^ 2 + x$ilab$data$coverageFactor[mss][os][i]^2 * x$fit$tau^2), rep(i, 2),
              lwd = 2, col = "red")
      }
    } else {
      for (i in 1:p) {
        rect(x$ilab$data$value[mss][os][i] - 
               x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i], i-0.25,
             x$ilab$data$value[mss][os][i] + 
               x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i], i+0.25,
             col = "pink", border = NA)
        lines(x$ilab$data$value[mss][os][i] +
                c(-1,1) * sqrt((x$ilab$data$expandedUnc[mss][os][i]/x$ilab$data$coverageFactor[mss][os][i]) ^ 2 + x$fit$tau^2), rep(i, 2),
              lwd = 2, col = "red")
      }
    }

if (FALSE) {
    for (i in 1:p) {
      lines(x$ilab$data$value[os][i] + c(-1,1) * x$ilab$data$coverageFactor[os][i] * 
			x$ilab$data$expandedUnc[os][i], rep(i, 2),
            lwd = 3, col = 4)
      lines(x$ilab$data$mean[os][i] + c(-1,1) * x$ilab$data$k[os][i] * 
			sqrt(x$ilab$data$sd[os][i] ^ 2 + x$fit$tau^2), rep(i, 2),
            lwd = 1, col = 4)
    }
    
    tcol<- rgb(118, 238, 0, maxColorValue = 255, alpha=127)
    
    polygon(x$fit$value+c(-1, -1, 1, 1, -1)*x$fit$U, c(1, p, p, 1, 1), 
            border = NA, col = tcol)
    lines(rep(x$fit$value, 2), c(1, p), col = "darkgreen")

#    lines(rep(x$fit$value - x$fit$U, 2), c(1, p), col = "lightgreen")
#    lines(rep(x$fit$value + x$fit$U, 2), c(1, p), col = "lightgreen")
    
    lines(xx.pdf, p + 1 + 2*yy.pdf/max(c(yy.eq,yy.pdf)), lwd = 2)
    lines(xx.pdf, p + 1 + 2*yy.eq/max(c(yy.eq,yy.pdf)), lwd = 2, col = 4)
}
  }

  points(xx, yy, pch = 20, cex = 1, col = 2 ) #my.color)
}
