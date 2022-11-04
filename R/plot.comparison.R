## 2021-11-19: se incluyo la opcion de mostrar la incertidumbre expandida o la incertidumbre estandar
plot.comparison <- function(x, ...) {
  # x is a data.frame with method:string, mean:numeric, u:numeric
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
  
  n <- length(x$fit$value)
  ci <- matrix(cbind(x$fit$value, x$fit$expandedUnc/x$fit$coverageFactor, x$fit$expandedUnc), n, 3, 
		byrow = FALSE)

  if (display.order == "location") {
    ss <- order(x$fit$value)
  } else if (display.order == "dispersion") {
    ss <- order(x$fit$expandedUnc)
  } else {
    ss <- order(x$fit$method) 
  }

  zz <- c(1:n)
  zlim <- c(1, n)
  ww <- x$fit$value[ss]
  wlim <- range(c(x$fit$value - x$fit$expandedUnc, x$fit$value + x$fit$expandedUnc))
  ssxlab <- "method"
  wlab <- paste0(x$gconsensus$exercise, ": ", x$gconsensus$measurand, " /(", x$fit$unit[1], ")")

  if (display.orientation == "horizontal") {
    xx <- zz
    yy <- ww
    xlim <- zlim
    ylim <- wlim
    xaxis <- 1
    yaxis <- 2
    xlab <- ssxlab
    ylab <- wlab
  } else {
    xx <- ww
    yy <- zz
    xlim <- wlim
    ylim <- zlim
    xaxis <- 2
    yaxis <- 1
    xlab <- wlab
    ylab <- ssxlab
  }

    plot(xx[ss], 
         yy[ss],
         xlim = xlim,
         ylim = ylim,
         axes = FALSE,
         xlab = xlab,
         main = expression("Comparing Consensus Values " * (1 - alpha) * "% CI"),
	   ylab = ylab
    )

  if (display.orientation == "horizontal") {
    if (display.expandedUncertainty) {
      for (i in 1:n) {
        lines(rep(i, 2), x$fit$value[ss][i] + c(-1, 1) * x$fit$expandedUnc[ss][i], lwd=1)
      }
    } else {
      for (i in 1:n) {
        lines(rep(i, 2), x$fit$value[ss][i] + c(-1, 1) * x$fit$expandedUnc[ss][i]/x$fit$coverageFactor[ss][i], lwd=1)
      }
    }
  } else {
    if (display.expandedUncertainty) {
      for (i in 1:n) {
        lines(x$fit$value[ss][i] + c(-1, 1) * x$fit$expandedUnc[ss][i], rep(i, 2), lwd=1)
      }
    } else {
      for (i in 1:n) {
        lines(x$fit$value[ss][i] + c(-1, 1) * x$fit$expandedUnc[ss][i]/x$fit$coverageFactor[ss][i], rep(i, 2), lwd=1)
      }
    }
  }

    box()
    axis(yaxis)
    axis(xaxis, at = c(1:n), labels = x$fit$method[ss])
}
