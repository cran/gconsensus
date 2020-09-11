plot.comparison <- function(x, ...) {
  # x is a data.frame with method:string, mean:numeric, u:numeric
  display.order <- getOption("display.order")
  display.shownames <- getOption("display.shownames")
  display.orientation <- getOption("display.orientation")
  display.length.out <- getOption("display.length.out")
  display.tab.size <- getOption("display.tab.size")
  display.signif.digits <- getOption("display.signif.digits")
  
  if (is.null(display.order)) display.order <- "location"
  if (is.null(display.shownames)) display.shownames <- FALSE
  if (is.null(display.orientation)) display.orientation <- "horizontal"
  if (is.null(display.length.out)) display.length.out <- 101
  if (is.null(display.tab.size)) display.tab.size <- 12
  if (is.null(display.signif.digits)) display.signif.digits <- 

  n <- length(x$fit$value)
  ci <- matrix(cbind(x$fit$value, x$fit$U/x$fit$k, x$fit$U), n, 3, 
		byrow = FALSE)

  if (display.order == "location") {
    ss <- order(x$fit$value)
  } else if (display.order == "dispersion") {
    ss <- order(x$fit$U)
  } else {
    ss <- order(x$fit$method) 
  }

  zz <- c(1:n)
  zlim <- c(1, n)
  ww <- x$fit$value[ss]
  wlim <- range(c(x$fit$value - x$fit$U, x$fit$value + x$fit$U))
  ssxlab <- "method"
  wlab <- paste0(x$gconsensus$study, ": ", x$gconsensus$measurand, " /(", x$fit$unit[1], ")")

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
    for (i in 1:n) {
      lines(rep(i, 2), x$fit$value[ss][i] + c(-1, 1) * x$fit$U[ss][i], lwd=1)
      lines(rep(i, 2), x$fit$value[ss][i] + c(-1, 1) * x$fit$U[ss][i]/x$fit$k[ss][i], lwd=3)
    }
  } else {
    for (i in 1:n) {
      lines(x$fit$value[ss][i] + c(-1, 1) * x$fit$U[ss][i], rep(i, 2), lwd=1)
      lines(x$fit$value[ss][i] + c(-1, 1) * x$fit$U[ss][i]/x$fit$k[ss][i], rep(i, 2), lwd=3)
    }
  }

    box()
    axis(yaxis)
    axis(xaxis, at = c(1:n), labels = x$fit$method[ss])
}
