plot.doe <- function(x, ...) {
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
  if (is.null(display.signif.digits)) display.signif.digits <- 2

  mss <- rep(TRUE, length(x$gconsensus$ilab$data$value))

  n <- length(x$fit$value)
  subset <- x$gconsensus$ilab$data$included[mss] == 1

  if (display.order == "location") {
    ss <- order(x$gconsensus$ilab$data$value[mss])
  } else if (display.order == "dispersion") {
    ss <- order(x$gconsensus$ilab$data$expandedUnc[mss])
  } else {
    if (display.shownames) {
      ss <- order(x$fit$lab[mss]) 
    } else {
      ss <- order(x$fit$code[mss]) 
    }
  }
  
  if (display.shownames) {
    xlab <- x$fit$lab[mss]
    ssxlab<- "Source name" 
  } else {
    xlab <- x$fit$code[mss]
    ssxlab<- "Source code"
  }
  ssxlab <- ""
  ssylab <- x$fit$unit[1]

  sxlab <- c(xlab)
  # mark excluded values writeing name or code between brackets
  for (i in 1:n) if (!x$gconsensus$subset[i]) {
    sxlab[i] <- paste0("<", xlab[i], ">")
  } else sxlab[i] <- paste(xlab[i])

  zz <- c(1:n)
  zlim <- c(1, n)
  ww <- x$fit$value[ss]
  wlim <- range(c(x$fit$value - x$fit$expandedUnc, x$fit$value + x$fit$expandedUnc))

  if (display.orientation == "horizontal") {
    xx <- zz
    yy <- ww
    xlim <- zlim
    ylim <- wlim
    xaxis <- 1
    yaxis <- 2
    xlab <- ssxlab
    ylab <- ssylab
  } else {
    xx <- ww
    yy <- zz
    xlim <- wlim
    ylim <- zlim
    xaxis <- 2
    yaxis <- 1
    xlab <- ssylab
    ylab <- ssxlab
  }

  # this shown all the participants either their values are used for 
  # the consensus of not
  plot(xx, yy,  
       xlab = xlab,
       axes = FALSE,
       ylab = ylab,
       pch = 19,
       main = paste0("Unilateral Degrees of Equivalence\n", 
			x$gconsensus$exercise, " - ", x$gconsensus$measurand),
       xlim = xlim, ylim = ylim)
  
  axis(yaxis)
  axis(xaxis, at = c(1:n), labels = sxlab[ss], las = 2)
  for (i in 1:n) if (!subset[ss][i]) {
    axis(xaxis, at = i, labels = sxlab[ss][i], las = 2, col.axis = 2)
  }

  if (display.orientation == "horizontal") {
    for (ii in 1:n) {
      lines(rep(ii, 2), x$fit$value[ss][ii] + c(-1, 1)*x$fit$expandedUnc[ss][ii])
    }
    for (i in 1:n) 
      if (!x$gconsensus$subset[ss][i])
        axis(xaxis, at = i, labels = sxlab[ss][i], las = 2, col.axis = 2)
    abline(h = 0)
  } else {
    for (ii in 1:n) {
      lines(x$fit$value[ss][ii] + c(-1, 1)*x$fit$expandedUnc[ss][ii], rep(ii, 2))
    }
    for (i in 1:n) 
      if (!x$gconsensus$subset[ss][i])
        axis(xaxis, at = i, labels = sxlab[ss][i], las = 2, col.axis = 2)
    abline(v = 0)
  }
}
