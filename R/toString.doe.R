toString.doe <- function(x, ...) {
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

  mss <- rep(TRUE, length(x$gconsensus$ilab$data$mean))
  
  if (display.order == "location") {
    ss <- order(x$gconsensus$ilab$data$mean[mss])
  } else if (display.order == "dispersion") {
    ss <- order(x$gconsensus$ilab$data$sd[mss])
  } else {
    if (display.shownames) {
      ss <- order(x$gconsensus$ilab$data$participant[mss])
    } else {
      ss <- order(x$gconsensus$ilab$data$code[mss])
    }
  }
  
  str <- paste0("Consensus method: ", x$gconsensus$method,".\n")
  str <- paste0(str, "Study: ", x$gconsensus$study,".\n")
  str <- paste0(str, "Measurand: ", x$gconsensus$measurand,".\n")
  str <- paste0(str, "Evaluation based on Degrees of Equivalence.\n\n")

  str2 <- capture.output(x$fit[ss, ], file = NULL)
  for (i in 1:length(str2)) str <- paste0(str, (str2[i]), "\n")

  return(str)
}
