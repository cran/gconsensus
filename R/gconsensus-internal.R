##
## formats the vector quantities up to the specified significant digits
##

#=================================
#
# Mean of Means Method
#
#=================================
.internal.mom <- function(o.mui, o.si2, o.ni, alpha=0.05) {
	mu.e <- mean(o.mui)
	sd.e <- sqrt((var(o.mui) + mean(o.si2/o.ni))/length(o.mui))
	ci.e <- mu.e + qt(c(alpha/2, 1 - alpha/2), length(o.mui) - 1)*sd.e
	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, sb2 = var(o.mui)) )
}

#=================================
#
# Median Method
#
#=================================
.internal.my.median <- function(o.mui, o.si2, o.ni, alpha=0.05) {
	mu.e <- median(o.mui)
	sd.e <- sqrt(pi/2)*sqrt((mad(o.mui)^2 + median(o.si2/o.ni))/length(o.mui))
	ci.e <- mu.e + qt(c(alpha/2, 1 - alpha/2), length(o.mui) - 1)*sd.e
	return(list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, 
	  sb2 = (pi/2*mad(o.mui))^2 ))
}

#=================================
#
# Grand Mean Method
#
#=================================
.internal.gmm <- function(o.mui, o.si2, o.ni, alpha = 0.05) {
	nT.e <- sum(o.ni)
	mu.e <- sum(o.ni*o.mui)/nT.e
	sd.e <- sqrt(sum(o.ni*(o.mui - mu.e)^2)/(nT.e - length(o.mui))/nT.e)
	ci.e <- mu.e + qt(c(alpha/2, 1 - alpha/2), nT.e - length(o.mui))*sd.e
	return(list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, dof = nT.e - 1, sb2 = 0) )
}


#=================================
#
# Weighted Mean Method
#
#=================================
.internal.wmm <- function(o.mui, o.si2, o.ni, alpha=0.05) {
	wi <- o.ni/o.si2
	wi <- wi/sum(wi)
	mu.e <- sum(wi*o.mui)
	sd.e <- sqrt(1/sum(o.ni/o.si2)*(1 + 4*sum(wi*(1 - wi)/(o.ni - 1))))
	ci.e <- mu.e + qnorm(c(alpha/2, 1 - alpha/2))*sd.e
	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e) )
}

##
## Huber estimators
##

.internal.hh <- function(x, s, n, alpha = 0.05, p = 15, max.iter = 500) {
  cnfg <- matrix(c(
    1.0, 0.516, 
    1.1, 0.578,
    1.2, 0.635,
    1.3, 0.688,
    1.4, 0.736,
    1.5, 0.778,
    1.6, 0.816,
    1.7, 0.849,
    1.8, 0.877,
    1.9, 0.900,
    2.0, 0.921
  ), 11, 2, byrow = TRUE)
  
  cc <- cnfg[p - 9, 1]
	bb <- cnfg[p - 9, 2]
	pp <- length(x)
	k <- 1
	xs <- x
	mu <- median(xs)
	ss <- mad(x)/0.6745
	while (k < max.iter) {
		delta <- cc*ss
		xs[xs < mu - delta] <- mu - delta
		xs[xs > mu + delta] <- mu + delta
		mu <- mean(xs)
		ss <- sqrt(sum((xs - mu)^2)/(pp - 1))/bb^2
		k <- k + 1
	}
	kp <- qnorm(alpha/2, 1 - alpha/2)
	return( list(mu = mu, u.mu = ss/sqrt(length(x)),
	             ci.mu = mu + kp*ss/sqrt(length(x))) )
}




#=================================
#
# Paule-Mandel Method
#
# pmm(d[,1], d[,2]^2, d[,3])
#
#=================================
.internal.pmm <- function(o.mu, o.s2, o.n, trace = FALSE, max.iter = 50000,
                max.tol = .Machine$double.eps^0.5, alpha = 0.05) {
	DL.res <- .internal.dersimonian.laird(o.mu, o.s2, o.n, alpha = alpha)
	sb2 <- DL.res$sb2
	curr.tol <- 1
	if (trace) print(c("initial value: ", sb2))
	if (sb2 > 0) delta.sb2 <- sb2/2 else delta.sb2 <- 1
	k <- 1
	while ((curr.tol > max.tol) && (sb2 > 0) && (k < max.iter)) {
		if (trace) print(c("cycle :", k))
		w <- 1/(o.s2/o.n + sb2)
		y.tilde <- sum(w * o.mu) / sum(w)
		F.sb2 <- sum(w*(o.mu - y.tilde)^2) - (length(o.mu) - 1)
		Fp.sb2 <- 2*sum(w*(o.mu - y.tilde)) *
		  (sum(w^2*o.mu)/sum(w) - sum(w*o.mu)*sum(w^2)/sum(w)^2) -
		  sum(w^2*(o.mu - y.tilde)^2)
		delta.sb2 <- -F.sb2/Fp.sb2
		if (trace) print(c("delta.sb2 :", delta.sb2))
		new.sb2 <- sb2 + delta.sb2
		new.w <- 1/(o.s2/o.n + new.sb2)
		new.y.tilde <- sum(w*o.mu)/sum(w)
		if (trace) print(paste("mu=",new.y.tilde,"sb2=", new.sb2))
		k <- k + 1
		curr.tol <- max(abs(y.tilde - new.y.tilde), abs(sb2 - new.sb2),
		                abs(delta.sb2/(abs(sb2) + abs(new.sb2))))
		sb2 <- new.sb2
		y.tilde <- new.y.tilde
	}
	if (trace) print(c("sb2:", sb2))
	if (sb2 < 0) sb2 <- 0
	#
	# the estimate depends on the initial values
	# starting at the mean usually diverges
	# starting at the meadian usually converges
	#
	w <- 1/(o.s2/o.n + sb2)
	y.tilde <- sum(w*o.mu)/sum(w)
	if (trace) print(c("w:", w))
	mu.e <- sum(w*o.mu)/sum(w)
	sd.e <- sqrt(sum(w^2*(o.mu - y.tilde)^2))/sum(w)
	ci.e <- mu.e + qnorm(c(alpha/2, 1 - alpha/2))*sd.e
	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, sb2 = sb2, iter = k,
	             sug = length(o.mu) > 5, wi = w) )
}


#=================================
#
# modified Paul-Mandel Method
#
#=================================
.internal.mpmm <- function(o.mu, o.s2, o.n, trace = FALSE, alpha = 0.05) {
	DL.res <- .internal.dersimonian.laird(o.mu, o.s2, o.n, alpha = alpha)
	sb2 <- DL.res$sb2
	if (trace) print(c("initial value: ", sb2))
	delta.sb2 <- 10
	if (sb2 > 0) delta.sb2 <- sb2/2 else delta.sb2 <- 1
	k <- 1
	while ((abs(delta.sb2/sb2) > .Machine$double.eps*2) && (sb2 > 0)) {
		if (trace) print(c("cycle :", k))
		w <- 1/(o.s2/o.n + sb2)
		y.tilde <- sum(w*o.mu)/sum(w)
		F.sb2 <- sum(w*(o.mu - y.tilde)^2) - (length(o.mu))
		Fp.sb2 <- 2*sum(w*(o.mu - y.tilde)) *
		  (sum(w^2*o.mu)/sum(w) - sum(w*o.mu)*sum(w^2)/sum(w)^2) -
		  sum(w^2*(o.mu - y.tilde)^2)
		delta.sb2 <- -F.sb2/Fp.sb2
		if (trace) print(c("delta.sb2 :", delta.sb2))
		sb2 <- sb2 + delta.sb2
		if (trace) print(c("sb2 :", sb2))
		k <- k + 1
	}
	if (trace) print(c("sb2:", sb2))
	if (sb2 < 0) sb2 <- 0

	#
	# the estimate depends on the initial values
	# starting at the mean usually diverges
	# starting at the meadian usually converges
	#

	w <- 1/(o.s2/o.n + sb2)
	y.tilde <- sum(w*o.mu)/sum(w)
	if (trace) print(c("w:", w))
	mu.e <- sum(w*o.mu)/sum(w)
	sd.e <- sqrt(sum(w^2*(o.mu - y.tilde)^2))/sum(w)
	ci.e <- mu.e + qnorm(c(alpha/2, 1 - alpha/2))*sd.e

	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, sb2 = sb2, iter = k,
	             sug = length(o.mu) > 5, wi = w) )
}


#=================================
#
# Graybill-Deal Method
#
# graybill.deal(d[,1], d[,2]^2, d[,3])
# 2017-12-13.
#   returns the weights.
# 2018-07-14.
#   accepts the type B uncertainty, by default it is zero.
#
#=================================
.internal.graybill.deal <- function(o.mu, o.s2, o.n, 
      o.s2.B = rep(0, length(o.mu)), method = "naive", alpha = 0.05) {
	o.u2 <- o.s2/o.n
	o.n[o.n == Inf] <- 1

	w <- o.n/o.s2
	mu.e <- sum(w*o.mu)/sum(w)

	if (method == "naive") {
		w <- o.n/o.s2
		sd.e <- sqrt(1/sum(w))
	}

	if (method == "sinha") {
		if (any(o.n < 2)) 
		  stop("the Sinha method requires 2 readings at least.")
		w <- o.n/o.s2/sum(o.n/o.s2)
		sd.e <- sqrt(1/sum(o.n/o.s2)*(1 + 4*sum(w*(1 - w)/(o.n - 1))))
	}

	if (method == "zhang1") {
		if (any(o.n < 3)) 
		  stop("the Zhang method requires 4 readings at least.")
	  w.i <- 1/((o.n-1)/(o.n-3)*o.s2/o.n + o.s2.B)
		w <- w.i/sum(w.i)
		sd.e <- sqrt(1/sum(w.i))
	}

	if (method == "zhang2") {
		if (any(o.n < 3)) 
		  stop("the Zhang method requires 4 readings at least.")
	  w.i <- 1/((o.n-1)/(o.n-3)*o.s2/o.n + o.s2.B)
		w <- w.i/sum(w.i)
		sd.e <- sqrt(1/sum(w.i))*(1 + 2*sum(w*(1 - w)/(o.n - 1)))
	}

	ci.e <- mu.e + qnorm(c(alpha/2, 1 - alpha/2))*sd.e
	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, iter = 1, w.i = w) )
}

#=================================
#
# DerSimonian-Laird Method
#
# dersimonian.laird(d[,1], d[,2]^2, d[,3])
#
#=================================
.internal.dersimonian.laird <- function(o.mu, o.s2, o.n, alpha = 0.05) {
	gd <- .internal.graybill.deal(o.mu, o.s2, o.n)
	k <- length(o.mu)
	term1 <- sum(o.n*(o.mu - gd$mu)^2/o.s2) - k + 1
	term2 <- sum(o.n/o.s2)
	term3 <- sum((o.n/o.s2)^2)
	Y.DL <- max(0, term1/(term2 - term3/term2))
	w.i <- 1/(Y.DL + o.s2/o.n)
	w <- w.i/sum(w.i)
	x.DL <- sum(w*o.mu)

	# Higgins et al
	sd.DL <- 1/sqrt(sum(w.i))

	# as stated in CONSENSU_KCRV_V10.pdf section 3.4
#	sd.DL <- sqrt(sum(w^2*(o.mu - x.DL)^2/(1 - w)))

	ci.e <- x.DL + qnorm(c(alpha/2, 1 - alpha/2))*sd.DL
	return(list(mu = x.DL, u.mu = sd.DL, ci.mu = ci.e, sb2 = Y.DL, iter = 1,
	            w.i = w.i))
}

.internal.ss.dersimonian.laird <- function(o.mu, o.s2, o.n, alpha = 0.05) {
	gb <- .internal.graybill.deal(o.mu, o.s2, o.n)
	k <- length(o.mu)
	term1 <- sum(o.n*(o.mu - gb$mu)^2/o.s2) - k + 1 
	term2 <- sum(o.n/o.s2)
	term3 <- sum((o.n/o.s2)^2)
	Y.DL <- max(0, term1/(term2 - term3/term2))
	w <- 1 / (Y.DL + o.s2/o.n)
	w <- w / sum(w)

	gb$mu <- sum(w*o.mu)
	term1 <- sum(o.n*(o.mu - gb$mu)^2/o.s2) - k + 1 
	term2 <- sum(o.n/o.s2)
	term3 <- sum((o.n/o.s2)^2)
	Y.DL <- max(0, term1/(term2 - term3/term2))
	w <- 1/(Y.DL + o.s2/o.n)
	w <- w/sum(w)

	gb$mu <- sum(w*o.mu)
	term1 <- sum(o.n*(o.mu - gb$mu)^2/o.s2) - k + 1 
	term2 <- sum(o.n/o.s2)
	term3 <- sum((o.n/o.s2)^2)
	Y.DL <- max(0, term1/(term2 - term3/term2))
	w.i <- 1/(Y.DL + o.s2/o.n)
	w <- w.i/sum(w.i)

	x.DL <- sum(w*o.mu)
	sd.DL <- sqrt(sum(w^2*(o.mu - x.DL)^2/(1 - w)))
	if (length(o.mu)<10) {
		sd.DL <- 1/sqrt(sum(w.i))
	}
	ci.e <- x.DL + qnorm(c(alpha/2, 1 - alpha/2))*sd.DL
	return( list(mu = x.DL, u.mu = sd.DL, ci.mu = ci.e, sb2 = Y.DL, 
	  w.i = w.i) )
}


#=================================
#
# Type B on Bias Method
#
# 2017-03-14: modified to include the use of alpha and expansion type
#
#=================================
.internal.bob <- function(o.mu, o.s2, o.n, alpha = 0.05, expansion = "naive") {
  mm <- length(o.mu)
  kp <- 2
  if (expansion == "naive") {
    kp <- 2
  } else {
    if (expansion == "large sample") {
      kp <- qnorm(1 - alpha/2)
    } else {
      kp <- qt(1 - alpha/2, mm - 1)
    }
  }
  
  mu.e <- mean(o.mu)
	sd.e <- sqrt( sum(o.s2/o.n)/length(o.mu)^2 + (max(o.mu) - min(o.mu))^2/12)
	ci.e <- mu.e + kp*c(-1,1)*sd.e
	return( list(mu = mu.e, u.mu = sd.e, ci.mu = ci.e, iter = 1) )
}



#=================================
#
# Schiller-Eberhardt Method [alpha]
#
# 2017-03-14: modified to include the use of alpha
#
#=================================
.internal.Schiller.Eberhardt <- function(o.mu, o.s2, o.n, alpha = 0.05, 
                               s2.h = 0, df.h = 1) {
#	s2.h, df.h: are estimated independent of the sample information
#
	pm.e <- .internal.pmm(o.mu, o.s2, o.n)
	sb2 <- pm.e$sb2

	W.i <- 1/(o.s2/o.n + sb2)
	w <- W.i/sum(W.i)

	mu.e <- sum(w*o.mu)

	omega.i <- o.n/o.s2
	wi <- omega.i/sum(omega.i)
	s2.x.tilde <- sum(w^2*o.s2/o.n)

	bias.allowance <- max(abs(o.mu - mu.e))

	edof <- sum(w^2*o.s2/o.n + s2.h)^2 /
	  (sum((w^2*o.s2/o.n)^2/(o.n - 1)) + s2.h^2/df.h)

	sd.e3 <- (qt(1 - alpha/2,edof)*sqrt(s2.x.tilde + s2.h) + bias.allowance)
	sd.e5 <- 1/sqrt(sum(1/(sb2 + o.s2/o.n + s2.h)))
	return( list(mu = mu.e, u.mu = sd.e5, U.mu = sd.e3, sb2 = sb2,
	             bias.allowance = bias.allowance, dof = edof, iter = pm.e$iter,
	             w.i = W.i) )
}

.internal.mm.pdf <- function(x, B = 1e5, new.seed = 12345) {
	set.seed(new.seed)
	p <- length(x$mean)
	Xb <- matrix(rnorm(B*p, x$mean, x$sd), B, p, byrow = TRUE )

	return(Xb)
}


.internal.mm.mean <- function(Xb, p, alpha = 0.05) {
  res <- apply(Xb, 1, mean)
  res2 <- apply(Xb, 1, var)
  tau2 <- var(res)
  ci.mu = res[order(res)][round(length(res)*c(alpha/2, 1 - alpha/2))]
	return( list(mu = mean(res), u.mu = sd(res), 
	  ci.mu = ci.mu, var.b = tau2 ) )
}

.internal.mm.linearpool <- function(Xb, p, alpha = 0.05) {
  res <- apply(Xb, 1, mean)
  tau2 <- var(res)
  ci.mu = res[order(res)][round(length(res)*c(alpha/2, 1 - alpha/2))]
  return( list(mu = mean(res), u.mu = sd(c(Xb)), 
               ci.mu = ci.mu, var.b = tau2 ) )
}

.internal.mm.median <- function(Xb, p, alpha = 0.05) {
  # Xb is a matrix(B, p)
#  print( dim(Xb) )
  res <- apply(Xb, 1, median)
  res2 <- apply(Xb, 1, mad)
#  print(length(res2))
#  hist(res2, prob = TRUE, col = 8, nclass = 64)

#  X11()
#  par(mfrow=c(5, 5))
#  for (i in 1:p) {
#    hist(Xb[, i], prob = TRUE, nclass = 64, col = 8)
#  }

  tau2 <- (pi/2/p)*(median(res2^2))
  ci.mu = Xb[order(Xb)][round(length(Xb)*c(alpha/2, 1 - alpha/2))]
  return( list(mu = median(Xb), 
               u.mu = sqrt((sqrt(pi/2/p)*mad(Xb))^2 + tau2),
               var.b = tau2, ci.mu = ci.mu))
}
