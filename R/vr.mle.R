#=========================================================
# ONE-WAY RANDOM EFFECTS MODEL MLE
# 1. USING SCORE ESTIMATING EQUATIONS
# 2. USING THE VANGEL-RUKHIN ITERATIVE ALGORITHM
#
# NIST
# Hugo Gasca-Aragon
# created: May 2009
# last update:
#	Apr 2010
#	The VR algorithm was updated with
#       a Gauss-Newton search for the roots
#       instead of the linear search
#
#	Jun 2010
#	The VR result was extended to include
#	if the convergence criteria was met,
#	if the reduced model is suggested
#	(no evidence of random effect)
#
#	Oct 2010
#	The effective degrees of freedom
#	  of the variance of the mean was 
#	  implemented using Satterthwaite approach
#
#	Oct 2015 @ CENAM
#	The trace with each iteration of the computation is 
#	  optionally provided as part of the result instead 
#	  of printing them out. This is included for validation 
#	  purposes.
#	  When compared to Dataplot one often gets different
#	  results. We have no control over Dataplot criteria
#	  to stop the maximum likelihood search. By imposing 
#	  stronger conditions on this algorithm (i.e., 
#	  max.iter=1000, tol=1e-12), the trace output should 
#	  contain the result provided by Dataplot.
#	  When compared between them vr.mle and mlr.1wre provide
#	  distinct results. The mle.1wre fails to produce an 
#	  estimate withint the parameter space (sigma2>0). While 
#	  the vr.mle fails to return the true maximum, it 
#	  continues the search until the estimated parameters 
#	  are closer at each iteration, passing the true maximum.
#	  Hence a search of the maximum likelihood is included
#	  on the traced results and then it is returned.
#
#=========================================================



.internal.mle.1wre <- function(xi, si2, ni, labi = c(1:length(xi)), 
	max.iter = 200, tol = .Machine$double.eps ^ 0.5, trace = FALSE, 
	init.mu = mean(xi), init.sigma2 = var(xi), lambda = 1) {
#
# parameters
#
# xi=reported sample means, si2=reported sample variances, 
# ni=reported sample sizes, labi=participant labels
#

#	require(MASS)

	# sort the datapoints by labi
	xi <- xi[order(labi)]
	si2 <- si2[order(labi)]
	ni <- ni[order(labi)]
	labi <- labi[order(labi)]

	# remove all datapoints with undefined si2
	xi <- xi[!is.na(si2)]
	ni <- ni[!is.na(si2)]
	labi <- labi[!is.na(si2)]
	si2 <- si2[!is.na(si2)]

	p <- length(xi)
	N <- sum(ni)
	Theta <- matrix(-Inf, max.iter + 1, p + 3)

	# set the initial values for sigma2 and sigmai2
	mu <- init.mu
	sigma2 <- init.sigma2
	sigmai2 <- si2
	t <- 1

	llh <- 0
	llh <- -N/2 * log(2 * pi) - sum((ni - 1)*log(sigmai2))/2 - 
		sum(log(sigmai2 + ni * sigma2))/2 - sum((ni - 1)*si2/sigmai2)/2 -
		sum(ni*(xi - mu)^2/(sigmai2 + ni*sigma2))/2

	Theta[t, ] <- c(mu, sigmai2, sigma2, llh)
	cur.rel.abs.error <- Inf

	while ((cur.rel.abs.error > tol) && (t < max.iter) && all(sigmai2 > 0) &&
		(sigma2 > 0)) {

		mu <- Theta[t,1]
		sigmai2 <- Theta[t, c(2:(p + 1))]
		sigma2 <- Theta[t, p + 2]
		wi <- ni/(sigmai2 + ni*sigma2)
		var.mu <- 1/sum(wi)

		llh<- -N/2*log(2*pi) - sum((ni - 1)*log(sigmai2))/2 - 
			sum(log(sigmai2 + ni*sigma2))/2 - sum((ni-1)*si2/sigmai2)/2 - 
			sum(ni*(xi - mu)^2/(sigmai2 + ni*sigma2))/2

		# get the score vector
		S<- c(sum((xi - mu)/(sigma2 + sigmai2/ni)), 
			-(ni - 1)/sigmai2/2 - 1/(ni*sigma2 + sigmai2)/2 + 
				(ni - 1)*si2/sigmai2^2/2 + 
				ni*(xi - mu)^2/(ni*sigma2 + sigmai2)^2/2,
			-sum(1/(sigma2 + sigmai2/ni))/2 + 
				sum((xi - mu)^2/(sigma2 + sigmai2/ni)^2)/2)

		# get the information matrix
		I <- matrix(0, p + 2, p + 2)

		I[1, 1] <- sum(1/(sigma2 + sigmai2/ni))
		for (j in 1:p) {
			I[j + 1, j + 1] <- (1/2)*(ni[j] - 1)/sigmai2[j]^2 + 
				1/(sigmai2[j] + ni[j]*sigma2)^2/2
			I[j + 1, p + 2] <- ni[j]/2/(sigmai2[j]+ni[j]*sigma2)^2
			I[p + 2, j + 1] <- I[j+1,p+2]
		}
		I[p + 2, p + 2] <- sum(1/(sigma2 + sigmai2/ni)^2)/2

		# get the inverse of the information matrix
		Iinv <- ginv(I)

		# update the estimates
		Theta[t + 1, 1:(p + 2)] <- Theta[t, 1:(p + 2)] + lambda*(Iinv %*% S)
		Theta[t + 1, p + 3] <- llh

		# if the new estimate of sigma is negative set to zero
		# compute the maximum relative absolute error
		# measure the relative error based on the estimates
		# the values of the used weigths and sigma2 are relative, the weigths 
		# and sigma2 keep moving
		old.theta <- Theta[t, ]
		new.theta <- Theta[t+1, ]
		map.old.theta <- old.theta
		map.old.theta[2:(p + 1)] <- 1/(old.theta[2:(p + 1)]/ni +
			old.theta[p + 2])
		map.old.theta[p + 2] <- 1/sum(map.old.theta[2:(p + 1)])

		map.new.theta <- new.theta
		map.new.theta[2:(p + 1)] <- 1/(new.theta[2:(p + 1)]/ni + 
			new.theta[p + 2])
		map.new.theta[p + 2] <- 1/sum(map.new.theta[2:(p + 1)])

		cur.rel.abs.error <- max(abs((map.old.theta - 
			map.new.theta)/map.new.theta))

		t <- t + 1

		mu <- Theta[t, 1]
		sigmai2 <- Theta[t, c(2:(p + 1))]
		sigma2 <- Theta[t, p + 2]

		if (is.na(sigma2)) stop("sigma2 became undefined.")
		if (any(is.na(sigmai2))) stop("some sigmai2 became undefined")
		if (is.na(mu)) stop("mu became undefined")
		if (is.na(cur.rel.abs.error)) 
			stop("current relative absolute error became undefined")
	} # while

	if ((t == max.iter) || (cur.rel.abs.error > tol) || any(sigmai2 <= 0) ||
		(sigma2 <= 0)) {
		warning("Non convergence or slow convergence condition was found.")
	} 

	Theta <- Theta[1:t,]
	mu <- Theta[t, 1]
	sigmai2 <- Theta[t, c(2:(p + 1))]
	sigma2 <- Theta[t, p + 2]
	wi <- ni/(sigmai2 + ni*sigma2)
	var.mu <- 1/sum(wi)
	llh <- -N/2*log(2*pi) - sum((ni - 1)*log(sigmai2))/2 - 
		sum(log(sigmai2 + ni*sigma2))/2 - sum((ni - 1)*si2/sigmai2)/2 - 
		sum(ni*(xi - mu)^2/(sigmai2 + ni*sigma2))/2

	if (!trace) Theta <-NULL
	result<- list( mu = as.vector(mu), 
		var.mu = as.vector(var.mu), 
		sigma2 = as.vector(sigma2), 
		llh = as.vector(llh), 
		tot.iter = as.vector(t), 
		max.rel.abs.error = as.vector(cur.rel.abs.error), 
		sigmai2 = as.vector(sigmai2),
		trace = Theta
	)
	class(result)<- "summary.mle.1wre"

	return(result)

} # function block



.internal.newton.raphson <- function(f, fp, init.value, 
		max.tol = .Machine$double.eps^0.5, trace = FALSE, range.tol = 0) {
	max.iter <- 40
	tol <- 1
	iter <- 0
	x <- init.value
	if (fp(x) != 0) {
		while ((iter < max.iter) & (tol > max.tol)) {
			xn <- x - f(x)/fp(x)
			if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
			iter <- iter + 1

			if (range.tol == 0) {
				tol <- abs(xn - x)
			} else 
			if (range.tol == 1) {
				tol <- abs(xn - x)/abs(max(x,xn))
			} else {
				tol <- abs(f(xn))
			}

			x <- xn
			if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
		}
		if (x > 1) x <- 1
		if (x < 0) x <- 0
		if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
	} else {
		x <- init.value
	}
	return(x)
}


.internal.find.roots <- function(mu, sigma2, gammai, xi, si2, ni, 
		tol = .Machine$double.eps^0.25, trace = FALSE)
{
	p <- length(xi)
	N <- sum(ni)
	wi <- gammai/sigma2
	sigmai2 <- ni*sigma2*(1 - gammai)/gammai
	var.mu <- 1/sum(wi)
	nroots <- rep(NA, p)

	llh <- -N/2*log(2*pi) - sum((ni - 1)*log(sigmai2))/2 - 
		sum(log(sigmai2 + ni*sigma2))/2 - sum((ni - 1)*si2/sigmai2)/2 - 
		sum(ni*(xi - mu)^2/(sigmai2 + ni*sigma2))/2

	if (trace) print( c(t, mu, var.mu, sigma2, sigmai2, llh) )

	# evaluate the maximum likelihood function on the estimated parameters
	s2 <- sigma2 

	# first find the solutions to the variance components
	nj <-1
	xj <-1
	sj2 <-1

	ai <- sigma2/(xi - mu)^2
	bi <- si2/(ni*(xi - mu)^2)
	bb <- -(ai + 2)
	cc <- ((ni + 1)*ai + (ni - 1)*bi + 1)
	dd <- -ni*ai

	pol <- function(x) {
		return( x^3 + bb*x^2 + cc*x + dd )
	}

	roots <- c(1:p)
	vroots <- matrix(NA, p, 3)

	if (s2 == 0) {
		# the polynomial reduces to a monomial and we have straight solutions.
		roots <- (xi - mu)^2 + (ni - 1)*si2/ni
	} else {

		delta <- roots
		delta.w <- roots
		if (trace) print( "list of positive roots" )
		for (i in 1:p) {
			# print( c("searching for component variance ", i) )
			nj <- ni[i]
			xj <- xi[i]
			sj2 <- si2[i]

			aj <- sigma2/(xj - mu)^2
			bj <- sj2/(nj*(xj - mu)^2)

			#	let's find how many roots the polynomial has.
			a <- (2*nj - 1)*s2 - (xj - mu)^2 - (nj - 1)*sj2/nj
			b <- (nj - 1)*s2*(nj*s2 - 2*sj2)
			c <- -(nj - 1)*nj*sj2*(s2^2)
			delta[i] <- -4*a^3*c + a^2*b^2 - 4*b^3 + 18*a*b*c - 27*c^2

			bb <- -(aj + 2)
			cc <- ((nj + 1)*aj+(nj - 1)*bj + 1)
			dd <- -nj*aj
			delta.w[i] <- -18*bb*cc*dd + 4*bb^3*dd + bb^2*cc^2 - 4*cc^3 +
				27*dd^2
	
			if (delta.w[i] < 0) {
				# there is one single real positive root
				# print("looking for one root")
				j <- 1
				root <- rep(0, j)

				root[1] <- uniroot( pol, interval=c(0, 1), tol = tol )$root

				if (trace) print( c(i, root) )
				vroots[i, 1] <- root[1]
				nroots[i] <- 1

				if (root[1] <= 0) {
					stop( c("single negative root was found for dataset ", i) )
				}

				pllh <- rep(0, nroots[i])

				for (j in 1:nroots[i])
					pllh[j] <- -nj/2*log(2*pi) - (nj - 1)*log(root[j])/2 - 
						log(root[j] + nj*s2)/2 - (nj - 1)*sj2/root[j]/2 - 
						nj*(xj - mu)^2/(root[j] + nj*s2)/2
				if (trace) print( c(i, pllh) )
				if (trace) print( c(i, roots[i]) )

			} else if (delta[i] == 0) {
				# there are at most 2 different real roots

				#=============================================
				# finding all the real roots
				# first obtain the local maximum and local minimum of the polynomial
				# by finding the roots of the derivate of the polynomial
				# then look for the roots in the partition formed by (0, root1, root2, Inf)
				# finally choose the root that maximize the likelihood function
				#=============================================

				lim.root1 <- (-2*a + sqrt(4*a^2 - 12*b))/6
				lim.root2 <- (-2*a - sqrt(4*a^2 - 12*b))/6
				lim <- c(1:4)

				lim[1] <- min(lim.root1, lim.root2)
				lim[2] <- max(lim.root1, lim.root2)

				if (lim[1] > 0) {
					j <- 0
					lim[3] <- lim[1] - 1
					while (pol(lim[3] - 10^j) > 0) j <- j + 1
					lim[3] <- lim[3] - 10^j
				} else {
					lim[3] <- 0
				}
				j <-0
				lim[4] <- lim[2] + 1
				while (pol(lim[4] + 10^j) < 0) j <- j + 1
				lim[4] <- lim[4] + 10^j

				lim<- sort(lim)

				if (lim[1] == lim[2]) lim <- lim[-1]

				root <- c(1:2)
				root[1] <- uniroot( pol, interval = c(lim[1], lim[2]), 
					tol=tol)$root
				root[2] <- uniroot( pol, interval = c(lim[2], lim[3]), 
					tol=tol)$root

				if (trace) print( c(i, root) )

				# there are 2 roots but some roots can be negative
				if (root[1] <= 0) {
					root <- root[2]
				}

				if (trace) print( c(i, root) )

				nroots[i] <- length(root)
				pllh <- rep(0, nroots[i])

				for (j in 1:nroots[i])
					pllh[j] <- -nj/2*log(2*pi) - (nj - 1)*log(root[j])/2 - 
						log(root[j] + nj*s2)/2 - (nj - 1)*sj2/root[j]/2 - 
						nj*(xj - mu)^2/(root[j] + nj*s2)/2

				if (any(is.na(pllh))) {
					print(pllh)
					print(c(N, mu, s2, root, xj, sj2, nj))
					stop("Abnormal condition found")
				}

				for (j in 1:nroots[i])
					if (pllh[j] == max(pllh)) roots[i] <- root[j]

				if (trace) print( c(i, pllh) )
				if (trace) print( c(i, roots[i]) )
			} else {
				#=============================================
				# finding all the real roots
				# first obtain the local maximum and local minimum of the polynomial
				# by finding the roots of the derivate of the polynomial
				# then look for the roots in the partition formed by (0, root1, root2, Inf)
				# finally choose the root that maximize the likelihood function
				#=============================================
				lim.root1 <- (-2*a + sqrt(4*a^2 - 12*b))/6
				lim.root2 <- (-2*a - sqrt(4*a^2 - 12*b))/6

				lim <- c(1:4)
				lim[1] <- min(lim.root1, lim.root2)
				lim[2] <- max(lim.root1, lim.root2)

				j <-0
				lim[3] <- lim[1] - 1
				while (pol(lim[3] - 10^j)>0) j <- j + 1
				lim[3] <- lim[3] - 10^j

				j <-0
				lim[4] <- lim[2] + 1
				while (pol(lim[4] + 10^j)<0) j <- j + 1
				lim[4] <- lim[4] + 10^j

				lim <- sort(lim)

				root <- c(1:3)
				root[1] <- uniroot( pol, interval = c(lim[1],lim[2]), 
					tol = tol)$root
				root[2] <- uniroot( pol, interval = c(lim[2],lim[3]), 
					tol = tol)$root
				root[3] <- uniroot( pol, interval = c(lim[3],lim[4]), 
					tol = tol)$root

				if (trace) print( c(i, root) )

				# there are 3 roots but some roots can be negative
				if (root[1] <= 0) {
					root <- root[2:3]
				}
				if (root[1] <= 0) {
					root <- root[2]
				}


				if (trace) print( c(i, root) )
				nroots[i] <- length(root)
				pllh <- rep(0, nroots[i])

				for (j in 1:nroots[i])
					pllh[j] <- -nj/2*log(2*pi) - (nj - 1)*log(root[j])/2 - 
						log(root[j] + nj*s2)/2 - (nj - 1)*sj2/root[j]/2 - 
						nj*(xj - mu)^2/(root[j] + nj*s2)/2

				if (any(is.na(pllh))) {
					print(pllh)
					print(c(N, mu, s2, root, xj, sj2, nj))
					stop("Abnormal condition found")
				}

				for (j in 1:nroots[i])
					if (pllh[j] == max(pllh, na.rm = TRUE)) roots[i] <- root[j]

				if (trace) print( c(i, pllh) )
				if (trace) print( c(i, roots[i]) )
			}
		} # for i in 1:p
	} # look for the roots that maximize the log likelihood function

	return( c(roots, nroots) )
}






vr.mle <- function(xi, si2, ni, labi = c(1:length(xi)), max.iter = 1000, 
	tol = .Machine$double.eps^0.5, init.mu = mean(xi), init.sigma2 = var(xi), 
	trace = FALSE, alpha = 0.05)
{
#
# parameters
#
# xi=reported sample means, si2=reported sample variances, ni=reported sample sizes, labi=participant labels
#

	# sort the datapoints by labi
	xi <- xi[order(labi)]
	si2 <- si2[order(labi)]
	ni <- ni[order(labi)]
	labi <- labi[order(labi)]

	# remove all datapoints with undefined si2
	xi <- xi[!is.na(si2)]
	ni <- ni[!is.na(si2)]
	labi <- labi[!is.na(si2)]
	si2 <- si2[!is.na(si2)]

	p <- length(xi)
	if (p <= 1) { stop("vr.mle requires 2 or more sources of information.") }
	N <- sum(ni)
	Theta <- matrix(-Inf, max.iter + 1, p + 3)
	Delta.theta <- matrix(0, max.iter + 1, p)

	# set the initial values for sigma2, sigmai2 and gammai
	mu <- init.mu
	sigma2 <- init.sigma2
	sigmai2 <- si2
	gammai <- sigma2/(sigma2 + sigmai2/ni)

	t <- 1

	var.mu <- sigma2/sum(gammai)
	llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2+sum(ni*log(gammai/sigma2))/2 - 
		sum((ni - 1)*log(1-gammai))/2 - sum(gammai*((xi - mu)^2 + 
		(ni - 1)*si2/ni/(1 - gammai)))/sigma2/2
	Theta[t, ] <- c(mu, gammai, sigma2, llh)
	cur.rel.abs.error <- Inf

	while ((cur.rel.abs.error > tol) && (t < max.iter) && all(gammai > 0) && 
			(sigma2 > 0)) {
		mu <- Theta[t, 1]
		gammai <- Theta[t, c(2:(p + 1))]
		sigma2 <- Theta[t, p + 2]
		var.mu <- sigma2/sum(gammai)

		llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2 + 
			sum(ni*log(gammai/sigma2))/2 - sum((ni - 1)*log(1 - gammai))/2 - 
			sum(gammai*((xi - mu)^2 + (ni - 1)*si2/ni/(1 - gammai)))/sigma2/2

		# Iterative update
		ai <- sigma2/(xi - mu)^2
		bi <- si2/ni/(xi - mu)^2

		aa <- rep(1, p)
		bb <- -(ai + 2)
		cc <- ((ni + 1)*ai + (ni - 1)*bi + 1)
		dd <- -ni*ai

		Delta <- -(18*aa*bb*cc*dd - 4*bb^3*dd + bb^2*cc^2 - 4*aa*cc^3 - 
			27*aa^2*dd^2)
		Delta.theta[t, ] <- Delta
		ss <- Delta.theta[t,] <= 0
		if (any(ss)) {
			# there are 3 positive roots 
			# search for the extreme points as the roots of p'=3*aa*x^2+2*bb*x+cc
			sss <- (bb^2 - 3*cc)[ss] >= 0
			if (any(sss)) {
				# there are real roots hence there are two local extreme points
				# and we can find their location as
				gamma.lim <- matrix(NA, sum(sss), 2)
				for(iii in c(1:length(sss))[sss]) {
					gamma.lim[iii, ] <- sort((-bb[ss][sss][iii] + 
						c(-1,1)*sqrt( (bb^2 - 3*cc)[ss][sss][iii] ))/3)
					if (min(gamma.lim[iii, ]) > 1) {
						# there is nothing to do in addition, the solution is out of the parameter space
					} else {
						# there is at least one additional point to evaluate the MLE
						# we need to evaluate the MLE for each combination of these values
						# and select the combination that maximized the response.
						warning(paste("There are additional points ",
							"where the MLE was not evaluated."))
					}
				}
			}
		}

		# solve for weights gammai, this just looks for the closer root, it does not evaluate the three posible solutions
		for (i in 1:p) {
			# this just looks for the closer solution.
			nlm.res <- .internal.newton.raphson(function(x) {
					x^3 - (ai[i] + 2)*x^2 + 
					((ni[i] + 1)*ai[i] + (ni[i] - 1)*bi[i] + 1)*x - 
					ni[i]*ai[i]}, 
				function(x) {3*x^2 - 2*(ai[i] + 2)*x + (ni[i] + 1)*ai[i] + 
					(ni[i] - 1)*bi[i] + 1}, 
				init.value = Theta[t, 1 + i], max.tol = tol) 

			# we must look at all the possible solutions.
			# we need to find all the roots of each parameter, then evaluate the combination that leads to the maximum log likelihood.

			gammai[i] <- nlm.res #$estimate
			if (gammai[i] > 1) stop("gamma[", i, "]>1")
			if (gammai[i] < 0) stop("gamma[", i, "]<0")
		}
		
##		increased accuracy, reducing error due to rounding and double underflow effects
		mu <- sum(gammai/sum(gammai)*xi)
		sigma2 <- sum(gammai/sum(gammai)*gammai*(xi - mu)^2)

		Theta[t+1, ] <- c(mu, gammai, sigma2, llh)

		# compute the maximum relative absolute error

		# measure the relative error based on the estimates
		# the values of the used weigths and sigma2 are relative, the weigths and sigma2 keep moving
		old.theta <- Theta[t, ]
		new.theta <- Theta[t + 1, ]
		map.old.theta <- old.theta
		map.old.theta[2:(p + 1)] <- old.theta[2:(p + 1)]/
			sum(old.theta[2:(p + 1)])
		map.old.theta[p + 2] <- old.theta[p + 2]/sum(old.theta[2:(p + 1)])
		map.new.theta <- new.theta
		map.new.theta[2:(p + 1)] <- new.theta[2:(p + 1)]/
			sum(new.theta[2:(p + 1)])
		map.new.theta[p + 2] <- new.theta[p + 2]/sum(new.theta[2:(p + 1)])

		cur.rel.abs.error<- max(abs((map.old.theta - map.new.theta)/
			map.new.theta))

		t <- t + 1

		mu <- Theta[t, 1]
		gammai <- Theta[t, c(2:(p + 1))]
		sigma2 <- Theta[t, p + 2]

		if (is.na(sigma2)) stop("sigma2 became undefined.")
		if (any(is.na(gammai))) stop("some gammai became undefined")
		if (is.na(mu)) stop("mu became undefined")
		if (is.na(cur.rel.abs.error)) 
			stop("current relative absolute error became undefined")
	} # while

	ccm<- TRUE
	if ((t == max.iter)||(cur.rel.abs.error > tol)|| any(gammai <= 0) || 
		(sigma2 <= 0)) {
		# convergence condition is not met
		ccm <- FALSE
	} 

	Theta <- Theta[1:t,]
	t <- max(c(1:t)[Theta[, p + 3] == max(Theta[, p + 3])])
	Theta <- Theta[1:t, ]
	Delta.theta <- Delta.theta[1:t, ]
	mu <- Theta[t, 1]
	gammai <- Theta[t, c(2:(p + 1))]
	sigma2 <- Theta[t, p + 2]
	var.mu <- sigma2/sum(gammai)

	old.theta <- Theta[t - 1, ]
	new.theta <- Theta[t, ]
	map.old.theta <- old.theta
	map.old.theta[2:(p + 1)] <- old.theta[2:(p + 1)]/sum(old.theta[2:(p + 1)])
	map.old.theta[p + 2] <- old.theta[p + 2]/sum(old.theta[2:(p + 1)])
	map.new.theta <- new.theta
	map.new.theta[2:(p + 1)] <- new.theta[2:(p + 1)]/sum(new.theta[2:(p + 1)])
	map.new.theta[p + 2] <- new.theta[p + 2]/sum(new.theta[2:(p + 1)])

	cur.rel.abs.error <- max(abs((map.old.theta - map.new.theta)/
		map.new.theta))

	# iterative reweighted mean if sigma2==0
	if (sigma2 <= 0) { 
		sigma2 <- 0
		si <- sqrt(si2)
		sigmai <- si
		iter <- 0
		tol <- 1
		mu <- mean(xi)
		max.tol <- .Machine$double.eps^0.5
		while ((iter < 100) & (tol > max.tol)) {
			mu.n <- sum(xi/sigmai^2)/sum(1/sigmai^2)
			sigmai.n <- sqrt(ni*(xi - mu)^2 + (ni - 1)*si^2)/ni
			tol <- abs(mu.n - mu)
			mu <- mu.n
			sigmai <- sigmai.n
			iter <- iter + 1
		}
		wi <- 1/sigmai^2
		u.mu <- 1/sqrt(sum(wi))
		N <- sum(ni)
		llh <- -N*log(2*pi)/2 - sum(log(sigmai^2))/2 - 
			sum((ni*(xi - mu)^2 + (ni - 1)*si^2)/(2*sigmai^2))
		var.mu <- u.mu^2

		gammai <- wi
	} else {
		# likelihood evaluated at the estimated parameters
		wi <- gammai/sigma2
		llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2 + 
			sum(ni*log(gammai/sigma2))/2 - sum((ni - 1)*log(1-gammai))/2 - 
			sum(gammai*((xi - mu)^2 + (ni - 1)*si2/ni/(1 - gammai)))/sigma2/2
	}

	if (sigma2 > 0) {
		# Welch-Satterthwaite equation for the One-Way Random Effects Model for the mean MSA/nT
		v <- sum(wi)^2/sum(wi^2/(ni - 1))

		# Welch-Satterthwaite equation for the One-Way Random Effects Model for the within variance MSE
		sigma2.wi <- (1/gammai - 1)*sigma2*ni
		v.w <- sum((ni - 1)*sigma2.wi)^2/sum((ni - 1)*sigma2.wi^2)

		var.w <- (sum(gammai)^2*var.mu - sum(gammai^2)*sigma2)/sum(gammai^2/ni)
		v.w <- (sum(ni)*var.mu - var.w)/sigma2 - 1
	} else {
		v <- N - 1
		v.w <- v

		var.w <- (sum(gammai)^2*var.mu)/sum(gammai^2/ni)
		v.w <- var.w/var.mu*sum(gammai^2)/sum(gammai)^2
	}

	if (!trace) {
		Theta<- NULL
	}
	result <- list( mu = as.vector(mu), 
     		u.mu = as.vector(sqrt(var.mu)),
     		ci.mu = as.vector(mu + qnorm(alpha/2, 1 - alpha/2)*sqrt(var.mu)),
		var.mu = as.vector(var.mu), 
		var.b = as.vector(sigma2), 
		var.w = as.vector(var.w),
		dof.w = v.w,
		llh = as.vector(llh), 
		tot.iter = as.vector(t), 
		max.rel.abs.error = as.vector(cur.rel.abs.error), 
		gammai = as.vector(gammai), 
		ccm = ccm, 
		reduced.model = (sigma2 == 0), 
		dof = v,
		trace = Theta,
		discriminant = Delta.theta
	)

	class(result) <- "summary.vr.mle"

	return(result)

} # function block


.internal.vr.mle.fixed <- function(xi, si2, ni, labi = c(1:length(xi)), 
	max.iter = 1000, tol = .Machine$double.eps^0.5, fixed.mu = mean(xi), 
	init.sigma2 = var(xi), trace = FALSE, alpha = 0.05)
{
#
# parameters
#
# xi=reported sample means, si2=reported sample variances, ni=reported sample sizes, labi=participant labels
#

	# sort the datapoints by labi
	xi <- xi[order(labi)]
	si2 <- si2[order(labi)]
	ni <- ni[order(labi)]
	labi <- labi[order(labi)]

	# remove all datapoints with undefined si2
	xi <- xi[!is.na(si2)]
	ni <- ni[!is.na(si2)]
	labi <- labi[!is.na(si2)]
	si2 <- si2[!is.na(si2)]

	p <- length(xi)
	if (p <= 1) stop("vr.mle requires 2 or more sources of information.")
	N <- sum(ni)
	Theta <- matrix(-Inf, max.iter + 1, p + 2)
	Delta.theta <- matrix(0, max.iter + 1, p)

	# set the initial values for sigma2, sigmai2 and gammai
	mu <- fixed.mu
	sigma2 <- init.sigma2
	sigmai2 <- si2
	gammai <- sigma2/(sigma2 + sigmai2/ni)

	t <- 1

	var.mu <- sigma2/sum(gammai)
	llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2 + sum(ni*log(gammai/sigma2))/2 - 
		sum((ni - 1)*log(1 - gammai))/2 - 
		sum(gammai*((xi - mu)^2 + (ni - 1)*si2/ni/(1 - gammai)))/sigma2/2
	Theta[t, ] <- c(gammai, sigma2, llh)
	cur.rel.abs.error <- Inf

	while ((cur.rel.abs.error > tol) && (t < max.iter) && all(gammai > 0) && 
			(sigma2 > 0)) {
		gammai <- Theta[t, c(1:p)]
		sigma2 <- Theta[t, p + 1]
		var.mu <- sigma2/sum(gammai)

		llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2 + sum(ni*log(gammai/sigma2))/2 - 
			sum((ni - 1)*log(1 - gammai))/2 - 
			sum(gammai*((xi - mu)^2 + (ni - 1)*si2/ni/(1 - gammai)))/sigma2/2

		# Iterative update
		ai <- sigma2/(xi - mu)^2
		bi <- si2/ni/(xi - mu)^2

		aa <- rep(1, p)
		bb <- -(ai + 2)
		cc <- ((ni + 1)*ai + (ni - 1)*bi + 1)
		dd <- -ni*ai

		Delta <- -(18*aa*bb*cc*dd - 4*bb^3*dd + bb^2*cc^2 - 4*aa*cc^3 - 
			27*aa^2*dd^2)
		Delta.theta[t, ] <- Delta
		ss <- Delta.theta[t, ] <= 0
		if (any(ss)) {
			# there are 3 positive roots 
			# search for the extreme points as the roots of p'=3*aa*x^2+2*bb*x+cc
			sss <- (bb^2 - 3*cc)[ss] >= 0
			if (any(sss)) {
				# there are real roots hence there are two local extreme points
				# and we can find their location as
				gamma.lim <- matrix(NA, sum(sss), 2)
				for(iii in c(1:length(sss))[sss]) {
					gamma.lim[iii, ] <- sort((-bb[ss][sss][iii] + 
						c(-1, 1)*sqrt( (bb^2 - 3*cc)[ss][sss][iii] ))/3)
					if (min(gamma.lim[iii, ]) > 1) {
						# there is nothing to do in addition, the solution is out of the parameter space
					} else {
						# there is at least one additional point to evaluate the MLE
						# we need to evaluate the MLE for each combination of these values
						# and select the combination that maximized the response.
						warning(paste("There are additional points where ",
							"the MLE was not evaluated."))
					}
				}
			}
		}

		# solve for weights gammai, this just looks for the closer root, it does not evaluate the three posible solutions
		for (i in 1:p) {
			# this just looks for the closer solution.
			nlm.res <- .internal.newton.raphson(function(x) {x^3 - 
					(ai[i] + 2)*x^2 + 
					((ni[i] + 1)*ai[i] + (ni[i] - 1)*bi[i] + 1)*x - 
					ni[i]*ai[i]}, 
				function(x) {3*x^2 - 2*(ai[i] + 2)*x + 
					(ni[i] + 1)*ai[i] + (ni[i] - 1)*bi[i] + 1}, 
				init.value = Theta[t, 1 + i], max.tol = tol) 

			# we must look at all the possible solutions.
			# we need to find all the roots of each parameter, then evaluate the combination that leads to the maximum log likelihood.

			# nlm(f=function(x) {x^4/4-(ai[i]+2)*x^3/3+((ni[i]+1)*ai[i]+(ni[i]-1)*bi[i]+1)*x^2/2-ni[i]*ai[i]*x}, p=c(Theta[t,1+i]), steptol=tol)
			gammai[i] <- nlm.res #$estimate
			if (gammai[i] > 1) stop("gamma[",i,"]>1")
			if (gammai[i] < 0) stop("gamma[",i,"]<0")
		}
		
##		increased accuracy, reducing error due to rounding and double underflow effects
#		mu<- sum(gammai/sum(gammai)*xi)
		sigma2 <- sum(gammai/sum(gammai)*gammai*(xi - mu)^2)

		Theta[t+1, ] <- c(gammai, sigma2, llh)

		# compute the maximum relative absolute error
		# measure the relative error based on the estimates
		# the values of the used weigths and sigma2 are relative, the weigths and sigma2 keep moving
		old.theta <- Theta[t, ]
		new.theta <- Theta[t + 1, ]
		map.old.theta <- old.theta
		map.old.theta[1:p] <- old.theta[1:p]/sum(old.theta[1:p])
		map.old.theta[p + 1] <- old.theta[p + 1]/sum(old.theta[1:p])
		map.new.theta <- new.theta
		map.new.theta[1:p] <- new.theta[1:p]/sum(new.theta[1:p])
		map.new.theta[p + 1] <- new.theta[p + 1]/sum(new.theta[1:p])

		cur.rel.abs.error <- max(abs((map.old.theta - map.new.theta)/
			map.new.theta))

		t <- t + 1
		gammai <- Theta[t, c(1:p)]
		sigma2 <- Theta[t, p + 1]

		if (is.na(sigma2)) stop("sigma2 became undefined.")
		if (any(is.na(gammai))) stop("some gammai became undefined")
		if (is.na(mu)) stop("mu became undefined")
		if (is.na(cur.rel.abs.error)) 
			stop("current relative absolute error became undefined")
	} # while

	ccm <- TRUE
	if ((t == max.iter) && (cur.rel.abs.error > tol)) {
		# convergence condition is not met
		ccm <- FALSE
	} 

	Theta <- Theta[1:t, ]
	t <- max(c(1:t)[Theta[, p + 2] == max(Theta[, p + 2])])
	Theta <- as.matrix(Theta[1:t, ], t, p + 2)
	Delta.theta <- Delta.theta[1:t, ]

	gammai <- Theta[t, c(1:p)]
	sigma2 <- Theta[t, p + 1]
	var.mu <- sigma2/sum(gammai)

	old.theta <- Theta[t - 1, ]
	new.theta <- Theta[t, ]
	map.old.theta <- old.theta
	map.old.theta[1:p] <- old.theta[1:p]/sum(old.theta[1:p])
	map.old.theta[p + 1] <- old.theta[p + 1]/sum(old.theta[1:p])
	map.new.theta <- new.theta
	map.new.theta[1:p] <- new.theta[1:p]/sum(new.theta[1:p])
	map.new.theta[p + 1] <- new.theta[p + 1]/sum(new.theta[1:p])

	cur.rel.abs.error <- max(abs((map.old.theta - map.new.theta)/
		map.new.theta))

	# iterative reweighted mean if sigma2==0
	if (sigma2 == 0) { 
		sigma2 <- 0
		si <- sqrt(si2)
		sigmai <- si
		iter <- 0
		tol <- 1
		mu <- mean(xi)
		max.tol <- .Machine$double.eps^0.5
		while ((iter < 100) & (tol > max.tol)) {
			mu.n <- sum(xi/sigmai^2)/sum(1/sigmai^2)
			sigmai.n <- sqrt(ni*(xi - mu)^2 + (ni - 1)*si^2)/ni
			tol <- abs(mu.n - mu)
			mu <- mu.n
			sigmai <- sigmai.n
			iter <- iter + 1
		}
		wi <- 1/sigmai^2
		u.mu <- 1/sqrt(sum(wi))
		N <- sum(ni)
		llh <- -N*log(2*pi)/2 - sum(log(sigmai^2))/2 - 
			sum((ni*(xi - mu)^2 + (ni - 1)*si^2)/(2*sigmai^2))
		var.mu <- u.mu^2

		gammai <- wi
	} else {
		# likelihood evaluated at the estimated parameters
		wi <- gammai/sigma2
		llh <- -N/2*log(2*pi) - sum(ni*log(ni))/2 + 
			sum(ni*log(gammai/sigma2))/2 - sum((ni - 1)*log(1 - gammai))/2 - 
			sum(gammai*((xi - mu)^2 + (ni - 1)*si2/ni/(1 - gammai)))/sigma2/2
	}

	if (sigma2 > 0) {
		# Welch-Satterthwaite equation for the One-Way Random Effects Model for the mean MSA/nT
		v <- sum(wi)^2/sum(wi^2/(ni - 1))

		# Welch-Satterthwaite equation for the One-Way Random Effects Model for the within variance MSE
		sigma2.wi <- (1/gammai - 1)*sigma2*ni
		v.w <- sum((ni - 1)*sigma2.wi)^2/sum((ni - 1)*sigma2.wi^2)

		var.w <- (sum(gammai)^2*var.mu - sum(gammai^2)*sigma2)/sum(gammai^2/ni)
		v.w <- (sum(ni)*var.mu - var.w)/sigma2 - 1
	} else {
		v <- N - 1
		v.w <- v

		var.w <- (sum(gammai)^2*var.mu)/sum(gammai^2/ni)
		v.w <- var.w/var.mu*sum(gammai^2)/sum(gammai)^2
	}

	if (!trace) {
		Theta <- NULL
	}
	result <- list( mu = as.vector(mu), 
     		u.mu = as.vector(sqrt(var.mu)),
     		ci.mu = as.vector(mu + qnorm(alpha/2, 1 - alpha/2)*sqrt(var.mu)),
		var.mu = as.vector(var.mu), 
		var.b = as.vector(sigma2), 
		var.w = as.vector(var.w),
		dof.w = v.w,
		llh = as.vector(llh), 
		tot.iter = as.vector(t), 
		max.rel.abs.error = as.vector(cur.rel.abs.error), 
		gammai = as.vector(gammai), 
		ccm = ccm, 
		reduced.model = (sigma2 == 0), 
		dof = v,
		trace = Theta,
		discriminant = Delta.theta
	)

	class(result) <- "summary.vr.mle.fixed"

	return(result)

} # function block
