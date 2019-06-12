OMiSALN <-
function(obstime, delta, X, total.reads = NULL, tree, cov = NULL, pow = c(1/4, 1/3, 1/2, 1), n.perm = 5000) {
	if (is.null(cov)) {
		cov <- as.matrix(rep(1,length(obstime)))
	}
	
	r <- coxph(Surv(obstime, delta) ~ ., data=as.data.frame(cov))$residuals
	r.s <- list()
	for (j in 1:n.perm) {
		r.s[[j]] <- r[shuffle(length(r))]
	}
	
	if (is.null(total.reads)) {
		X.trans.misaln <- lapply(as.list(pow), function(y) apply(t(apply(X,1,function(x) (x/sum(x))^y)),2,scale))
	} else {
		X.trans.misaln <- lapply(as.list(pow), function(y) apply((X/total.reads)^y,2,scale))
	}
	R.misaln <- lapply(X.trans.misaln, function(y) y%*%t(y))

	U.func <- function(y) t(r)%*%y%*%r
	Ts.misaln <- unlist(lapply(R.misaln, U.func))
	
	U0.func <- function(y) unlist(lapply(r.s, function(z) t(z)%*%y%*%z))
	T0s.misaln <- lapply(R.misaln,U0.func)
	
	pvs.misaln <- rep(NA, length(pow))
	for (i in 1:length(pow)) {
		pvs.misaln[i] <- length(which(abs(T0s.misaln[[i]]) >= abs(Ts.misaln[i])))/n.perm
	}
	Q.misaln <- min(pvs.misaln)
	
	Q0.misaln <- rep(NA, n.perm)
	U.a.func <- function(y) return(t(r.s[[l]])%*%y%*%r.s[[l]])
	for (l in 1:n.perm) {
		Q0.misaln.s.n <- list()
		for (m in 1:length(pow)) {
			Q0.misaln.s.n[[m]] <- T0s.misaln[[m]][-l]
		}
		a.Ts.misaln <-  unlist(lapply(R.misaln,U.a.func))
		Q0.misaln[l] <- min(unlist(mapply(function(x,y)length(which(abs(x) >= abs(y)))/(n.perm-1),Q0.misaln.s.n,a.Ts.misaln)))
	}
	p.omisaln <- length(which(Q0.misaln <= Q.misaln))/n.perm

	names(Ts.misaln) <- paste("MiSALN(", pow, ")", sep = "")
	names(Q.misaln) <- "OMiSALN"
	names(pvs.misaln) <- paste("MiSALN(", pow, ")", sep = "")
	names(p.omisaln) <- "OMiSALN"
	return(list(pvs.misaln=pvs.misaln, p.omisaln=p.omisaln))
}
