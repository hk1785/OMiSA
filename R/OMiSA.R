OMiSA <-
function(obstime, delta, X, total.reads = NULL, tree, cov = NULL, pow = c(1/4, 1/3, 1/2, 1), g.unif.alpha=c(0.50), n.perm = 5000) {
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

	if (is.null(total.reads)) {
		unifs <- GUniFrac(X, tree, alpha = c(g.unif.alpha, 1))$unifracs
		bray.curtis <- as.matrix(bcdist(X))
		u.unif <- unifs[,,"d_UW"]
		w.unif <- unifs[,,"d_1"]
		g.unif <- list()
	} else {
		unifs <- GUniFrac2(X, tree, alpha = c(g.unif.alpha, 1), total.reads=total.reads)$unifracs
		bray.curtis <- as.matrix(bcdist(X))
		u.unif <- unifs[,,"d_UW"]
		w.unif <- unifs[,,"d_1"]
		g.unif <- list()
	}
	for (k in 1:length(g.unif.alpha)) {
		g.unif[[k]] <- unifs[,,paste("d_",g.unif.alpha[k],sep="")]
	}
	for (j in 1:length(obstime)) {
		ind <- is.na(u.unif[j,])
		if (sum(ind) != 0)
		u.unif[,ind] <- 0
		
		ind <- is.na(w.unif[j,])
		if (sum(ind) != 0)
		w.unif[,ind] <- 0		
	}
	for (k in 1:length(g.unif.alpha)) {
		for (j in 1:length(obstime)) {
			g.unif.ind <- g.unif[[k]]
			ind <- is.na(g.unif.ind[j,])
			if (sum(ind) != 0)
			g.unif.ind[,ind] <- 0	
			g.unif[[k]] <- g.unif.ind
		}
	}
	bray.curtis.kern <- D2K(bray.curtis)
	u.unif.kern <- D2K(u.unif)
	w.unif.kern <- D2K(w.unif)
	g.unif.kern <- list()
	for (k in 1:length(g.unif.alpha)) {
		g.unif.kern[[k]] <- D2K(g.unif[[k]])
	}
	list.kernels <- c(list(bray.curtis.kern = bray.curtis.kern, u.unif.kern = u.unif.kern, w.unif.kern = w.unif.kern), g.unif.kern)
	Ts.mirkats <- rep(NA, length(list.kernels))
	for (j in 1:length(list.kernels)) {
		Ts.mirkats[j] <- t(r)%*%list.kernels[[j]]%*%r
	}

	T0s.mirkats <- list()
	for (j in 1:length(list.kernels)) {
		T0s.mirkats.inv <- rep(NA, n.perm)
		for (k in 1:n.perm) {
			T0s.mirkats.inv[k] <- t(r.s[[k]])%*%list.kernels[[j]]%*%r.s[[k]]
		}
		T0s.mirkats[[j]] <- T0s.mirkats.inv
	}
	pvs.mirkats <- rep(NA, length(list.kernels))
	for (j in 1:length(list.kernels)) {
		pvs.mirkats[j] <- length(which(abs(T0s.mirkats[[j]]) >= abs(Ts.mirkats[[j]])))/n.perm
	}
	Q.mirkats <- min(pvs.mirkats)
	Q0.mirkats <- rep(NA, n.perm)
	for (l in 1:n.perm) {
		Q0.mirkats.s.n <- list()
		for (m in 1:length(list.kernels)) {
			Q0.mirkats.s.n[[m]] <- T0s.mirkats[[m]][-l]
		}
		a.Qs.mirkats <- unlist(lapply(list.kernels,function(x) return(t(r.s[[l]])%*%x%*%r.s[[l]])))
		a.pvs <- unlist(mapply(function(x,y)length(which(abs(x) >= abs(y)))/(n.perm-1),Q0.mirkats.s.n,a.Qs.mirkats))
		Q0.mirkats[l] <- min(a.pvs)
	}
	p.omirkats <- length(which(Q0.mirkats <= Q.mirkats))/n.perm
	Q.omisa <- min(Q.misaln, Q.mirkats)
	Q0.omisa <- apply(cbind(Q0.misaln, Q0.mirkats),1,min)
	p.omisa <- length(which(Q0.omisa <= Q.omisa))/n.perm
	names(Ts.misaln) <- paste("MiSALN(", pow, ")", sep = "")
	names(Ts.mirkats) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("Q.MiRKATS(", g.unif.alpha, ")", sep = ""))
	names(Q.misaln) <- "OMiSALN"
	names(Q.mirkats) <- "OMiRKAT-S"
	names(Q.omisa) <- "OMiSA"
	names(pvs.misaln) <- paste("MiSALN(", pow, ")", sep = "")
	names(pvs.mirkats) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("G.UniFrac(", g.unif.alpha, ")", sep = ""))
	names(p.omisaln) <- "OMiSALN"
	names(p.omirkats) <- "OMiRKAT-S"
	names(p.omisa) <- "OMiSA"
	return(list(pvs.misaln=pvs.misaln, pvs.mirkats=pvs.mirkats, p.omisaln=p.omisaln, p.omirkats=p.omirkats, p.omisa=p.omisa))
}
