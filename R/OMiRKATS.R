OMiRKATS <-
function(obstime, delta, X, total.reads = NULL, tree, cov = NULL, g.unif.alpha=c(0.5), n.perm = 5000) {
	if (is.null(cov)) {
		cov <- as.matrix(rep(1,length(obstime)))
	}
	
	r <- coxph(Surv(obstime, delta) ~ ., data=as.data.frame(cov))$residuals
	r.s <- list()
	for (j in 1:n.perm) {
		r.s[[j]] <- r[shuffle(length(r))]
	}

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
	names(Ts.mirkats) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("Q.MiRKATS(", g.unif.alpha, ")", sep = ""))
	names(Q.mirkats) <- "OMiRKAT-S"
	names(pvs.mirkats) <- c("Bray-Curtis", "U.UniFrac", "W.UniFrac", paste("G.UniFrac(", g.unif.alpha, ")", sep = ""))
	names(p.omirkats) <- "OMiRKAT-S"
	return(list(pvs.mirkats=pvs.mirkats, p.omirkats=p.omirkats))
}
