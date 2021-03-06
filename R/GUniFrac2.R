GUniFrac2 <-
function (otu.tab, tree, alpha = c(0, 0.5, 1), total.reads=total.reads) {
    if (!is.rooted(tree)) 
        stop("Rooted phylogenetic tree required!")
    otu.tab <- as.matrix(otu.tab)
    row.sum <- total.reads
    otu.tab <- otu.tab/row.sum
    n <- nrow(otu.tab)
    if (is.null(rownames(otu.tab))) {
        rownames(otu.tab) <- paste("comm", 1:n, sep = "_")
    }
    dimname3 <- c(paste("d", alpha, sep = "_"), "d_UW", "d_VAW")
    unifracs <- array(NA, c(n, n, length(alpha) + 2), dimnames = list(rownames(otu.tab), 
        rownames(otu.tab), dimname3))
    for (i in 1:(length(alpha) + 2)) {
        for (j in 1:n) {
            unifracs[j, j, i] <- 0
        }
    }
    if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
        stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
    }
    absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
    if (length(absent) != 0) {
        tree <- drop.tip(tree, absent)
        warning("The tree has more OTU than the OTU table!")
    }
    tip.label <- tree$tip.label
    otu.tab <- otu.tab[, tip.label]
    ntip <- length(tip.label)
    nbr <- nrow(tree$edge)
    edge <- tree$edge
    edge2 <- edge[, 2]
    br.len <- tree$edge.length
    cum <- matrix(0, nbr, n)
    for (i in 1:ntip) {
        tip.loc <- which(edge2 == i)
        cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
        node <- edge[tip.loc, 1]
        node.loc <- which(edge2 == node)
        while (length(node.loc)) {
            cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
            node <- edge[node.loc, 1]
            node.loc <- which(edge2 == node)
        }
    }
    cum.ct <- round(t(t(cum) * row.sum))
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            cum1 <- cum[, i]
            cum2 <- cum[, j]
            ind <- (cum1 + cum2) != 0
            cum1 <- cum1[ind]
            cum2 <- cum2[ind]
            br.len2 <- br.len[ind]
            mi <- cum.ct[ind, i] + cum.ct[ind, j]
            mt <- row.sum[i] + row.sum[j]
            diff <- abs(cum1 - cum2)/(cum1 + cum2)
            for (k in 1:length(alpha)) {
                w <- br.len2 * (cum1 + cum2)^alpha[k]
                unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * 
                  w)/sum(w)
            }
            ind2 <- (mt != mi)
            w <- br.len2 * (cum1 + cum2)/sqrt(mi * (mt - mi))
            unifracs[i, j, (k + 2)] <- unifracs[j, i, (k + 2)] <- sum(diff[ind2] * 
                w[ind2])/sum(w[ind2])
            cum1 <- (cum1 != 0)
            cum2 <- (cum2 != 0)
            unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <- sum(abs(cum1 - 
                cum2)/(cum1 + cum2) * br.len2)/sum(br.len2)
        }
    }
    return(list(unifracs = unifracs))
}
