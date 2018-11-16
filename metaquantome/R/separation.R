
# calculate the avg distance between n clusters
# divided by the total variation within clusters
sep_n <- function(clust){
    # clust is a list of dataframes, where the columns in each dataframe are the
    #
    nclust <- length(clust)
    means <- lapply(clust, colMeans)
    possible_dists <- combn(1:nclust, 2)
    n_possible_dists <- ncol(possible_dists)
    dists <- rep(0, n_possible_dists)
    for (i in 1:n_possible_dists){
        comb <- possible_dists[, i]
        dists[i] <- sum((means[[comb[1]]] - means[[comb[2]]])^2)
    }
    avg_dist <- mean(dists)
    within_variance <- sapply(1:nclust, function(i) mean((clust[[i]] - means[[i]])^2))
    avg_dist / sum(within_variance)
}
