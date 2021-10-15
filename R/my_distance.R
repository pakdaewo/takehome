my_distance <- function(X, group, type = "euclidean") {

  group_name <- unique(group)
  ngroup <- length(group_name)

  distance <- matrix(0, nrow = 5, ncol = 5)
  cov_X <- cov(X)
  V <- diag(cov_X)
  p <- ncol(X)

  for (i in 1:ngroup) {

    for (j in i:ngroup) {

      mui <- colMeans(X[group == group_name[i], ])
      muj <- colMeans(X[group == group_name[j], ])

      if (type == "euclidean") distance[i, j] <- sqrt(sum((mui - muj)^2)) # euclidean
      if (type == "penrose") distance[i, j] <- sqrt(sum((mui - muj)^2/(p*V))) #penrose
      if (type == "mahalanobis") distance[i, j] <- t(mui - muj) %*% solve(cov_X) %*% (mui - muj) # mahalanobis

    }

  }

  colnames(distance) <- rownames(distance) <- group_name

  return(t(distance))

}
