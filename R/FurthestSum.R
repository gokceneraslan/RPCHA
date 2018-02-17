
mysqrt <- function(x) {
  if(any(x<0)) {
    att <- attributes(x)
    x <- as.complex(x)
    attributes(x) <- att
  }
  sqrt(x)
}

mymax <- function(x) {
 if(is.complex(x))x[which.max(abs(x))] else max(x)
}

FurthestSum <- function(K, noc, i, exclude=c()) {

  I <- nrow(K)
  J <- ncol(K)
  index <- seq_len(J)
  index[exclude] <- 0
  index[i] <- 0

  ind_t <- i
  sum_dist <- rep(0, J) #row vec.

  if (J>noc*I) {

    Kt <- K
    Kt2 <- colSums(Kt^2)
    for (k in seq_len(noc+10)) {
        if (k>noc-1) { # Remove initial seed
           Kq <- t(Kt[,i[1]]) %*% Kt
           sum_dist <- sum_dist - mysqrt(Kt2 - 2*Kq + Kt2[i[1]])
           index[i[1]] <- i[1]
           i <- i[-1]
        }
        t <- which(index!=0)
        Kq <- t(Kt[, ind_t]) %*% Kt
        sum_dist <- sum_dist + mysqrt(Kt2 - 2*Kq + Kt2[ind_t])
        val <- mymax(sum_dist[t])
        ind <- which(sum_dist[t]==val)[1]
        ind_t <- t[ind[1]]
        i <- c(i, t[ind[1]])
        index[t[ind[1]]] <- 0
    }

  } else {
     if (I!=J || sum(K-t(K)) !=0) { # Generate kernel if K not a kernel matrix
        Kt=K
        K=t(Kt) %*% Kt
        #K=sqrt(repmat(diag(K)',J,1)-2*K+repmat(diag(K),1,J));
        K = mysqrt(matrix(rep(diag(K), J), nrow=J, byrow = T) - 2*K +
                   matrix(rep(diag(K), J), ncol=J))
     }
    Kt2 = diag(K)
    for (k in seq_len(noc+10)) {
      if (k>noc-1) {
        sum_dist <- sum_dist - mysqrt(Kt2 - 2*K[i[1],] + Kt2[i[1]])
        index[i[1]] <- i[1]
        i <- i[-1]
      }
      t <- which(index!=0)
      sum_dist <- sum_dist + mysqrt(Kt2 - 2*K[ind_t,] + Kt2[ind_t])
      val <- mymax(sum_dist[t])
      ind <- which(sum_dist[t]==val)[1]
      ind_t <- t[ind[1]]
      i <- c(i, ind_t)
      index[t[ind[1]]] <- 0
    }
  }
  return(i)
}

test.furthest <- function() {
  X = as.matrix(read.csv('X.csv', header=F, stringsAsFactors = F))

  noc=3 # Number of archetypes

  U=1:ncol(X) #Entries in X used that is modelled by the AA model
  I=1:ncol(X) #Entries in X used to define archetypes
  # if two expensive to useall entries for I find N relevant observations by
  # the following procedure:
  N=50
  I=FurthestSum(X, N, ceiling(0.2*ncol(X)))

  I.correct <- read.csv('FurthestSum.csv', header=F, stringsAsFactors = F)
  all(I == I.correct)

}