`initBeta` <-
function(y, X, n, w)
{
  X = as.matrix(X)
  owv <- expandWeights(n, w) # vector of observed weight vectors
  z <- y + y - 1
  B <- solve(  t(X)%*%(X*owv), t(X)%*%(z*owv)  )  
  for (i in 1:2)
  {
    u <- X %*% B
    u <- 1 / ( 1 + exp(-u) )
    v <- u * (1-u) * owv
    z <- t(X) %*% ( (y-u)*owv )
    d <- solve( t(X) %*% ( X*(matrix( rep(v,ncol(X)),ncol=ncol(X), byrow=F) ) ), z )
    B <- B + d
    # Note that the MATRIX function inside the "d <- solve(...)" line is so
    # R does not complain that X and v do not have the same dimension.
    # This because X is a matrix and v is a vector.  To perform elementwise
    # multiplication, they must have the same dimension.
  }
  return(as.vector(B))
}

