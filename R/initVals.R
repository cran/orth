`initVals` <-
function(y, X, n, w, q, iA, iB)
{
  if (missing(iA)) { iA <- rep(0.01, q) }
  if (missing(iB)) { iB <- initBeta(y,X,n,w) }
  return( list( alpha = iA, beta = iB ) )
}

