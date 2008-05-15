`orthVar` <-
function(score, p1, p2)
{
  if (!score$success) 
  { 
    print("Covariance cannot be computed since SCORE() failed")
    return( list(robust=NULL, naiveA=NULL, naiveB=NULL) )
  }
  else
  {
    invUstar <- solve(score$Ustar)
    robust   <- invUstar %*% score$UUtran %*% t(invUstar)
    naiveB   <- solve( score$Ustar[1:p1,1:p1] )
    naiveA   <- solve( score$Ustar[(p1+1):(p1+p2),(p1+1):(p1+p2)] )
    return( list(robust=robust, naiveA=naiveA, naiveB=naiveB) )
  }
}

