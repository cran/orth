`R.inv.A` <-
function(ni, lambda, A)
{
  if ( is.vector(A) ) { A <- matrix( A, nrow = length(A) ) }
  ab <- initConst(ni, lambda)
  c1 <- ab$c1
  c2 <- ab$c2
  return( c1*A + c2*rep(1, ch2(ni)) %*% matrix(apply(A, 2, sum), nrow=1) )
}

