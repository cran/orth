`abc` <-
function(ni, mu.i, gamma.i, x.i, y.i)
{
  CC <- createC(x.i, mu.i)
  B <- createB(mu.i, ni, gamma.i)
  A <- createA(mu.i, y.i)
  return( list( A=A, B=B, CC=CC, invB = solve(B) ) )
}

