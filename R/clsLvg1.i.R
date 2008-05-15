`clsLvg1.i` <-
function(x.i, y.i, mu.i, naiveB)
{
  v1       <- mu.i*(1-mu.i)    # 1 by 1
  v2       <- x.i * as.vector(v1)         # 1 by p1
  H.bet.i  <- v2 %*% naiveB %*% t(x.i)
  return(H.bet.i)
}

