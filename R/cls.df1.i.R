`cls.df1.i` <-
function(x.i, y.i, mu.i, hbet.i, naiveB)
{
  U.i <- t(x.i) %*% (y.i - mu.i) / ( 1 - as.vector(hbet.i) )  # p1 by 1
  dfbetcls <- naiveB %*% U.i                                  # p1 by 1
  return(dfbetcls)
}

