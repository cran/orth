`invBig` <-
function(aic, aim, M, cc, first, last)
{
  aic <- aic
  aim <- aim
  for (i in seq(first, last) )
  {
    b <- aim[,i]
    btm <- t(b) %*% M
    gam <- 1 - btm[,i]
    bg <- b/gam
    aic <- aic + bg %*% ( t(b)%*%cc )
    if (i < last) { aim <- aim + bg %*% btm }
  }
  return( list(aic=aic, aim=aim) )
}

