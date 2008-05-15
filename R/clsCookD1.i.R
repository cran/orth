`clsCookD1.i` <-
function(df.bet, invRobust, invRobust.b, Ustar, p1, p2)
{
  cookD     <- t(df.bet) %*% invRobust[1:p1, 1:p1] %*% df.bet / (p1 + p2)
  cookD.bet <- t(df.bet) %*% invRobust.b %*% df.bet  / p1
  cookD.alp <- 0
  
  cookD.robust <- c(cookD, cookD.bet, cookD.alp)
  
  naive.cookD     <- t(df.bet) %*% Ustar[1:p1, 1:p1] %*% df.bet  /(p1+p2)
  naive.cookD.bet <- t(df.bet) %*% Ustar[1:p1, 1:p1] %*% df.bet  /(p1)
  naive.cookD.alp <- 0
  
  cookD.naive     <- c(naive.cookD, naive.cookD.bet, naive.cookD.alp)
   
  return ( list( cookD.robust = cookD.robust, cookD.naive = cookD.naive ) )
}

