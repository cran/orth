`clsCookD2.i` <-
function(df.alp, df.bet, invRobust, invRobust.b, invRobust.a, Ustar, p1, p2)
{
  dfboth     <- c(df.bet, df.alp)

  # These are Cook's distance based on the robust covariance.
  # cookD is overall, cookD.bet is for beta, cookD.alp is for alpha
  cookD      <- ( t(dfboth) %*% invRobust %*% dfboth ) / (p1 + p2)
  cookD.bet  <- ( t(df.bet) %*% invRobust.b %*% df.bet ) / p1
  cookD.alp  <- ( t(df.alp) %*% invRobust.a %*% df.alp ) / p2

  # These are Cook's distance based on the model-based covariance.
  # naive.cookD is overall, naive.cookD.bet is for beta, naive.cookD.alp is for alpha
  naive.cookD <- ( t(dfboth) %*% Ustar %*% dfboth ) / (p1 + p2)
  naive.cookD.bet <- ( t(df.bet) %*% Ustar[1:p1, 1:p1] %*% df.bet ) / p1
  naive.cookD.alp <- ( t(df.alp) %*% Ustar[(p1+1):(p1+p2), (p1+1):(p1+p2)] %*% df.alp ) / p2

  cookD.robust <- c(cookD, cookD.bet, cookD.alp)
  cookD.naive  <- c(naive.cookD, naive.cookD.bet, naive.cookD.alp)
  return (list( cookD.robust = cookD.robust, cookD.naive=cookD.naive ) )
}

