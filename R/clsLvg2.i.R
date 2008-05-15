`clsLvg2.i` <-
function(naiveA, naiveB, CC, invB, DA, ni, mi, c1, c2, rootVEE)
{
  H.bet   <- CC %*% naiveB %*% ( t(CC) %*% invB )
  H.bet.i <- sum( diag(H.bet) )

  colA     <- DA/rootVEE
  p.inv.Da <- ( c1*colA + c2*rep(1,mi) %*% matrix(apply(colA,2,sum), nrow=1) ) /rootVEE
  Hi2Diag  <- rep(0,mi)
  l <- 1
  for (j in seq(1,ni-1) )
  {
    for (k in seq(j+1, ni) )
    {
      Hi2Diag[l] <- DA[l, ] %*% naiveA %*% p.inv.Da[l,]
      l <- l + 1
    }
  }
  H.alp.i <- sum( Hi2Diag )
  return(list(Hbet.i = H.bet.i, Halp.i = H.alp.i) )
}

