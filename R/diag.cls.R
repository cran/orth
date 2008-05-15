`diag.cls` <-
function(object)
{
  temp <- x.mat.orth(object)
  y <- temp$y
  x <- temp$X
  z <- z.mat.orth(object)
  
  x.names <- colnames(x)
  z.names <- colnames(z)
  
  beta   <- object@betas
  alpha  <- object@alphas
  lambda <- object@lambda
  

  #cluster size vector
  n <- cls.size(temp$id)
  
  xindx <- beginEnd(n)
  zindx <- beginEnd(ch2(n))

  first.x <- xindx$first
  last.x  <- xindx$last

  first.z <- zindx$first
  last.z  <- zindx$last

  m <- ch2(n)
  K <- length(n)

  p1 <- length(as.vector(object@betas))
  p2 <- length(as.vector(object@alphas))

  inv.robust   <- solve(object@variance$robust)
  inv.robust.b <- solve(object@variance$robust[1:p1, 1:p1])
  inv.robust.a <- solve(object@variance$robust[(p1+1):(p1+p2), (p1+1):(p1+p2)])

  # These objects store cluster level dffits statistics
  dfBeta.cls  <- matrix( rep(0, K*p1),ncol=p1 )
  dfAlpha.cls <- matrix( rep(0, K*p2),ncol=p2 )

  # These objects store cluster level leverage statistics
  H.bet <- rep(0, K)
  H.alp <- rep(0, K)

  # These objects store cluster level Cook's distance
  cookD.robust <- matrix(rep(0,K*3), ncol=3)
  cookD.naive  <- matrix(rep(0,K*3), ncol=3)


  sqee <- sqe( naiveA = object@variance$naiveA, naiveB = object@variance$naiveB )

  for (i in 1:K)  ## BEGIN: Main Loop ##
  {
    x.i  <- x[first.x[i]:last.x[i],]
    y.i  <- y[first.x[i]:last.x[i]]

    if (n[i]==1)
    {
      z.i <- matrix(rep(0,p2),nrow=1)
      x.i <- matrix(x.i, nrow=1)
    }
    else
    {
      # This pair of if-else statement is to prevent R from converting 
      # matrices back into vectors.  This happens when when only have
      # one column for the Z matrix.  When ever you extract a single row
      # or a single column from a matrix, R converts it to a vector.  This
      # causes no end of head-ache if not caught.
      if (p2 == 1) { z.i <- matrix( z[ first.z[i]:last.z[i], ], ncol=1 ) }
      else { z.i <- z[ first.z[i]:last.z[i], ] }
      
      if (n[i]==2) { z.i <- matrix(z.i,nrow=1) }
      ## This last "if" is needed because when you extract     ##
      ## a single row from a matrix, R converts it to a column ##
      ## vector resulting in problems with matrix operations   ##
      ## later on.  When n[i]==2, we only extract a single row ##
      ## from the association matrix Z.                        ##
    }

    mu.i <- 1/(  1 + exp(-x.i %*% beta)  )

    gamma.i <- z.i %*% alpha

    if (n[i] == 1)
    {
      H.bet[i] <-  clsLvg1.i(x.i, y.i, mu.i, object@variance$naiveB)
      df.beta  <- cls.df1.i(x.i, y.i, mu.i, H.bet[i], object@variance$naiveB)  
      dfBeta.cls[i, ]  <-   t( df.beta  )
      cD <- clsCookD1.i(df.beta, inv.robust, inv.robust.b, object@score$Ustar, p1, p2)
      cookD.robust[i, ] <- t(cD$cookD.robust)
      cookD.naive[i, ] <- t(cD$cookD.naive)
    }
    else
    {
      M <- abc(n[i], mu.i, gamma.i, x.i, y.i)
      vabr <- vee.dAdB.R(n[i], mu.i, gamma.i, p1, p2, x.i, y.i, z.i, i)

      sqVEE <- sqrt(vabr$VEE)
      const1 <- 1/(1-lambda)
      const2 <- -lambda / ( (1-lambda) * ( 1 + (m[i]-1)*lambda ) )

      df.fits <- cls.df2.i(vabr$DA, sqVEE, sqee$sqe1, sqee$sqe2, vabr$R, M$A, M$invB, M$CC,
                           const1, const2, m[i], p1, p2, object@variance)

      dfBeta.cls[i, ]  <- t(df.fits$dfBeta.cls)
      dfAlpha.cls[i, ] <- t(df.fits$dfAlpha.cls)

      H <- clsLvg2.i(object@variance$naiveA, object@variance$naiveB, M$CC, M$invB, vabr$DA, n[i], m[i],
                     const1, const2, sqVEE)
      H.bet[i] <- H$Hbet.i
      H.alp[i] <- H$Halp.i

      cD <- clsCookD2.i(df.fits$dfAlpha.cls, df.fits$dfBeta.cls, inv.robust, inv.robust.b, inv.robust.a, object@score$Ustar, p1, p2)
      cookD.robust[i, ] <- t(cD$cookD.robust)
      cookD.naive[i, ]  <- t(cD$cookD.naive)

    }

  } ## END:  Main Loop ##
  
  dfBeta.cls   <- cbind(unique(temp$id), n, dfBeta.cls)
  rownames(dfBeta.cls) <- seq(1:length(n))
  colnames(dfBeta.cls) <- c("cluster", "cluster.size", x.names)
  dfAlpha.cls  <- cbind(unique(temp$id), n, dfAlpha.cls)
  rownames(dfAlpha.cls) <- seq(1:length(n))
  colnames(dfAlpha.cls) <- c("cluster", "cluster.size", z.names)
  H.bet        <- cbind(unique(temp$id), n, H.bet)
  rownames(H.bet) <- seq(1:length(n))
  colnames(H.bet) <- c("cluster", "cluster.size", "leverage")
  H.alp        <- cbind(unique(temp$id), n, H.alp)
  rownames(H.alp) <- seq(1:length(n))
  colnames(H.alp) <- c("cluster", "cluster.size", "leverage")
  cookD.robust <- cbind(unique(temp$id), n, cookD.robust)
  rownames(cookD.robust) <- seq(1:length(n))
  colnames(cookD.robust) <- c("cluster", "cluster.size", "overall.cook", "cook.beta", "cook.alpha")
  cookD.naive  <- cbind(unique(temp$id), n, cookD.naive)
  rownames(cookD.naive) <- seq(1:length(n))
  colnames(cookD.naive) <- c("cluster", "cluster.size", "overall.cook", "cook.beta", "cook.alpha")
  
  cls <- new("cluster.diag.orth", dfbet=dfBeta.cls, dfalp=dfAlpha.cls, h.bet=H.bet, h.alp=H.alp,
             cook.robust = cookD.robust, cook.naive=cookD.naive)
  
  return( cls )
}

