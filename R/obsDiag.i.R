`obsDiag.i` <-
function(ni, abc.obj, orthVar, invRobust, invRobust.b, invRobust.a, p1, p2, lambda, vabr.obj, Ustar, sqe2)
{
  dfbetobs.i <- matrix( rep(0, ni*p1), ncol=p1)
  dfalpobs.i <- matrix( rep(0, ni*p2), ncol=p2)

  cookDobs.robust <- matrix( rep(0, ni*3), ncol=3)
  cookDobs.naive  <- matrix( rep(0, ni*3), ncol=3)

  for (ii in 1:ni)
  {
    vii <- abc.obj$B[-ii, ii]      # (ni-1) by 1    
    wii <- abc.obj$invB[-ii, ii]   # (ni-1) by 1
    
    ## if cluster size is 2, we have to force vii and wii to become a 1 by 1 matrix
    if (ni == 2)
    {
      vii <- matrix(vii,nrow=1)
      wii <- matrix(wii, nrow=1)
    }
    
    Wti <- abc.obj$invB[-ii, -ii] - (wii %*% t(wii)) / as.vector(abc.obj$invB[ii,ii]) # (ni-1) by (ni-1) #
    WV  <- Wti %*% vii  # This is of dimension (ni-1) by 1

    ## The matrix CC is of dimension (ni by p1). Therefore, CC[ii,] is (1 by p1).
    ## Similarly, CC[-ii,] is of dimension (ni-1) by p1. 

#    print("abc.obj$CC[-ii] is"); print(abc.obj$CC[-ii,])
#    print("WV is"); print(WV)
    
    if (ni > 2) { Ctilde <- abc.obj$CC[ii,] - t(abc.obj$CC[-ii,]) %*% WV }
    else { Ctilde <- abc.obj$CC[ii,] - abc.obj$CC[-ii,] %*% WV  }
    
    Htilde <- t(Ctilde) %*% orthVar$naiveB %*% Ctilde        # 1 by 1
    
    Stilde <- abc.obj$B[ii,ii] - t(vii)%*% WV                    # 1 by 1
    Atilde <- abc.obj$A[ii] - t(abc.obj$A[-ii]) %*% WV           # 1 by 1
    
    hdobs   <- (Ctilde %*% Atilde) / as.vector(Stilde - Htilde)  # p1 by 1
    dfbetii <- orthVar$naiveB %*% hdobs                      # p1 by 1
    
    dfbetobs.i[ii, ] <- t(dfbetii) ## DFFITS statistic for beta for cluster i
    
    ## Cook's D based on robust variance
    cook.bet.ii <- ( t(dfbetii)%*%invRobust.b %*% dfbetii ) / p1

    oc <- initIndex1(ni)
    comp <- oc$comp
    obs  <- oc$obs
    
    count <- 0
    ncomp <- 1
    nii   <- 1
    oc2 <- initIndex2 (ni, count, ncomp, nii, comp, obs, ii) # returns obs and comp again.
    
    cc <- initConst(ni, lambda)   # initializes constant
    c1 <- cc$c1                   
    c2 <- cc$c2
    
#    print("oc2 is "); print(oc2)
    
    vabr2 <- vee.dAdB.R2(oc2, vabr.obj, ni)  # Returns  modified DA and modified sqrt VEE
    DAt   <- vabr2$DA
    sqVEE <- vabr2$sqVEE
    Rt    <- vabr2$R
    
    mm2t <- DAt %*% sqe2
    colmm2t <- mm2t/sqVEE
    
    ai2Rt <- ( c1*(Rt/sqVEE) + c2*sum(Rt/sqVEE) )/ sqVEE
    ai2m2t <- ( c1*colmm2t + c2*( matrix(rep(1,nrow(colmm2t)),ncol=1) %*% matrix(apply(colmm2t,2,sum),nrow=1)  ) ) /sqVEE
    
    ai2Rt2 <- invBig(ai2Rt, ai2m2t, mm2t, Rt, 1, p2)
    
#    print("DAt is"); print(DAt)
#    print("ai2Rt2$aic"); print(ai2Rt2$aic)
    
    if (ni == 2) {U.alp <- DAt %*% ai2Rt2$aic}
    else if (ni > 1) { U.alp <- t(DAt) %*% ai2Rt2$aic }
    
    dfalpii <- orthVar$naiveA %*% U.alp
    
    dfalpobs.i[ii, ] <- t(dfalpii)  ## DFFITS statistis for alpha for cluster i
    
    dfboth <- c(dfbetii, dfalpii)
    cookDobs.ii <- t(dfboth) %*% invRobust %*% dfboth / (p1 + p2)
    cook.alp.ii <- t(dfalpii) %*% invRobust.a %*% dfalpii / p2
    
    naive.cookD.ii <- t(dfboth) %*% Ustar %*% dfboth / (p1 + p2)
    naive.cookD.bet.ii <- t(dfbetii) %*% Ustar[1:p1, 1:p1] %*% dfbetii / p1
    naive.cookD.alp.ii <- t(dfalpii) %*% Ustar[(p1+1):(p1+p2),(p1+1):(p1+p2)] %*% dfalpii / p2
    
    cookDobs.robust[ii,] <- t( c(cookDobs.ii, cook.bet.ii, cook.alp.ii) )
    cookDobs.naive[ii,]  <- t( c(naive.cookD.ii, naive.cookD.bet.ii, naive.cookD.alp.ii) )
  }
  return( list(df.bet.obs.i = dfbetobs.i, df.alp.obs.i = dfalpobs.i, cookD.robust.i = cookDobs.robust,
               cookD.naive.i = cookDobs.naive) )
  
}

