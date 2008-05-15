`diag.obs` <-
function(object)
{
    temp <- x.mat.orth(object)
    y <- temp$y
    x <- temp$X
    z <- z.mat.orth(object)
    
    x.names <- colnames(x)
    z.names <- colnames(z)
    
    beta  <- object@betas
    alpha <- object@alphas

    n <- cls.size(temp$id)

    p1 <- length(beta)
    p2 <- length(alpha)
    
    inv.robust   <- solve(object@variance$robust)
    inv.robust.b <- solve(object@variance$robust[1:p1, 1:p1])
    inv.robust.a <- solve(object@variance$robust[(p1+1):(p1+p2), (p1+1):(p1+p2)])
    
    x.indx  <- beginEnd(n)
    z.indx  <- beginEnd(ch2(n))
    
    first.x <- x.indx$first
    last.x  <- x.indx$last
    first.z <- z.indx$first
    last.z  <- z.indx$last
    
    K <- length(n)  # Number of clusters
    nobs <- sum(n)  # Number of observations
    
    dfBet.obs <- matrix( rep(0,nobs*p1), ncol=p1 )
    dfAlp.obs <- matrix( rep(0,nobs*p2), ncol=p2 )
    cookD.robust.obs <- matrix( rep(0,nobs*3), ncol=3 )
    cookD.naive.obs  <- matrix( rep(0,nobs*3), ncol=3 )
    
    sqee <- sqe(naiveA = object@variance$naiveA)
    sqe2 <- sqee$sqe2
    
    for (i in 1:K)
    {
        x.i <- x[ first.x[i]:last.x[i], ]
        y.i <- y[ first.x[i]:last.x[i] ]
        
        if ( n[i] == 1 )
        {
            z.i <- matrix( rep(0,p2), nrow=1 )
            x.i <- matrix( x.i, nrow=1 )
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

            if (n[i] == 2) { z.i <- matrix(z.i, nrow=1) }
        }
        
        mu.i    <- 1 / ( 1 + exp(-x.i %*% beta) )
        gamma.i <- z.i %*% alpha
        
        if (n[i] == 1)
        {
            B   <- mu.i * (1 - mu.i)
            CC  <- x.i * matrix(B, ncol=p1)
            Hi1 <- as.vector( CC %*% object@variance$naiveB %*% t(x.i) )
            U.i <- t(x.i) %*% (y.i - mu.i) / (1 - Hi1)
            dfBet.obs[first.x[i], ] <- t( object@variance$naiveB %*% U.i )
            
            # Computing Cook's distance
            H.bet.i <- clsLvg1.i(x.i, y.i, mu.i, object@variance$naiveB)
            df.beta <- cls.df1.i(x.i, y.i, mu.i, H.bet.i, object@variance$naiveB)
            cD      <- clsCookD1.i(df.beta, inv.robust, inv.robust.b, object@score$Ustar,p1,p2)
            cookD.robust.obs[first.x[i], ] <- t(cD$cookD.robust)
            cookD.naive.obs[first.x[i],  ] <- t(cD$cookD.naive )
        }# End of n[i] == 1 IF statement
        else
        {
            vabr <- vee.dAdB.R( n[i], mu.i, gamma.i, p1, p2, x.i, y.i, z.i, i )
            abcc <- abc( n[i], mu.i, gamma.i, x.i, y.i )
            diag.i <- obsDiag.i( n[i], abcc, object@variance, inv.robust, inv.robust.b, inv.robust.a, p1, p2, object@score$lambda,
                                 vabr, object@score$Ustar, sqe2 )
            dfBet.obs[first.x[i]:last.x[i], ]        <- diag.i$df.bet.obs.i
            dfAlp.obs[first.x[i]:last.x[i], ]        <- diag.i$df.alp.obs.i
            cookD.robust.obs[first.x[i]:last.x[i], ] <- diag.i$cookD.robust.i
            cookD.naive.obs[first.x[i]:last.x[i], ]  <- diag.i$cookD.naive.i
        }
    }
    
    long.n <- rep(n, n)
    
    cook.names <- c("cluster", "cluster.size", "overall.cook", "cook.beta", "cook.alpha")
    cookD.robust.obs <- cbind(temp$id, long.n, cookD.robust.obs)
    colnames(cookD.robust.obs) <- cook.names
    rownames(cookD.robust.obs) <- seq(1, length(temp$id))
    
    cookD.naive.obs <- cbind(temp$id, long.n, cookD.naive.obs)
    colnames(cookD.naive.obs)  <- cook.names
    rownames(cookD.naive.obs)  <- seq(1, length(temp$id))
    
    dfBet.obs <- cbind(temp$id, long.n,  dfBet.obs)
    colnames(dfBet.obs) <- c("cluster", "cluster.size", x.names)
    rownames(dfBet.obs) <- seq(1, length(temp$id))
    
    dfAlp.obs <- cbind(temp$id, long.n,  dfAlp.obs)
    colnames(dfAlp.obs) <- c("cluster", "cluster.size",z.names)
    rownames(dfAlp.obs) <- seq(1, length(temp$id))

    obs <- new("obsvnl.diag.orth", dfbet=dfBet.obs, dfalp=dfAlp.obs, cook.robust=cookD.robust.obs,
                                   cook.naive = cookD.naive.obs)

    return(obs)
}

