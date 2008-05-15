`orth.fit` <-
function(y, X, Z, n, w, iAlp, iBet, lambda0=0,
                     maxiter=20, tol=1e-6, estLam=FALSE, monitor=FALSE)
{
  p1 <- ncol(X)
  p2 <- ncol(Z)
  delta <- rep(2*tol, p1+p2)

  if (is.null(iAlp) & is.null(iBet))        { theta.0 <- initVals(y, X, n, w, p2) }
  else if (is.null(iAlp)  & !is.null(iBet)) { theta.0 <- initVals(y, X, n, w, p2, iB=iBet) }
  else if (!is.null(iAlp) &  is.null(iBet)) { theta.0 <- initVals(y, X, n, w, p2, iA=iAlp) }
  else if (!is.null(iAlp) & !is.null(iBet)) { theta.0 <- initVals(y, X, n, w, p2, iA=iAlp, iB=iBet) }

  bet <- theta.0$beta
  alp <- theta.0$alpha

  for(niter in 1:maxiter)
  {
    if ( max(abs(delta)) <= tol  ) { break }
    score.obj <- score(bet, alp, y, X, Z, n, w, p1, p2, estLambda=estLam, lambda=lambda0, niter=niter)
    if (!score.obj$success) 
    {
      return( list( betas=c(NA), alphas=c(NA),  variance=list(robust=NULL, naiveA=NULL, naiveB=NULL),
                    num.iter=0, converge=FALSE, lambda=0, score=score.obj ) )
    }
    theta <- c(bet, alp)
    delta <- solve(score.obj$Ustar, score.obj$U)
    theta <- theta + delta
    bet <- theta[1:p1];
    alp <- theta[(p1+1):(p1+p2)]
    covar <- orthVar(score.obj, p1, p2) ## A list: (1)Robust Cov, (2) Naive cov(A), (3) Naive cov(B)

    if (monitor)
    {
    print(paste("iteration number", niter))
    print("beta:")
    print(bet)
    print("alpha:")
    print(alp)
    }
    # This line was added in September 17, 2007
    lambda0 <- score.obj$lambda

    converge <- ( max( abs(delta) ) <= tol )
  }
    niter <- niter - 1

    #========================================================================#
    ### BEGIN:  Testing the dffits.cls() routine                           ###
    ## First run the SCORE() routine one last time
#    score.obj <- score(bet, alp, y, X, Z, n, w, p1, p2, estLambda=F, lambda = score.obj$lambda, flag=T, niter=niter)
#    theta <- c(bet, alp)
#    delta <- solve(score.obj$Ustar, score.obj$U)
#    theta <- theta + delta
#    bet <- theta[1:p1];
#    alp <- theta[(p1+1):(p1+p2)]
#    orth.var.obj <- orthVar(score.obj, p1, p2)

    ### END: Testing the dffits.cls() routine                              ###
    #========================================================================#

                
    return( list(betas=bet, alphas=alp, variance = covar, num.iter=niter, converge=converge, 
                lambda=score.obj$lambda, score=score.obj ) )

}

