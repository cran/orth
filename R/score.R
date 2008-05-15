`score` <-
function(bet, alp, y, X, Z, n, w, p1, p2, estLambda=FALSE, lambda=0, flag=FALSE, niter)
{
  U      <- rep(0,p1+p2)
  UUtran <- matrix( rep(0,(p1+p2)^2),nrow=p1+p2 )
  Ustar  <- matrix( rep(0,(p1+p2)^2),nrow=p1+p2 )

  if (estLambda) 
  {
    sumres <- 0
    denom  <- 0
  }
    
  xIndex  <- beginEnd(n)
  first.x <- xIndex$first
  last.x  <- xIndex$last
  
  zIndex  <- beginEnd( ch2(n) )
  first.z <- zIndex$first
  last.z  <- zIndex$last
    
  for ( i in seq(1,length(n)) )  # \label{Loop A}
  {# Begin of loop A
    X_i <- X[first.x[i]:last.x[i],]  
    y_i <- y[first.x[i]:last.x[i]] 
     
    if (n[i]==1) 
    { 
      Z_i <- matrix(rep(0,p2),nrow=1) 
      X_i <- matrix(X_i, nrow=1)
    }
    else 
    { 
      # This pair of if-else statement is to prevent R from converting 
      # matrices back into vectors.  This happens when when only have
      # one column for the Z matrix.  When ever you extract a single row
      # or a single column from a matrix, R converts it to a vector.  This
      # causes no end of head-ache if not caught.  
      if (p2 == 1) { Z_i <- matrix( Z[first.z[i]:last.z[i],], ncol=1 ) }
      else { Z_i <- Z[ first.z[i]:last.z[i], ] } 
      
      if (n[i]==2) { Z_i <- matrix(Z_i,nrow=1) }
      ## This last "if" is needed because when you extract     ##
      ## a single row from a matrix, R converts it to a column ##
      ## vector resulting in problems with matrix operations   ##
      ## later on.  When n[i]==2, we only extract a single row ##
      ## from the association matrix Z.                        ##
    }
    
    if (n[i]==1)    # \label{if.else.A}
    {
      mu_i      <- 1 / (1 + exp(-X_i %*% bet) )    # dimension 1 x 1
      if (is.infinite(mu_i))
      {
        print ("Your BETAs are diverging in the score() function")
        return( list(U = NULL, UUtran=NULL, Ustar=NULL, lambda=NULL, success=FALSE) )
      }
      
      U_i       <- t(X_i) %*% (y_i - mu_i)   # dimension p1 x 1 
      
      UUtran_i  <- U_i %*% t(U_i)            # dimension p1 x p1
      Ustar_i   <- t(X_i) %*% ( mu_i*(1-mu_i) ) %*% X_i  # dimension p1 by p1
      # Note: the operation in Ustar works because n_i = 1. Otherwise, it
      # will bomb. So don't freak out.
      U[1:p1]            <- U[1:p1] + U_i * w[i]
      UUtran[1:p1, 1:p1] <- UUtran[1:p1, 1:p1] + UUtran_i * w[i]
      Ustar[1:p1, 1:p1]  <- Ustar[1:p1, 1:p1] + Ustar_i * w[i]
    }
    else if (n[i] > 1)
    {
      U_i     <- rep(0, p1 + p2)  
      Ustar_i <- matrix( rep(0,(p1+p2)^2), nrow=p1+p2 )
      mu_i    <- 1/( 1+ exp(-X_i %*% bet) ) # dimension n_i x 1
      if ( any( is.infinite(mu_i) ) )
      {
        print ("Your BETAs are diverging in the score() function")
        return( list(U = NULL, UUtran=NULL, Ustar=NULL, lambda=NULL, success=FALSE) )
      }
      gamma_i <- Z_i %*% alp;               # dimension ch2(n_i) x 1
      
      # Computes VEE, dA, dB, and R; see the vee.dA.dB.R() function for a 
      # description of what these are.
      vabr    <- vee.dAdB.R(n[i], mu_i, gamma_i, p1, p2, X_i, y_i, Z_i, i) 
      if (!vabr$success) 
      {
        print( paste("Failed to compute VEE, dA, dB, and R for cluster number ", i) )
        print( paste("This occured in iteration number", niter) )
        return( list(U = NULL, UUtran=NULL, Ustar=NULL, lambda=NULL, success=FALSE) ) 
      } 
      

      if (any(vabr$VEE < 0) )
      {
        print("The SCORE() routine failed because negative values were obtained for VEE.")
        print("The square root of VEE is needed in the routine which cannot be done if ")
        print("VEE < 0.") 
        return( list(U = NULL, UUtran=NULL, Ustar=NULL, lambda=NULL, success=FALSE) )
      }
      
      colA <- vabr$DA /sqrt(vabr$VEE) 
      colB <- vabr$DB /sqrt(vabr$VEE) 

      Cp <- t( createC(X_i, mu_i) )
      B  <- createB(mu_i, n[i], gamma_i)
      A  <- createA(mu_i, y_i)
      
      # Checks for positive definiteness of B
      if (flag)
      {
        if ( min(eigen(B, symmetric=T)$values) <= 0 ) 
        {
          print( paste("Cluster = ", i) )
          print("Var(Y) is not positive definite")
          print("Joint distribution does not exists")
          print("Your results may not be valid")
        } 
      } 
      
#      # If cluster size is 2, lambda does not appear in that cluster's contribution
#      cnstnt <- initConst(n[i], lambda)
#      const1 <- cnstnt$c1; const2 <- cnstnt$c2
      U_i[1:p1]           <- Cp %*% solve(B) %*% A
      U_i[(p1+1):(p1+p2)] <- t(vabr$DA) %*% ( R.inv.A( n[i], lambda, vabr$R/sqrt(vabr$VEE) ) / sqrt(vabr$VEE) )
      UUtran_i            <- U_i %*% t(U_i)
       
      Ustar_i[1:p1, 1:p1]   <- Cp %*% solve(B) %*% t(Cp)
      Ustar_i[(p1+1):(p1+p2), 1:p1] <- t(vabr$DA) %*% ( R.inv.A (n[i], lambda, colB) / sqrt(vabr$VEE) )
      Ustar_i[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- t(vabr$DA) %*% ( R.inv.A( n[i], lambda, colA ) / sqrt(vabr$VEE) )

      U      <- U + U_i*w[i]
      UUtran <- UUtran + UUtran_i * w[i]
      Ustar  <- Ustar + Ustar_i * w[i] 
      
      if ( (ch2(n[i]) > 1) & (estLambda) )
      {
        res.part   <- ( ( sum(vabr$R/sqrt(vabr$VEE)) )^2  - sum(vabr$R^2 / vabr$VEE) )/2
        denom.part <- ch2(ch2(n[i]))
        sumres     <- sumres + res.part * w[i] 
        denom      <- denom + denom.part * w[i]
      }                        
    } # End of if.else.A
  }# End of loop A
  if (estLambda) { lambda <- lambda.est(max(n), sumres, denom, niter) } 

  return( list(U = U, UUtran=UUtran, Ustar=Ustar, lambda=lambda, success=TRUE) )
}

