`vee.dAdB.R` <-
function(ni, mu.i, gamma.i, p1, p2, x.i, y.i, z.i, i)
{  
  mi    <- ch2(ni)
  vee.i <- rep(0,mi)
  R.i   <- rep(0,mi)
  DA    <- matrix( rep(0,p2*mi), ncol=p2 )
  DB    <- matrix( rep(0,p1*mi), ncol=p1 )

  l <- 1
  for ( j in seq(1,ni-1) )   #label = loopj
  {
    for ( k in seq(j+1,ni) ) #label = loopk
    {
      mujk <- mardia( mu.i[j], mu.i[k], gamma.i[l] )
      b1   <- b.1( mu.i, mujk, j, k )
      b2   <- b.2( mu.i, mujk, j, k )
      
      vee.i[l] <- vee(mu.i, mujk, j, k)
      temp1 <- getRowBeta(mu.i, mujk, j, k, b1, b2, x.i)
      temp2 <- getRowAlpha(mu.i[j], mu.i[k], mujk, z.i[l,])
      if ( !temp1$success | !temp2$success )
      {
        print("Failure in the vee.dAdB.R() function ")
        return( list( VEE=NULL, DB=NULL, DA=NULL, R=NULL, success=F ) )
      }
      else
      {
        DB[l,] <- matrix(temp1$row, nrow=1)
        DA[l,] <- matrix(temp2$row, nrow=1)
        R.i[l]  <- y.i[j]*y.i[k] - mujk - b1*(y.i[j] - mu.i[j]) - b2*( y.i[k] - mu.i[k] )
      }
      l <- l+1
    } # End of loopk
  }   # End of loopj
  if (any( vee.i < 0 ))
  {
    print("Negative VEEs in the vee.dAdB.R() function")
    return( list( VEE=NULL, DB=NULL, DA=NULL, R=NULL, success=F ) )
  }
  else
  {
    return( list( VEE = vee.i, DB = DB, DA = DA, R = R.i, success=T ) )
  }

}

