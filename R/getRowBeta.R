`getRowBeta` <-
function(u, u_jk, j, k, b1, b2, X_i)
{
  a31 <- -1/(1-u[j]-u[k] + u_jk) - 1/(u[j] - u_jk) 
  a32 <- -1/(1-u[j]-u[k] + u_jk) - 1/(u[k] - u_jk)
  a33 <- 1/u_jk + 1/(1-u[j]-u[k] + u_jk) +  1/(u[j] - u_jk) +  1/(u[k] - u_jk)
  if ( any( is.infinite( c(a31, a32, a33) ) ) )
  {
    print("Division by 0 has occurred in the function 'getRowBeta()'")
    return(list(row=NULL, success=FALSE))
    # It is a good idea to have a mechanism that does some kind of
    # checking prior to calling this function.  The reason we want
    # this is because in a simulation setting, we want the program to
    # continue without throwing a tantrum and quitting. 
  }
  row <- ( -(a31/a33 + b1) * u[j] * ( 1-u[j] ) * X_i[j,] ) - ( (a32/a33 + b2) * u[k] * ( 1-u[k] ) * X_i[k,] )
  return( list(row=row, success=TRUE) )
}

