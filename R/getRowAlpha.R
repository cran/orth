`getRowAlpha` <-
function(u_j, u_k, u_jk, Z)
{
  a33 <- 1/u_jk + 1/(1-u_j-u_k + u_jk) +  1/(u_j - u_jk) +  1/(u_k - u_jk)
  if ( is.infinite(a33) )
  {
    print("Division by zero has occured in the function 'getRowAlpha()' ");
    return( list(row=NULL, success=FALSE) )
  }
  else { return( list(row=Z/a33, success=TRUE) ) }
}

