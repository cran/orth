`createB` <-
function(u, n, lpsi)
{
  v = u*(1-u) #;print("v"); print(v);
  B <- diag( as.vector(u*(1-u)) ) # variance
#  print(B);
  l <- 1
  for ( j in seq(1,n-1) )
  {
    for ( k in seq(j+1 , n ) )
    {
      B[j,k] <- mardia( u[j], u[k], lpsi[l]) - u[j]*u[k] 
      l <- l+1
    }
  }
  return( u21(B) ) # The function u21(v) is in the file CLF.R;
                   # equivalent to SAS sqrsym(symsqr(B`)).
}

