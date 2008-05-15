`lambda.est` <-
function(N, nmr, dnm, niter)    
{
  lambda <- nmr/dnm
  if (lambda <= -1/( ch2(max(N))-1 ))
  {
    print ("LAMBDA OUT OF RANGE: NEGATIVE")
    print ("LAMBDA RESET TO 0.8 TIMES LOWER LIMIT")
    print ( paste("ITERATION", niter) )
    lambda <- 0.8 * ( -1/(ch2(max(N))-1) )
  }
  else if (lambda >= 1) 
  {
    print ("LAMBDA OUT OF RANGE: POSITIVE .")
    print ("LAMBDA RESET TO 0.8 .")
    print ( paste("ITERATION", niter) )
    lambda <- 0.8 
  }
  return(lambda) 
}

