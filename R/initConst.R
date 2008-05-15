`initConst` <-
function(ni, lambda)
{
  # When n[i] equals 2, the contribution of lambda is 0
  if (ni == 2) { c1 <- 1; c2 <- 0 }
  else
  {
    c1 <- 1 / (1-lambda)
    c2 <- -lambda / ( (1-lambda)*(1 + (ch2(ni) - 1)* lambda ) ) 
  }
  return(list(c1 = c1, c2 = c2))
}

