`mardia` <-
function(uj, uk, g)
{
  psi <- exp(g)
  if (psi == 1) { return(uj*uk) }
  a  <- uj + uk
  b  <- psi - 1
  cc <- 1 + a*b
  d  <- cc - sqrt(cc^2 - 4*psi*b*uj*uk)
  return( (0.5*d)/b )
}

