`u21` <-
function(v)
{
  v[lower.tri(v)] <- 0
  return(v + t(v) - diag(diag(v)))
}

