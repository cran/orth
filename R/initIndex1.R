`initIndex1` <-
function(ni)
{
  if (ni > 2)
  {
    obs  <- rep(0,ni-1)
    comp <- rep(0, ch2(ni-1))
  }
  else if (ni==2)
  {
    comp=NULL
    obs=1
  }
  return(list( obs=obs, comp=comp ) )
}

