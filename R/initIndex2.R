`initIndex2` <-
function(ni, count, ncomp, nii, comp, obs, ii)
{
  for (j in 1:(ni-1))
  {
    for (k in (j+1):ni)
    {
      count <- count + 1
      if ( (ii != j) & (ii != k) )
      {
        comp[ncomp] <- count
        ncomp <- ncomp + 1
      }
      else
      {
        obs[nii] <- count
        nii <- nii + 1
      }
    }
  }
  return( list(obs=obs, comp=comp) )
}

