`expandWeights` <-
function(n,w)
{
  owv <- rep(0,sum(n))
  fl <- beginEnd(n)
  for (i in 1:length(n))
  {
    owv[fl$first[i]:fl$last[i]] <- w[i]
  }
  return(owv)
}

