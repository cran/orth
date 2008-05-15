`vee` <-
function(mui, mujk, j, k)
{
  num <- (mui[j] - mujk) * mujk * (mui[k] - mujk) * (-1 + mui[j] + mui[k] - mujk)
  den <- mujk^2 - 2*mui[j]*mui[k]*mujk + mui[j]*mui[k]*(mui[j] + mui[k] - 1)
  return(num/den)
}

