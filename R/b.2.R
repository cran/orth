`b.2` <-
function(mui, mujk, j, k)
{
  numerator <- mujk*(1-mui[j])*(mui[j] - mujk)
  denominator <- mui[j]*(1-mui[j])*mui[k]*(1-mui[k]) - (mujk - mui[j]*mui[k])^2
  return(numerator/denominator)
}

