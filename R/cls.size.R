`cls.size` <-
function(id)
{
  return( as.vector( table( as.factor(id) ) ) )
}

