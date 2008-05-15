`beginEnd` <-
function(n)
{
  last <- cumsum(n)
  first <- last - n+1
  return( list(first=first, last=last) )
}

