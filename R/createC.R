`createC` <-
function(x.i, u.i)
{
  return( x.i * matrix( rep(u.i,ncol(x.i)),ncol=ncol(x.i) ) * (1 - matrix( rep(u.i,ncol(x.i)),ncol=ncol(x.i) )) )
}

