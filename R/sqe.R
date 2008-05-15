`sqe` <-
function(naiveA, naiveB)
{
  if (!missing(naiveA) & !missing(naiveB))
  {
    ee1 <- eigen(naiveB, symmetric=T)
    ee2 <- eigen(naiveA, symmetric=T)
    rooteval1 <- sqrt(ee1$values)
    rooteval2 <- sqrt(ee2$values)
    
    if (length(rooteval1) == 1) { sqe1 <- ee1$vectors %*% rooteval1 }
    else { sqe1 <- ee1$vectors %*% diag(rooteval1) }

    if ( length(rooteval2) == 1 ) { sqe2 <- ee2$vectors %*% rooteval2 }
    else { sqe2 <- ee2$vectors %*% diag(rooteval2) }
    
    return( list( sqe1 = sqe1, sqe2 = sqe2 ) )
  }
  else if (missing(naiveB))
  {

    ee2 <- eigen(naiveA, symmetric=T)
    rooteval2 <- sqrt(ee2$values)
    
    if ( length(rooteval2) == 1 ) { sqe2 <- ee2$vectors %*% rooteval2 }
    else { sqe2 <- ee2$vectors %*% diag(rooteval2) }
    
    return( list( sqe1 = NULL, sqe2 = sqe2 ) ) 
  }
}

