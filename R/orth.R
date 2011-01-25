`orth` <-
function(formula, data, weights, formula.z, dataz, id, contrasts=NULL,
                 alp.start=NULL, bet.start=NULL, lambda0=0,
                 maxiter=20, tol=1e-6, estLam=FALSE, monitor=FALSE)
{

  if (missing(formula))   { stop("You need to provide a model formula for the mean.")        }
  if (missing(data))      { stop("You need to provide data for marginal mean effects.")      }
  if (missing(dataz))     { stop("You need to provide data for marginal association.")       }
  if (missing(formula.z)) { stop("You need to provide a model formula for the association.") }
  if (missing(id))        { stop("You need to provide the cluster id.") }

  call <- match.call(expand.dots=F)
  
  if ( !missing(weights) )
  {
    temp1 <- from.call.to.xyw2(call)

  }
  else
  {
    temp1 <- from.call.to.xyw1(call)
    warning("Cluster weights were assumed 1 since no \n", "privision was made by the user. \n\n")
  }


  x <- temp1$X  # Extract design matrix for mean parameters
  y <- temp1$y  # Extract the response vector
  w <- temp1$w  # Extract the weight vector
  z <- from.call.to.z(call)  # Extract design matrix for association parameters
  
  obs.x <- seq(1, nrow(x))
  obs.z <- seq(1, nrow(z))

  rownames(x) <- obs.x
  rownames(z) <- obs.z
  rownames(y) <- obs.x

  x.names <- colnames(x)
  z.names <- colnames(z)
  
  n <- as.vector( cls.size(temp1$id) )

  fit <- orth.fit(y, x, z, n, w=w, iAlp=alp.start, iBet=bet.start, lambda0=lambda0,
                  maxiter=maxiter, tol=tol, estLam=estLam, monitor=monitor)

  names(fit$betas)  <- x.names
  names(fit$alphas) <- z.names

  return( new("orth", call=call, betas=fit$betas, alphas=fit$alphas, variance=fit$variance,
                         num.iter=fit$num.iter, converge=fit$converge, lambda=fit$lambda, score=fit$score) )


}

