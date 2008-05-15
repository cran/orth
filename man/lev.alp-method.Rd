\name{lev.alp-methods}
\docType{methods}
\alias{lev.alp}

\title{ Extracts Leverage Values For The Association Parameters }
\description{
    This functions extract cluster level leverage.  Works only on objects of
    class "cluster.diag.orth".  If this function is applied to an object of 
    class "obsvnl.diag.orth", a complaint is registered.
}
\section{Methods}{
    \describe{
    \item{object = "cluster.diag.orth"}{ Cluster level leverage }
    
    \item{object = "obsvnl.diag.orth"}{ Gives an error message }
    }
}

\seealso{\code{\link{lev.bet}}}
\keyword{methods}


