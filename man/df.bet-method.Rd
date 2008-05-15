\name{df.bet-methods}
\docType{methods}
\alias{df.bet}
\title{ Extracts Deletion Diagnostics For Mean Parameters}
\description{
    Extracts the deletion diagnostics for the mean parameters
    which we call "dfbetas".  Works on both objects of 
    class cluster.diag.orth and obsvnl.diag.orth.
}
\section{Methods}{
    \describe{
    \item{object = "cluster.diag.orth"}{ Extracts cluster level dfbetas }
    
    \item{object = "obsvnl.diag.orth"}{ Extracts observational level dfalphas }
    }
}
\keyword{methods}
\seealso{\code{\link{df.alp}}}
