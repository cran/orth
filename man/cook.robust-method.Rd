\name{cook.robust-methods}
\docType{methods}
\alias{cook.robust}
\title{ Extracts Robust Cook's Distance}
\description{
    This function works on both objects of class (S4)
    \bold{cluster.diag.orth} and \bold{obsvnl.diag.orth}. Using 
    these as inputs, it computes Cook's distance based on robust
    covariance. 
}
\section{Methods}{
\describe{

\item{object = "cluster.diag.orth"}{ Extracts cluster level robust Cook's distance }

\item{object = "obsvnl.diag.orth"}{ Extracts observation level robust Cook's distance }
}}
\seealso{\code{\link{cook.naive}}}

\keyword{methods}
