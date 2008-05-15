\name{cook.naive-methods}
\docType{methods}
\alias{cook.naive}
\title{ Extracts Naive Cook's Distance.}
\description{
    This function works on both objects of class (S4)
    \emph{cluster.diag.orth} and \emph{obsvnl.diag.orth}. Using 
    these as inputs, it computes Cook's distance based on naive
    covariance. 
}
\section{Methods}{
\describe{

\item{object = "cluster.diag.orth"}{ Extracts cluster level naive (model-based) Cook's distance }

\item{object = "obsvnl.diag.orth"}{ Extracts observation level naive (model-based) Cook's distance }
}}

\seealso{\code{\link{cook.robust}}}

\keyword{methods}
