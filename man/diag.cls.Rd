\name{diag.cls}
\alias{diag.cls}

\title{ Computes Cluster Level Diagnostics }
\description{
    Computes cluster level deletion diagnostics using as input
    an object of \bold{class} "orth" which is created by running the "orth" \bold{function}.
}

\usage{
diag.cls(object)
}

\arguments{
  \item{object}{ An object of class "orth" }
}


\value{
  The function returns an object of class "cluster.diag.orth".
}

\references{
    Preisser, J. S., Perin, J., Deletion diagnostics for marginal mean and correlation \cr
    model parameters in estimating equations, \emph{Statistics and Computing}, 17:381-393, 2007.
}

\keyword{classes}


\seealso{ \code{\link{cluster.diag.orth-class}}, \code{\link{diag.obs}} }

