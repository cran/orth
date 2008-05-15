\name{summary.orth-class}
\docType{class}
\alias{summary.orth-class}
\alias{show,summary.orth-method}
\alias{summary,summary.orth-method}
\title{Class summary.orth}
\description{An S4 class that contains the summary of an orth object.  This class
is created as a result of running the code \code{summary()} on an orth class. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("summary.orth", ...)}.
However, the user should never do this.
}
\section{Slots}{
	 \describe{
    \item{\code{betas}:}{Object of class \code{"matrix"}}
    \item{\code{alphas}:}{Object of class \code{"matrix"} }
    \item{\code{betas.naive}:}{Object of class \code{"matrix"} }
    \item{\code{alphas.naive}:}{Object of class \code{"matrix"} }
  }
}
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object = "summary.orth")}: Provides the summary
            from model fitting}
        \item{summary}{\code{signature(object = "summary.orth")}: Provides the summary
            of the summary.  This is done because a summary on an orth object gives everything
            in terms of robust covariance.  Because we do not want to clutter the output, the
            user has to perform a summary on the summary if he/she wants to see things in terms
            of naive (or model-based) standard errors.}
    }
}

\keyword{classes}
