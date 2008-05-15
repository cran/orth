\name{obsvnl.diag.orth-class}
\docType{class}
\alias{obsvnl.diag.orth-class}
\alias{cook.naive,obsvnl.diag.orth-method}
\alias{cook.robust,obsvnl.diag.orth-method}
\alias{df.alp,obsvnl.diag.orth-method}
\alias{df.bet,obsvnl.diag.orth-method}
\alias{lev.alp,obsvnl.diag.orth-method}
\alias{lev.bet,obsvnl.diag.orth-method}
\title{Class "obsvnl.diag.orth" }
\description{	An S4 class that contains the observational level deletion diagnostics.  It
    is created as a result of running the function \code{diag.cls(orth.obj)} }
    
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("obsvnl.diag.orth", ...)}.
    However, the user rarely (if ever) has to do this since the function
    \code{diag.obs(orth.obj)} is used to create this class. 
}

\section{Slots}{
	 \describe{
    \item{\code{dfbet}:}{Object of class \code{"matrix"}  }
    \item{\code{dfalp}:}{Object of class \code{"matrix"}  }
    \item{\code{cook.robust}:}{Object of class \code{"matrix"}  }
    \item{\code{cook.naive}:}{Object of class \code{"matrix"}  }
  }
}
\section{Methods}{
  \describe{
    \item{cook.naive}{\code{signature(object = "obsvnl.diag.orth")}: Cook's distance
        based on naive covariance matrix }
    \item{cook.robust}{\code{signature(object = "obsvnl.diag.orth")}: Cook's distance
        based onrobust covariance matrix }
    \item{df.alp}{\code{signature(object = "obsvnl.diag.orth")}: Computes DFALPHA }
    \item{df.bet}{\code{signature(object = "obsvnl.diag.orth")}: Computes DFBETA }
    \item{lev.alp}{\code{signature(object = "obsvnl.diag.orth")}: Computes leverage
        for the association parameters }
    \item{lev.bet}{\code{signature(object = "obsvnl.diag.orth")}: Computes leverage
        for the mean parameters }
	 }
}

\seealso{

	\code{\link{cluster.diag.orth-class} }
}

\keyword{classes}