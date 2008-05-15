\name{cluster.diag.orth-class}
\docType{class}
\alias{cluster.diag.orth-class}
\alias{cook.naive,cluster.diag.orth-method}
\alias{cook.robust,cluster.diag.orth-method}
\alias{df.alp,cluster.diag.orth-method}
\alias{df.bet,cluster.diag.orth-method}
\alias{lev.alp,cluster.diag.orth-method}
\alias{lev.bet,cluster.diag.orth-method}
\title{Class "cluster.diag.orth" }
\description{ An S4 class that contains the cluster level deletion diagnostics.  It
    is created as a result of running the function \code{diag.cls(orth.obj)}	}
    
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("cluster.diag.orth", ...)}.
    However, the user rarely (if ever) has to do this since the function
    \code{diag.cls(orth.obj)} is used to create this class. 
}

\section{Slots}{
	 \describe{
    \item{\code{dfbet}:}{Object of class \code{"matrix"}: Contains cluster level \bold{dfbeta}  }
    \item{\code{dfalp}:}{Object of class \code{"matrix"}: Contains cluster level \bold{dfalpha}  }
    \item{\code{h.bet}:}{Object of class \code{"matrix"}: Contains cluster level \bold{leverage} for the mean parameters  }
    \item{\code{h.alp}:}{Object of class \code{"matrix"}: Contains cluster level \bold{leverage} for the association parameters  }
    \item{\code{cook.robust}:}{Object of class \code{"matrix"}: Contains cluster level \bold{Cook's distance} based on
        robust covariance   }
    \item{\code{cook.naive}:}{Object of class \code{"matrix"}: Contains cluster level \bold{Cook's distance} based on
        naive covariance  }
  }
}
\section{Methods}{
  \describe{
    \item{cook.naive}{\code{signature(object = "cluster.diag.orth")}: Cook's distance
        based on naive covariance matrix }
    \item{cook.robust}{\code{signature(object = "cluster.diag.orth")}: Cook's distance
        based onrobust covariance matrix }
    \item{df.alp}{\code{signature(object = "cluster.diag.orth")}: Computes DFALPHA }
    \item{df.bet}{\code{signature(object = "cluster.diag.orth")}: Computes DFBETA }
    \item{lev.alp}{\code{signature(object = "cluster.diag.orth")}: Computes leverage
        for the association parameters }
    \item{lev.bet}{\code{signature(object = "cluster.diag.orth")}: Computes leverage
        for the mean parameters }
	 }
}

\seealso{
	\code{\link{obsvnl.diag.orth-class}}
}

\keyword{classes}
