\name{orth-class}
\docType{class}
\alias{orth-class}
\alias{summary,orth-method}
\alias{update,orth-method}
\alias{show,orth-method}
\alias{formula,orth-method}
\title{Class "orth"  }
\description{ An S4 class.  The class is returned from the "orth()" function}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("orth", ...)}.
    However, the user can pretty much ignore this since it is returned
    automatically from fitting a regression model.
}
\section{Slots}{
	 \describe{
    \item{\code{call}:}{Object of class \code{"call"}  }
    \item{\code{betas}:}{Object of class \code{"vector"}.  Estimates of the mean parameters.  }
    \item{\code{alphas}:}{Object of class \code{"vector"}.  Estimates of the association parameters }
    \item{\code{variance}:}{Object of class \code{"list"}.  Contains covariance matrices }
    \item{\code{num.iter}:}{Object of class \code{"numeric"}. Number of iterations  }
    \item{\code{converge}:}{Object of class \code{"logical"}.  Convergence status}
    \item{\code{lambda}:}{Object of class \code{"numeric"}.  Value of the correlation parameter if estimated }
    \item{\code{score}:}{Object of class \code{"list"}.   }
  }
}
\section{Methods}{
    \describe{
        \item{formula}{Extracts formula from the orth class.  User need not be concerned with this}
        \item{update}{Updates orth class when new model formula is provided}
        \item{show}{Prints the results from fitting regression model}
        \item{summary}{Gives summary from fitting regression model}
    }
}

\keyword{classes}