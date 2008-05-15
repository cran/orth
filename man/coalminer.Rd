\name{coalminer}
\alias{coalminer}
\docType{data}
\title{ Coal Miner Data}
\description{
    This data is from chapter 6 (Table 6.6) of McCullagh and Nelder (1989).  
}
\usage{data(coalminer)}
\format{
    A dataframe with 72 observations on the following 6 columns.
  \describe{
    \item{\code{id}}{Cluster id:  This identifies the age group}
    \item{\code{y}}{Binary variate indicating whether the subject has the event in question or not.  
                    Since each subject is measured on two-variates, breathlessness (B) and wheezing (W), 
                    it means that \eqn{Y_i} for a subject in the \eqn{i}-th age group can be thought of as a
                    2 by 1 vector of the form \eqn{Y_i = (Y_B, Y_W)} where \eqn{Y_B}, \eqn{Y_W = 0,1}.  For example,
                    \eqn{Y_B=0} and \eqn{Y_W=0} means no breathlessness and no wheezing.  }
    \item{\code{age}}{age group}
    \item{\code{ageCntrd}}{age centered at 5}
    \item{\code{wheeze}}{indicator for whether the response refers to wheezing or breathlessness.
                         1 if wheezing and 0 if breathlessness. }
    \item{\code{count}}{Cluster weights:  these number refer to how many subjects in a particular age group who
                        has the event (\eqn{y_B, y_W})}
  }
}

\references{
  P. McCullagh and J. A. Nelder. \emph{Generalized Linear Models (2nd ed.)}. \cr
  London: Chapman and Hall, 1989.
}

\keyword{datasets}
