\name{ui}
\alias{ui}
\docType{data}

\title{ Urinary Incontinence Study}
\description{
  Data from urinary incontinence study analyzed by Preisser and Qaqish (1999).
}
\usage{data(ui)}

\format{
  A data frame with 137 observations on the following 12 variables.
  \describe{
    \item{\code{pract_id}}{cluster id: used for identifying the practice}
    \item{\code{doct_id}}{doctor's id: used for identifying the doctors in the practice}
    \item{\code{pat_id}}{patient's id: used for identifying the patients in the practice}
    \item{\code{bothered}}{response:  whether patient's condition is bothersome}
    \item{\code{female}}{indicator for female}
    \item{\code{ageyrs}}{patient's age}
    \item{\code{age}}{ standardized age: (age - 76)/10 }
    \item{\code{weekacc}}{ patient report of number of leaking accidents they experienced in an average week }
    \item{\code{dayacc}}{ patient report of number of leaking accidents they experienced in an average day }
    \item{\code{severe}}{ severity of loss of urine: 1 if there is only moisture, 2 if the patient wet the underwear,
                          3 if the urine trickled down the thigh, and 4 if the patient wet the floor }
    \item{\code{toilet}}{ number of times patient frequents the toilet to urinate  (in a day) }
    \item{\code{mdage}}{age of the doctors in the practice}
  }
}

\references{
    Hammil, B. G., Preisser, J. S., A SAS/IML software program for GEE and regression diagnostics,
    \emph{Computational Statistics and Data Analysis}, 51: 1197 - 1212, 2006.\cr
    
    Preisser, J. S., Qaqish, B. F., Robust regression for clustered data with application to binary responses, 
    \emph{Biometrics}, 55: 574-579, 1999.   
}

\keyword{datasets}
