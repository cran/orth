\name{orth}
\alias{orth}

\title{ Logistic Regression For Correlated Binary Data By Way of Orthogonalized Residuals or ALR }
\description{
    Fits mean and association models for correlated binary data based
    on orthogonalized residuals.  The link function for the mean model
    is the logit link and the link function for the association model is
    the log-odds ratio.  
}

\usage{
orth(formula, data, weights , formula.z, dataz, id, contrasts = NULL, 
     alp.start = NULL, bet.start = NULL, lambda0 = 0, maxiter = 20, tol = 1e-06, 
     estLam = FALSE, monitor = FALSE, ...)
}

\arguments{
  \item{formula}{ Two-sided formula for mean model }
  \item{data}{ Data where the response and the covariates for the mean model are stored }
  \item{weights}{ Cluster weights }
  \item{formula.z}{ One-sided formula for the association model }
  \item{dataz}{ Data where the covariates for the association model are stored }
  \item{id}{ Cluster id }
  \item{contrasts}{ model contrasts to be used }
  \item{alp.start}{ Starting values for association parameters }
  \item{bet.start}{ Starting values for mean parameters }
  \item{lambda0}{ Starting values for lambda }
  \item{maxiter}{ Maximum number of iterations }
  \item{tol}{ Size of precision }
  \item{estLam}{ logical: 1 if we want to estimate lambda and 0 if no }
  \item{monitor}{ logical: 1 if we want to monitor the progress of the estimation }
  \item{\dots}{  }
}

\details{If \code{estLam = FALSE} then \code{lambda} is set to 0 and  
    alternating logistic regression based on orthogonalized residuals is fit to the data.
}

\section{Warning}{
    If your data set contains \bold{many} clusters with cluster sizes of \bold{20 or more}\cr
    then you will experience a noticeably long time for the estimation routine to finish.
}

\value{
  Returns an S4 class of class "orth" with the following slots
  \item{call }{ The call of the ``orth'' class }
  \item{betas }{ The mean parameters -- belongs to vector class }
  \item{alphas }{ The association parameters -- belongs to vector class }
  \item{variance }{ A list containing both the naive and robust covariance matrices }
  \item{num.iter }{ Number of iterations performed by the procedure }
  \item{converge }{ Logical:  1 if procedure converge }
  \item{lambda }{ Estimated lambda if \code{estLam=TRUE} is specified in \code{orth} }
  \item{score }{ A list containing information about the score vector and Fisher information }
}

\references{ R. C. Zink. \emph{Correlated Binary Regression Using Orthogonalized Residuals}.\cr
             PhD thesis, University of North Carolina, Chapel Hill, 2003. }
             
\keyword{models}
\keyword{regression}


\examples{

library(orth)

## EXAMPLE 1:  Six cities study.
##   
## Model formula for the mean:         1 + SMOKE + age
## Model formula for the association : 1 + SMOKE + agediff
##
## The link function for the mean model is the logit link.
## For the association model, we are modelling the log of the pairwise odds ratio 
## within each cluster.  
data(six.city)

## Creating the Z-matrix for the association model
n = as.vector( table( as.factor(six.city$id) ) )
last <- cumsum(n);
first <- last - n + 1;
z.6city <- NULL
for ( i in 1:length(n) )
{
    n.i <- n[i]
    id.i <- six.city$id[ first[i] ]
    age.i <- six.city$age[ first[i]:last[i] ]
    smoke.i <- six.city$smoke[ first[i]:last[i] ]
    
    l <- 1
    
    if (n.i == 1) {z.i <- cbind(0,0,0)}
    else
    {
        ## Note: ch2(m) = m(m - 1)/2
        id.i <- rep(id.i, choose(n.i, 2))
        z.i1 <- rep(NA, choose(n.i, 2) )
        z.i2 <- rep(NA, choose(n.i, 2) )
        for( j in seq(1, n.i - 1) )
        {
            for( k in seq(j+1, n.i) )
            {
                z.i1[l] <- abs( age.i[j] - age.i[k] )
                z.i2[l] <- smoke.i[1]
                l <- l+1
            }
        }
        z.i <- cbind(id.i, z.i1, z.i2)
    }

    z.6city <- rbind(z.6city, z.i) 
}

## Remove any row of zeros.  These correspond to cluster size 1 which do  ##
## not exists in the association matrix.                                  ##
z.6city <- data.frame(z.6city)
names(z.6city) <- c("id", "agediff", "smoke")
z.6city <- z.6city[(z.6city$id != 0), ]

orth1 =  orth(y~age+smoke,data=six.city, weights=count, formula.z=~agediff + smoke, dataz=z.6city, 
              id=id, maxiter=20, tol=0.001, estLam=TRUE)
orth1
              
o1.s <- summary(orth1)
o1.s
              
## Fits a mean model without intercept but same association model.
orth2 = update(orth1, . ~ . -1)

## Computes cluster diagnostics ##
cd = diag.cls(orth1)

## Cluster level Cook's D based on robust and naive covariance matrices
cook.robust(cd)
cook.naive(cd)

## Cluster level dbetas and dalphas
df.bet(cd)
df.alp(cd)

## Cluster level leverage for betas and alphas
lev.alp(cd)
lev.bet(cd)

## Computes observation level diagnostics.  We do not recomend computing
## these values if you have a bunch of clusters that have cluster size
## of 7 or more.  Side effects include high blood pressure and excessive
## swearing caused by excessive waiting time for the observation level
## diagnostics routine to finish.  If you have taken pills for these, go
## for it.  
od = diag.obs(orth1)
cook.robust(od)
cook.naive(od)
df.bet(od)
df.alp(od)


#------------------------------------------------------------------------------#
## EXAMPLE 2:  Urinary incontinence study.
##   
## Model formula for the mean: 1 + FEMALE + age + toilet + toilet2 + severe 
#                                + dayacc + dayacc2
## Model formula for the association : sameMD 
##
## The link function for the mean model is the logit link.
## For the association model, we are modelling the log of the pairwise odds ratio 
## within each cluster.
  
data(ui)

ui$toilet2 = ui$toilet^2
ui$dayacc2 = ui$dayacc^2

## Creates the Z matrix for the association model.

n = as.vector( table(as.factor(ui$pract_id)) )
last = cumsum(n)
first = last - n + 1
z.ui <- NULL

for ( i in 1:length(n) )
{
    n.i <- n[i]
    practID.i <- ui$pract_id[ first[i]]
    doctrID.i <- ui$doct_id[ first[i]:last[i] ]
    
    l <- 1

    if (n.i == 1) { z.i <- cbind( 0,0) }
    else
    {
        practID.i <- rep( practID.i, choose(n.i, 2) )
        z.i1 <- rep( NA, choose(n.i, 2) )
        for( j in seq(1, n.i-1) )
        {
            for (k in seq(j+1, n.i) )
            {
                z.i1[l] <- ( doctrID.i[j] == doctrID.i[k] )
                l <- l + 1
            }
        }
        z.i <- cbind(practID.i, z.i1)
    }
    z.ui <- rbind(z.ui, z.i)
}
z.ui <- data.frame(z.ui)
## Remove any row of zeros.  These correspond to cluster size 1 which do  ##
## not exists in the association matrix.                                  ##
names(z.ui) <- c("PracticeID", "sameMD")
z.ui <- z.ui[(z.ui$PracticeID != 0),]
z.ui$sameMD <- as.factor(z.ui$sameMD)

orth3 = orth(bothered~factor(female)+age+toilet+toilet2+severe+dayacc+dayacc2,data = ui,
             formula.z=~-1 + sameMD, dataz=z.ui, id=pract_id, maxiter=55, tol=0.001, estLam=TRUE)


## This data contains several clusters with huge cluster size.  So observation
## level diagnostics are not recommended.  Cluster level diagnostics are fast
## to compute.
cd = diag.cls(orth3)
df.bet(cd)
df.alp(cd)
lev.bet(cd)
lev.alp(cd)
cook.robust(cd)
cook.naive(cd)



#------------------------------------------------------------------------------#
## EXAMPLE 3:  Coalminer data from McCullagh and Nelder (1989).

## Model formula for the marginal mean:  1 + age
##
## Model formula for the marginal association : 1 + age
##
## The link for the mean model is the LOGIT link.  For the association, we
## are modelling the pairwise log-odds ratio of an event.

data(coalminer)

n = as.vector( table(as.factor(coalminer$id)) )
last = cumsum(n)
first = last - n + 1

age.i <- rep( NA, length(n) )

for(i in 1:length(n) )
{
    age.i[i] <- coalminer$ageCntrd[first[i]]
}


z.coalminer <- data.frame( id = as.vector(unique(coalminer$id) ), age=age.i )

orth4 = orth(y~-1 + factor(wheeze) + factor(wheeze):ageCntrd, data=coalminer, weights=count,
             formula.z = ~age, dataz=z.coalminer, id=id, maxiter=10, tol=0.001, estLam=FALSE)
             
## Looks at estimates and standard errors based on robust covariance
summary(orth4)

## Looks at estimates and standard errors based on model-based covariance
summary( summary(orth4) )

## Cluster level diagnostics 
cd = diag.cls(orth4)
df.bet(cd)
df.alp(cd)
lev.bet(cd)
lev.alp(cd)
cook.robust(cd)
cook.naive(cd)

#------------------------------------------------------------------------------#

}
