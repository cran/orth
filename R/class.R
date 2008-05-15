setClass(
  "orth",
  representation(call="call", betas="vector",alphas="vector", variance="list",
                  num.iter = "numeric", converge="logical", lambda="numeric", score="list"),
  prototype(betas=c(NA), alphas=c(NA), variance=list(robust=NULL, naiveA=NULL,
             naiveB=NULL), num.iter=0, converge=FALSE, lambda=0 )
)

# The slot "score" will be used when we want diagnostics.


setClass(
  "summary.orth",
  representation(betas="matrix", alphas="matrix", betas.naive="matrix", alphas.naive="matrix")
)

setClass(
  "cluster.diag.orth",
  representation(
    dfbet = "matrix",
    dfalp = "matrix",
    h.bet = "matrix",
    h.alp = "matrix",
    cook.robust = "matrix",
    cook.naive  = "matrix"
  )
)


setClass(
    "obsvnl.diag.orth",
    representation(
        dfbet="matrix",
        dfalp="matrix",
        cook.robust = "matrix",
        cook.naive  = "matrix"
    )
)

#setGeneric("formula",
#    function(object, ...)
#    {
#        standardGeneric("formula")
#    }
#)

#setMethod("formula", "orth",
#    function(object,...)
#    {
#        object@call$formula
#    }
#)

setMethod("formula", signature(x = "orth"),
          function(x, ...) x@call[["formula"]])

# Updating an orth class.  The update will allow changing the formulas for
# the marginal mean model only.  If the user wants to update the association
# formula, then re-fit using a different formula.z.
#setGeneric("update",
#    function(object, ...)
#    {
#        standardGeneric("update")
#    }
#)

setMethod("update", "orth",
    function (object, ...)
    {
        .local <- function (object, formula., ..., evaluate = TRUE)
        {
            call <- object@call
            if (is.null(call)) { stop("need an object with call slot") }
            extras <- match.call(expand.dots = FALSE)$...
            if (!missing(formula.))
                call$formula <- update.formula(formula(object), formula.)
            if (length(extras) > 0)
            {
                existing <- !is.na(match(names(extras), names(call)))
                for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
                if (any(!existing)) {
                    call <- c(as.list(call), extras[!existing])
                    call <- as.call(call)
                }
            }
            if (evaluate){ eval(call, parent.frame()) }
            else call
        }
        .local(object, ...)
    }
)

#=============================================================
# This one belongs to Douglas Bates
#
#setMethod("update", signature(object = "orth"),
#          function(object, formula., ..., evaluate = TRUE)
#      {
#          call <- object@call
#          if (is.null(call))
#              stop("need an object with call slot")
#          extras <- match.call(expand.dots = FALSE)$...
#          if (!missing(formula.))
#              call$formula <- update.formula(formula(object), formula.)
#          if (length(extras) > 0) {
#              existing <- !is.na(match(names(extras), names(call)))
#              for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
#              if (any(!existing)) {
#                  call <- c(as.list(call), extras[!existing])
#                  call <- as.call(call)
#              }
#          }
#          if (evaluate)
#              eval(call, parent.frame())
#          else call
#      })
#=============================================================

setMethod(
  "show", "orth",
  function(object)
  {
    cat("Class: ", class(object), "\n")
    cat("\nCall:\n"); print(object@call)
    cat("\nMarginal Mean Parameters:\n"); print(object@betas)
    cat("\nAssociation Parameters (log-odds):\n"); print(object@alphas)
    cat("\nConvergence Status:  ", object@converge, "\n")
    cat("Lambda:                ", object@lambda, "\n")
    cat("Number of iterations:  ", object@num.iter, "\n")
  }
)


setMethod(
  "summary", "orth",
  function(object)
  {
    betas <- object@betas;    p1 <- length(betas)
    alphas <- object@alphas;  p2 <- length(alphas)

    robust  <- object@variance$robust
    cov.bet <- matrix(robust[1:p1, 1:p1], nrow=p1)
    cov.alp <- matrix(robust[(p1+1):(p1+p2), (p1+1):(p1+p2)], nrow=p2)
    
    var.bet <- diag( cov.bet)
    var.alp <- diag(cov.alp)
    
    se.bet  <- sqrt( var.bet )
    se.alp  <- sqrt( var.alp )

    se.bet.model <- sqrt( diag(object@variance$naiveB) )
    se.alp.model <- sqrt( diag(object@variance$naiveA) )

    chi.bet <- (betas)^2/var.bet
    chi.alp <- (alphas)^2/var.alp
    p.bet   <- 1-pchisq(chi.bet, 1)
    p.alp   <- 1-pchisq(chi.alp,1)

    chi.bet.model <- (betas)^2/(se.bet.model)^2
    chi.alp.model <- (alphas)^2/(se.alp.model)^2
    p.bet.model   <- 1-pchisq(chi.bet.model,1)
    p.alp.model   <- 1-pchisq(chi.alp.model,1)

    coef.bet <- cbind(betas, se.bet, chi.bet, p.bet)
    dimnames(coef.bet)[[2]] <- c('  Estimate', '   Std. Error',  '   Chi Square',  '     Pr(>Chi)')
    coef.alp <- cbind(alphas, se.alp, chi.alp, p.alp)
    dimnames(coef.alp)[[2]] <- c('  Estimate', '   Std. Error',  '   Chi Square',  '     Pr(>Chi)')

    coef.bet.model <- cbind(betas, se.bet.model, chi.bet.model, p.bet.model)
    dimnames(coef.bet.model)[[2]] <- c('  Estimate', '   Std. Error',  '   Chi Square',  '     Pr(>Chi)')
    coef.alp.model <- cbind(alphas, se.alp.model, chi.alp.model, p.alp.model)
    dimnames(coef.alp.model)[[2]] <- c('  Estimate', '   Std. Error',  '   Chi Square',  '     Pr(>Chi)')

    my.summary <- new("summary.orth", betas=coef.bet, alphas=coef.alp, betas.naive=coef.bet.model, alphas.naive=coef.alp.model)
    return(my.summary)

  }
)


setMethod(
  "show", "summary.orth",
  function(object)
  {
    cat("\nClass: ", class(object), "\n")
    cat("\nSummary values based on robust covariance.  Those interested in\n")
    cat("model-based covariance may use the 'SUMMARY()' method on this \n")
    cat("summary.orth class.\n")

    cat("\nMarginal Mean Parameters:\n"); print(object@betas)
    cat("\n\nAssociation Parameters (log-odds):\n"); print(object@alphas); cat("\n")
  }
)

setMethod(
  "summary", "summary.orth",
  function(object)
  {
    cat("\nClass: ", class(object), "\n")
    cat("\nBased on robust standard errors.\n")
    cat("\nMarginal Mean Parameters:\n"); print(object@betas)
    cat("\nAssociation Parameters (log-odds):\n"); print(object@alphas)
    cat("\n\nBased on mode-based standard errors.\n")
    cat("\nMarginal Mean Parameters:\n"); print(object@betas.naive)
    cat("\nAssociation Parameters (log-odds):\n"); print(object@alphas.naive); cat("\n")
  }
)


## Begin:  Accessor function to extract dfbeta from the cluster.diag.orth class ##.
setGeneric("df.bet",
    function(object)
    {
        standardGeneric("df.bet")
    }
)

setMethod(
    "df.bet", "cluster.diag.orth",
    function(object)
    {
        object@dfbet
    }
)
## End:  Accessor function to extract dfbeta from the cluster.diag.orth class ##.

## Begin:  Accessor function to extract dfalpha from the cluster.diag.orth class ##.
setGeneric("df.alp",
    function(object)
    {
        standardGeneric("df.alp")
    }
)

setMethod(
    "df.alp", "cluster.diag.orth",
    function(object)
    {
        object@dfalp
    }
)

## End:  Accessor function to extract dfalpha from the cluster.diag.orth class ##.

## Begin:  Accessor function to extract leverage for beta from the cluster.diag.orth class ##.
setGeneric("lev.bet",
    function(object)
    {
        standardGeneric("lev.bet")
    }
)

setMethod(
    "lev.bet", "cluster.diag.orth",
    function(object)
    {
        object@h.bet
    }
)

## End:  Accessor function to extract leverage for beta from the cluster.diag.orth class ##.


## Begin:  Accessor function to extract leverage for alpha from the cluster.diag.orth class ##.
setGeneric("lev.alp",
    function(object)
    {
        standardGeneric("lev.alp")
    }
)

setMethod(
    "lev.alp", "cluster.diag.orth",
    function(object)
    {
        object@h.alp
    }
)

## End:  Accessor function to extract leverage for alpha from the cluster.diag.orth class ##.


## Begin:  Accessor function to extract robust Cook's distance from the cluster.diag.orth class ##.
setGeneric("cook.robust",
    function(object)
    {
        standardGeneric("cook.robust")
    }
)

setMethod(
    "cook.robust", "cluster.diag.orth",
    function(object)
    {
        object@cook.robust
    }
)

## End:  Accessor function to extract robust Cook's distance from the cluster.diag.orth class ##.


## Begin:  Accessor function to extract naive Cook's distance from the cluster.diag.orth class ##.
setGeneric("cook.naive",
    function(object)
    {
        standardGeneric("cook.naive")
    }
)

setMethod(
    "cook.naive", "cluster.diag.orth",
    function(object)
    {
        object@cook.naive
    }
)

## End:  Accessor function to extract naive Cook's distance from the cluster.diag.orth class ##.


## Function overloading for observational level diagnostics                    ##.
setMethod(
    "df.bet", "obsvnl.diag.orth",
    function(object)
    {
        object@dfbet
    }
)

setMethod(
    "df.alp", "obsvnl.diag.orth",
    function(object)
    {
        object@dfalp
    }
)

setMethod(
    "cook.robust", "obsvnl.diag.orth",
    function(object)
    {
        object@cook.robust
    }
)

setMethod(
    "cook.naive", "obsvnl.diag.orth",
    function(object)
    {
        object@cook.naive
    }
)

setMethod(
    "lev.bet", "obsvnl.diag.orth",
    function(object)
    {
        stop("Not defined for class 'obsvnl.diag.orth' ")
    }
)


setMethod(
    "lev.alp", "obsvnl.diag.orth",
    function(object)
    {
        stop("Not defined for class 'obsvnl.diag.orth' ")
    }
)
