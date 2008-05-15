`from.call.to.z` <-
function(object)
{
    z.index <- match( c("formula.z", "dataz"), names(object), 0 )
    z.match <- object[ c(1, z.index[1:2]) ]
    z.list <- as.list(z.match)
    z.list$formula <- z.list$formula.z
    z.list$data <- z.list$dataz
    call2 <- as.call(z.list)
    
    index2        <- match( c("formula", "data"), names(call2), 0 )
    z.match2      <- call2[ c(1,index2[1:2]) ]
    z.match2[[1]] <- as.name("model.frame")
    z.match2      <- eval( z.match2, parent.frame() )
    z.terms       <- attr(z.match2, "terms")
    
    z <- model.matrix( z.terms, z.match2, contrasts )
    rownames(z) <- NULL
    return(z)
}

