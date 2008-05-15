`from.call.to.xyw2` <-
function(object)
{
    match.index <- match( c("formula", "data", "id", "weights"), names(object), 0 )
    x.match <- object[c(1,match.index[1:4])]
    x.match[[1]] <- as.name("model.frame")
    x.match <- eval(x.match, parent.frame())
    id <- as.vector( model.id(x.match) )
    model.terms <- attr(x.match, "terms")
    x <- model.matrix(model.terms, x.match, contrasts)
    w <- as.vector( model.weights(x.match) )
    
    n <- cls.size(id)
    w.index <- beginEnd(n)
    first.w <- w.index$first
    
    weights <- rep(NA, length(n)) 
    
    for( i in seq(1,length(n)) )
    {
        weights[i] <-  w[ first.w[i] ]   
    }
    rownames(x) <- NULL
    return( list(id=id,  X = x, y=matrix(model.response(x.match, "any"), ncol=1), w = weights ) )   
}

