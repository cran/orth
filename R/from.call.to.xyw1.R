`from.call.to.xyw1` <-
function(object)
{
    match.index <- match( c("formula", "data", "id"), names(object), 0 )
    x.match <- object[c(1,match.index[1:3])]
    x.match[[1]] <- as.name("model.frame")
    x.match <- eval(x.match, parent.frame())
    id <- as.vector( model.id(x.match) )
    model.terms <- attr(x.match, "terms")
    x <- model.matrix(model.terms, x.match, contrasts)
    w <- rep( 1, length( cls.size(id) ) )
    rownames(x) <- NULL
    return( list(id=id,  X = x, y=matrix(model.response(x.match, "any"), ncol=1), w = w ) )   
}

