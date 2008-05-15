`cls.df2.i` <-
function(DA, rootVEE, sqe1, sqe2, R, A, invB, CC, c1, c2, mi,
                       p1, p2, orthVar)
{
    # associated with BETAs
    mm1   <- CC %*% sqe1
    ai1A  <- invB %*% A
    ai1m1 <- invB %*% mm1

    inv.cls.b <- invBig(ai1A, ai1m1, mm1, A, 1, p1)

    U.i.b <- t(CC) %*% inv.cls.b$aic
    dfbetclsi <- orthVar$naiveB %*% U.i.b

    # associated with ALPHAs
    mm2    <- DA %*% sqe2
    colmm2 <- mm2/rootVEE

#    print("mm2 is "); print(mm2);
#    print("colmm2 is"); print(colmm2)

    ai2R   <- (c1*(R/rootVEE) + c2*sum(R/rootVEE)) / rootVEE
    ai2m2  <- ( c1*colmm2 + c2 * ( matrix(rep(1,mi),ncol=1) %*% matrix(apply(colmm2,2,sum),nrow=1) ) ) / rootVEE

#    print("ai2R before invBig()"); print(ai2R)
#    print("ai2m2 before invBig()"); print(ai2m2)
#    print("const1 is"); print(c1)
#    print("const2 is"); print(c2)

    inv.cls.a <- invBig(ai2R, ai2m2, mm2, R, 1, p2)

#    print("ai2R"); print(inv.cls.a$aic)
#    print("ai1m2"); print(inv.cls.a$aim)

    U.i.a     <- t(DA) %*% inv.cls.a$aic
    dfalpclsi <- orthVar$naiveA %*% U.i.a

  return( list( dfBeta.cls = dfbetclsi, dfAlpha.cls = dfalpclsi ) )

}

