`vee.dAdB.R2` <-
function(iI2.obj, vabr.obj, ni)
{
#  print("vabr.obj is"); print(vabr.obj)
  vee <- vabr.obj$VEE
  
  if (ni > 2)
  {
    veeDiag <- diag(vee)
#    print("veeDiag is"); print(veeDiag)
#    print("vee dimension is"); print(dim(veeDiag))
#    print("iI2.obj$comp is"); print(iI2.obj$comp)
#    print("vabr.obj$VEE[iI2.obj$comp is"); print(vabr.obj$VEE[iI2.obj$comp])
    vee.inv <- 1/vabr.obj$VEE[iI2.obj$comp]
    if (length(vee.inv) > 1) { veeDcomp.inv <- diag(vee.inv) }
    else { veeDcomp.inv <- matrix(1/vee.inv,nrow=1) }
#    print("vee.inv is"); print(vee.inv)
#    print("veeDcomp.inv is"); print(veeDcomp.inv)
    off.set      <- veeDiag[iI2.obj$obs, iI2.obj$comp] %*% veeDcomp.inv
    DAt          <- vabr.obj$DA[iI2.obj$obs, ] - off.set %*% vabr.obj$DA[iI2.obj$comp,]
    Rt           <- vabr.obj$R[iI2.obj$obs] - off.set %*% vabr.obj$R[iI2.obj$comp]
    VEEt         <- veeDiag[iI2.obj$obs, iI2.obj$obs] - off.set%*%veeDiag[iI2.obj$comp, iI2.obj$obs]
    sqVEEt       <- sqrt(diag(VEEt)) 
  }
  else if (ni == 2)
  {
    DAt <- vabr.obj$DA[iI2.obj$obs, ]
    Rt  <- vabr.obj$R[iI2.obj$obs]
    VEEt <- vee
    sqVEEt <- sqrt(VEEt) 
#    VEEt <- veeDiag[iI2.obj$obs, iI2.obj$obs]
#    sqVEEt <- sqrt(diag(VEEt))    
  }
  return( list(DA = DAt, sqVEE = sqVEEt, R=Rt) )
}

