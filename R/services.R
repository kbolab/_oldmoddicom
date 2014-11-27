# ------------------------------------------------
# SV.getPointPlaneDistance
# ------------------------------------------------
SV.getPointPlaneDistance<-function(Punto,Piano) {
  return(abs(Piano[1]*Punto[1]+Piano[2]*Punto[2]+Piano[3]*Punto[3]+Piano[4])/sqrt(Piano[1]^2+Piano[2]^2+Piano[3]^2))
}
# ------------------------------------------------
# SV.get3DPosFromNxNy
# ------------------------------------------------
SV.get3DPosFromNxNy<-function(Nx,Ny,oM) {
  return(xy<-t(oM%*%c(Nx,Ny,0,1)))
}
# ------------------------------------------------
# SV.getPlaneEquationBetween3Points
# ------------------------------------------------
SV.getPlaneEquationBetween3Points<-function(Pa,Pb,Pc) {
  ac<-(Pb[2]-Pa[2])*(Pc[3]-Pa[3])-(Pc[2]-Pa[2])*(Pb[3]-Pa[3])
  bc<-(Pb[3]-Pa[3])*(Pc[1]-Pa[1])-(Pc[3]-Pa[3])*(Pb[1]-Pa[1])
  cc<-(Pb[1]-Pa[1])*(Pc[2]-Pa[2])-(Pc[1]-Pa[1])*(Pb[2]-Pa[2])
  dc<--(ac*Pa[1]+bc*Pa[2]+cc*Pa[3])
  return(c(ac,bc,cc,dc))
}