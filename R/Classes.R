#' dvmmatrix S4 class
#' This class handles the matrices containing the DVHs
setClass("dvhmatrix", 
         representation(
           dvh="matrix",                # matrix of DVHs, 1st column is dose, other columns are the volume
           dvh.type="character",        # cumulative or differential DVH
           vol.distr="character",       # absolute or relative distribution of the volume
           volume="numeric"             # vector of absolute volumes of DVH
           )
         )
