#!/usr/bin/env Rscript
# there seems to be a big problem
# actually it wasn't so big. plinkbin was coming out as object not found,
# but actually it *was* found in the first function, but I hadn't definted it in the second.
# That was the problem.
rawToMaster <- function(a=3){

  if( a== 5){
  plinkbin <- Sys.which("plink1")
  recodeCom <- paste(plinkbin,
                       " --noweb --silent --horse --file ", sep = "")    

  cat(c(recodeCom, "\n"))
  }
}
a<-5
rawToMaster(a)
rawToMaster()
