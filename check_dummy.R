library(rhdf5) # note: dev version to support complex
library(purrr) # mapping functions


f <- 'am.tmat.h5'

t <- rhdf5::h5read(f,"tmatrix", native = TRUE)

a <- rhdf5::h5ls(f, recursive = F, native = T, all = F)
a <- h5dump(f,read.attributes = TRUE,native=TRUE)
str(a)


check_h5 <- function(f){
  d <- h5dump(f, read.attributes = T, native=TRUE,drop=TRUE )
  
  a <- list(root = rhdf5::h5readAttributes(f,'/'),
            vacuum_wavelength = rhdf5::h5readAttributes(f,'vacuum_wavelength'),
            computation = rhdf5::h5readAttributes(f,'computation'),
            embedding = rhdf5::h5readAttributes(f,'embedding'),
            material = rhdf5::h5readAttributes(f,'scatterer/material'),
            geometry = rhdf5::h5readAttributes(f,'scatterer/geometry'))
  
  list(data = d, attributes = a)
}

library(dplyr)
om <- check_h5('am.tmat.h5')
lobstr::tree(om, show_attributes = F)
op <- check_h5('ap.tmat.h5')
# lobstr::tree(op, show_attributes = F)
oj <- check_h5('aj.tmat.h5')
# lobstr::tree(oj, show_attributes = F)
or <- check_h5('ar.tmat.h5')
lobstr::tree(or, show_attributes = F)

# str(oj$data$tmatrix)

library(diffobj)
# diffPrint(target=op, current=om)
# diffPrint(target=op, current=oj)
# diffPrint(target=op, current=or)
