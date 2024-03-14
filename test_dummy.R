library(rhdf5) # note: dev version to support complex
library(uuid) # uuid
library(glue) # string interpolation
library(purrr) # mapping functions


## dummy data

# possibly multiple wavelengths
wavelength <- seq(400, 800, by=50)
# wavelength <- 633.0
Nl <- length(wavelength)
# dummy 30x30 matrix values for each wavelength
tmatrix <- array(1:900 + 1i*(1:900), c(30,30,Nl))

modes <- list('l' = 1:30,'m' = 1:30,
              'polarization' =  rep(c('electric','magnetic'),15))

# dummy 'analytical zeros' for e.g. EBCM methods
zeros <- expand.grid(q=seq(1,30,by=2), qp = seq(1,30,by=2))

embedding <- list('relative_permeability' = 1.0, 'relative_permittivity' = 1.33^2)
materials <- list('embedding' = embedding, 
                  'Au' = list('relative_permeability' = 1.0, 'relative_permittivity' = -11.4+1.181i))

geometry <- list('shape' = 'spheroid','radiusxy' = 20.0, 'radiusz' = 40.0)

# details about computation, including full script
computation <- list('method' ='EBCM',
                    'software' = 'SMARTIES',
                    'version' = '1.1',
                    'Ntheta' = 40, 
                    'accuracy' = 1e-10, 
                    'analytical_zeros' = zeros,
                    'script' = glue::glue_collapse(readLines('test_dummy.R'), sep = "\n"))

uuid <- uuid::UUIDgenerate()

## create file
f <- 'ar.tmat.h5'

unlink(f)
h5createFile(f)
h5closeAll()

## write file
h5write(tmatrix, file=f, name="/tmatrix", native=TRUE) # store rowwise
# Using native = TRUE increases HDF5 file portability between programming languages
h5write(modes, file=f, name='/modes')
h5write(wavelength, file=f, name="/vacuum_wavelength")
h5write(embedding, file=f, name="/embedding")
h5write(materials, file=f, name="/materials")
h5write(geometry, file=f, name="/geometry")
h5write(computation, file=f, name='/computation')
h5write(uuid, file=f, name='/uuid')

## write attributes

write_attributes <- function(object, names, attributes, 
                             type=c("root","dataset","group")){
  
  switch(type,
         "root" = {
           stopifnot(object == "/")
           purrr::walk2(names, attributes, \(n,a) h5writeAttribute(n, fid, a))
         },         
         "dataset" = {
           did <- H5Dopen(fid, object)
           purrr::walk2(names, attributes, \(n,a) h5writeAttribute(n, did, a))
           H5Dclose(did)
         },
         "group" = {
           did <- H5Gopen(fid, object)
           purrr::walk2(names, attributes, \(n,a) h5writeAttribute(n, did, a))
           H5Gclose(did)
         },
         error("needs to be a group, dataset, or root"))
  
}


fid <- H5Fopen(f)
# ## root attributes
write_attributes("/",
                 names = list("name",
                              "created_with",
                              "storage_format_version",
                              "description", 
                              'keywords'), 
                 attributes = list("Au prolate spheroid in water",
                                   "rhdf5",
                                   glue::glue_collapse(H5get_libversion(),"."),
                                   "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids", 
                                   'gold, spheroid, ebcm'), 
                 type = "root")

## object attributes
write_attributes("tmatrix",
                 names=c("units","created_with"),
                 attributes = c("nm","rhdf5"), 
                 type = "dataset")

write_attributes("vacuum_wavelength",
                 names=c("unit"), attributes = c("nm"), 
                 type = "dataset")

write_attributes("uuid",
                 names=c("version"), attributes = c("4"), 
                 type = "dataset")

write_attributes("geometry",
                 names=c("name"), attributes = c("prolate spheroid"), 
                 type = "group")

## test 
h5ls(fid)

H5Fclose(fid)

