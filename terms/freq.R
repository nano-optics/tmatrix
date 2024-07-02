setwd("~/Documents/nano-optics/tmatrix/terms")

## generate wavelengths for TERMS
freq <- seq(240,400,by=1.6)

nm <- round(299792.458 / freq, 8)

cat(length(nm),"\n", file='wavelength_tio2')
write.table( nm, file='wavelength_tio2', append = TRUE,
             row.names = FALSE, quote = FALSE,col.names = FALSE)


## rewrite h5 file to use nanometers

library(rhdf5)
unlink('cylinder_nm.h5')
fs::file_copy('cylinder_tio2.h5', new_path = 'cylinder_nm.h5')
f <- 'cylinder_nm.h5'
fid <- H5Fopen(f)
fid
freq <- h5read(fid,"frequency")
h5delete(fid,'frequency')
# h5write(nm, fid, 'vacuum_wavelength',storage.mode='FLOAT')
h5createDataset(fid, "vacuum_wavelength", length(nm), storage.mode = "double")
h5write(nm, fid, 'vacuum_wavelength')
H5close()


unlink('spheroid_nm.h5')
fs::file_copy('tm_spheroid_au_h2o_radius_20.00_40.00-all.tmat.h5', new_path = 'spheroid_nm.h5')
f <- 'spheroid_nm.h5'
fid <- H5Fopen(f)
fid
nm <- h5read(fid,"vacuum_wavelength")
h5delete(fid,'vacuum_wavelength')
h5createDataset(fid, "vacuum_wavelength", length(nm), storage.mode = "double")
h5write(round(nm*1e9), fid, 'vacuum_wavelength')
did <- H5Dopen(fid, 'vacuum_wavelength')
h5writeAttribute('nm', did, 'unit')
H5Dclose(did)
H5close()


