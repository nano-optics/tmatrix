library(rhdf5) # note: dev version to support complex
library(purrr) # mapping functions
library(lobstr) # viewing functions

# custom_inline_hook <- function (x) 
# {
#   if (is.numeric(x)) 
#     x = round_digits(x)
#   paste(as.character(x), collapse = ", ")
# }


check_h5 <- function(f){
  d <- h5dump(f, read.attributes = T, native=TRUE,drop=TRUE )
  
  a <- list(root = rhdf5::h5readAttributes(f,'/'),
            vacuum_wavelength = rhdf5::h5readAttributes(f,'vacuum_wavelength'),
            computation = rhdf5::h5readAttributes(f,'computation'),
            embedding = rhdf5::h5readAttributes(f,'embedding'),
            material = rhdf5::h5readAttributes(f,'scatterer/material'),
            geometry = rhdf5::h5readAttributes(f,'scatterer/geometry'))
  
  r <- list(data = d, attributes = a)
  # 
  # tmp <- tempfile()
  # zz <- file(tmp, open = "wt")
  # sink(zz)
  # print()
  # sink()
  # close(zz)
  # txt <- paste(readLines(tmp))
  # unlink(tmp)
  # txt
  r
}


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

## ---- tmatrix_combinedindex
p_index <- function(l, m){
  p <- l*(l+1) + m
  return(p)
}

q_index <- function(p, s, pmax){
  
  q <- (s-1)*pmax + p
  return(q)
  
}

tmatrix_combinedindex <- function(d, lmax=max(d$l)){
  
  pmax <- lmax*(lmax + 1) + lmax
  mutate(d, p = p_index(l, m), 
         pp = p_index(lp, mp),
         q = q_index(p,  s,  pmax),
         qp = q_index(pp,  sp,  pmax))
  
}


## ---- tmatrix_breaks
tmatrix_breaks <- function(lmax){
  
  l <- seq.int(lmax)
  qmax <- 2*(lmax*(lmax + 1)+lmax)
  list(breaks = cumsum(rep((2*l+1), 2)),
  labels = rep(cumsum((2*l+1)), 2),
  minor_breaks = seq.int(qmax))
}

## ---- read_treams

read_tmat <- function(f){
  h5closeAll()
  
  tmatrix <- rhdf5::h5read(f, 'tmatrix', compoundAsDataFrame = FALSE, native = TRUE)
  modes <- rhdf5::h5read(f, 'modes', compoundAsDataFrame = FALSE, native = TRUE)
  wavelengths <- rhdf5::h5read(f, 'vacuum_wavelength', compoundAsDataFrame = FALSE, native = TRUE)
  # ind633 <- which(abs(wavelengths - 6.3e-7 ) < 0.01e-9)
  tmat_dims <- dim(tmatrix)
  message("tmatrix dimensions: ",paste(tmat_dims, collapse="x"))
  # tmatrixcan be either a matrix (single wavelength), or an array with 1 dimension wavelengths
  # so we need to handle both cases separately (it's awkward to slice by first index)
  if(length(tmat_dims) == 2){
    n <- tmat_dims[1]
    
  } else if(length(tmat_dims) == 3){
    stopifnot(tmat_dims[1] == length(wavelengths))
    n <- tmat_dims[2] # size of matrix
    
  } else {
    error("dimensions of t-matrix should be 2 or 3")
  }
  
  modes$s <- ifelse(modes$polarization == 'magnetic', 1, 2)
  l <- matrix(modes$l, nrow=n, ncol=n, byrow=FALSE); lp=t(l)
  m <- matrix(modes$m, nrow=n, ncol=n, byrow=FALSE); mp=t(m)
  s <- matrix(modes$s, nrow=n, ncol=n, byrow=FALSE); sp=t(s)
  
  # process a single wavelength
  single_wavelength <- function(tmat){
    
    nonzero <- is.finite(log10(Mod(tmat)))
    ind <- which(nonzero, arr.ind = TRUE)
    
    long_tmat<- data.frame(s = s[ind], sp = sp[ind], 
                              l = l[ind], lp = lp[ind],  
                              m = m[ind], mp = mp[ind], 
                              value = tmat[ind], 
                              Tr = Re(tmat[ind]), Ti = Im(tmat[ind])) |> 
      arrange(s,sp,l,lp,m,mp)
    
    long_tmat$mod <- Mod(long_tmat$Tr + 1i*long_tmat$Ti)
    
    long_tmat$p <- p_index(long_tmat$l, long_tmat$m)
    long_tmat$pp <- p_index(long_tmat$lp, long_tmat$mp)
    long_tmat$q <- q_index(long_tmat$p, long_tmat$s, max(long_tmat$p))
    long_tmat$qp <- q_index(long_tmat$pp, long_tmat$sp, max(long_tmat$pp))
    
    return(long_tmat)
  }
  
  if(length(tmat_dims) == 2){
    result <- single_wavelength(tmatrix)
  } else {
    # this will slice by first index, and return a list with Nl entries
    result <- apply(tmatrix, 1, single_wavelength, simplify = FALSE)
  } 
  return(result)
}


## ---- tmatrix_wide
tmatrix_wide <- function(d, lmax = max(d$l)){
  
  N <- 2*(lmax*(lmax+1) + lmax)
  m <- matrix(0.0+0.0i, N, N)
  d$p <- p_index(l = d$l, m = d$m)
  d$pp <- p_index(l = d$lp, m = d$mp)
  pmax <- lmax*(lmax + 1) + lmax
  d$q <- q_index(p = d$p, s = d$s, pmax)
  d$qp <- q_index(p = d$pp, s = d$sp, pmax)
  
  # fill values
  m[cbind(d$q,d$qp)] <- d$Tr + 1i*d$Ti
  
  m
}



## ---- tmatGrob
library(grid)
library(scales)
tmatGrob <- function(m){ #, min_value=min(m),max_value=max(m)
  rc <- dim(m)
  # m <- oob_censor(m, range=c(min_value, max_value), only.finite = FALSE)
  d <- scales::rescale(m, from = range(m, na.rm = TRUE, finite = T), to=c(1,0))
  d[is.infinite(d)] <- NA
  dim(d) <- rc
  g1 <- rectGrob(width=unit(1,'snpc'),height=unit(1,'snpc'),gp=gpar(fill='cornsilk'))
  g2 <- rasterGrob(d, interpolate = FALSE)
  # n <- nrow(d)
  # xx <- seq(-0.5+1/n,0.5-1/n,length=n)
  # grid.segments(x0=unit(0.5,'npc') - unit(xx,'snpc'),
  #               x1=unit(0.5,'npc') - unit(xx,'snpc'),
  #           y0=unit(0,'snpc'), y1= unit(1,'snpc'))
  # grid.segments(y0=unit(0.5,'snpc') - unit(xx,'snpc'),
  #               y1=unit(0.5,'snpc') - unit(xx,'snpc'),
  #               x0=unit(0.25,'npc'),x1= unit(0.75,'npc'))
  grobTree(g1,g2)
}


grid.tmat <- function(...) grid.draw(tmatGrob(...))


prettymat <- function(m, digits=2,...){
  m[] <- format(m, digits=digits)
  ramify::pprint(m, ...)
}



