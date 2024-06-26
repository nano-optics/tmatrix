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
read_treams <- function(f){
  tmat <- rhdf5::h5read(f, 'tmatrix', compoundAsDataFrame = FALSE)
  modes <- rhdf5::h5read(f, 'modes', compoundAsDataFrame = FALSE)
  params <- rhdf5::h5read(f, 'computation')$method_parameters
  
  nonzero <- is.finite(log10(Mod(tmat)))
  ind <- which(nonzero, arr.ind = TRUE)
  
  n <- nrow(tmat)
  modes$s <- ifelse(modes$polarization == 'electric', 1, 2)
  l <- matrix(modes$l, nrow=n, ncol=n, byrow=FALSE); lp=t(l)
  m <- matrix(modes$m, nrow=n, ncol=n, byrow=FALSE); mp=t(m)
  s <- matrix(modes$s, nrow=n, ncol=n, byrow=FALSE); sp=t(s)
  
  treams <- data.frame(s = s[ind], sp = sp[ind], 
                       l = l[ind], lp = lp[ind],  
                       m = m[ind], mp = mp[ind], 
                       Tr = Re(tmat[ind]), Ti = Im(tmat[ind])) |> 
    arrange(s,sp,l,lp,m,mp)
  
  treams$mod <- Mod(treams$Tr + 1i*treams$Ti)
  
  treams$p <- p_index(treams$l, treams$m)
  treams$pp <- p_index(treams$lp, treams$mp)
  treams$q <- q_index(treams$p, treams$s, max(treams$p))
  treams$qp <- q_index(treams$pp, treams$sp, max(treams$pp))
  
  return(cbind(treams, params))
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



