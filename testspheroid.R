setwd("~/Documents/nano-optics/tmatrix")
source("_fun.R")
library(ggplot2)
library(dplyr)

d <- read.table('data/tmat_Au20x40_Nmax3.tmat')
names(d) <- c('s','sp','l','lp','m','mp','Tr','Ti')
d <- tmatrix_combinedindex(d)

f <- 'data/tm_spheroid_au_h2o_radius_20.00_40.00-all.tmat.h5'
slices <- read_tmat(f)
# d2 <- single_wavelength(tmat[1,,]) 
d2 <- slices[[ind633]]

f <- 'data/smarties_spectrum.tmat.h5'
slices2 <- read_long_tmat(f)
d3 <- slices2[[1]]
# 
# h5closeAll()
# 
# tmatrix<- rhdf5::h5read(f, 'tmatrix', compoundAsDataFrame = FALSE, native = TRUE)
# modes <- rhdf5::h5read(f, 'modes', compoundAsDataFrame = FALSE, native = TRUE)
# wavelengths <- rhdf5::h5read(f, 'vacuum_wavelength', compoundAsDataFrame = FALSE, native = TRUE)
# wavelengths
# # ind633 <- which(abs(wavelengths - 6.3e-7 ) < 0.01e-9)
# dim(tmat)
# 
# n <- nrow(tmat[1,,])
# modes$s <- ifelse(modes$polarization == 'magnetic', 1, 2)
# l <- matrix(modes$l, nrow=n, ncol=n, byrow=FALSE); lp=t(l)
# m <- matrix(modes$m, nrow=n, ncol=n, byrow=FALSE); mp=t(m)
# s <- matrix(modes$s, nrow=n, ncol=n, byrow=FALSE); sp=t(s)
# 
# single_wavelength <- function(tmat){
#   
#   nonzero <- is.finite(log10(Mod(tmat)))
#   ind <- which(nonzero, arr.ind = TRUE)
#   
#   long_tmatrix<- data.frame(s = s[ind], sp = sp[ind], 
#                        l = l[ind], lp = lp[ind],  
#                        m = m[ind], mp = mp[ind], 
#                        value = tmat[ind], 
#                        Tr = Re(tmat[ind]), Ti = Im(tmat[ind])) |> 
#     arrange(s,sp,l,lp,m,mp)
#   
#   long_tmat$mod <- Mod(long_tmat$Tr + 1i*long_tmat$Ti)
#   
#   long_tmat$p <- p_index(long_tmat$l, long_tmat$m)
#   long_tmat$pp <- p_index(long_tmat$lp, long_tmat$mp)
#   long_tmat$q <- q_index(long_tmat$p, long_tmat$s, max(long_tmat$p))
#   long_tmat$qp <- q_index(long_tmat$pp, long_tmat$sp, max(long_tmat$pp))
#   
#   return(long_tmat)
# }
# 
# slices <- apply(tmat, 1, single_wavelength)
# slices

dm <- rbind(mutate(d, type = "SMARTIES TXT"),
            mutate(d3[,names(d)], type = "SMARTIES H5"),
           mutate(d2[,names(d)], type = "long_tmatrixH5"))

# update the plot with these data, and facet by type

lmax <- max(d$l)
breaks <- tmatrix_breaks(lmax)

p <- ggplot(dm, aes(q, qp, fill= log10(Mod(Tr + 1i*Ti)))) +
  geom_raster() +
  coord_equal() + facet_wrap(~type) +
  scale_fill_viridis_c(option = 'A', direction = -1) +
  annotate('segment',x=0.5,xend=max(breaks$breaks)+0.5,y=max(breaks$breaks)/2+0.5,
           yend=max(breaks$breaks)/2+0.5,colour='white')+
  annotate('segment',y=0.5,yend=max(breaks$breaks)+0.5,x=max(breaks$breaks)/2+0.5,
           xend=max(breaks$breaks)/2+0.5,colour='white')+
  scale_y_reverse(expand=c(0,0), breaks= breaks$breaks+0.5, minor_breaks=breaks$minor_breaks+0.5, labels=breaks$labels) +
  scale_x_continuous(expand=c(0,0), breaks= breaks$breaks+0.5, minor_breaks=breaks$minor_breaks+0.5, labels=breaks$labels) +
  theme_minimal() +
  theme(panel.grid = element_line(colour = 'white'), 
        panel.background = element_rect(fill='grey90',colour='white'),
        panel.border = element_rect(colour='black',fill=NA,linewidth = 0.2),
        axis.text.x = element_text(hjust=1),
        axis.text.y = element_text(vjust=0)) +
  labs(x="p",y="p'",fill=expression(log~"|T|"))

p
