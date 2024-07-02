library(here)
setwd("~/Documents/nano-optics/tmatrix/terms")

## ----load----
suppressPackageStartupMessages(require(terms))
theme_set(theme_grey())

## ----run----
# system("terms input_test > log")
# 
# 
# 
# system("terms input_smarties > log1")
# system("terms input_h5 > log2")

## ----read----
# grab the results
# xs <- terms::consolidate_xsec('cross_sections_tio2.h5')
xs <- terms::consolidate_xsec('cross_sections.h5')

# plot OA data
## ----oa----

glimpse(xs$mCOA)
mCOAt <- subset(xs$mCOA, variable == 'total')

p1 <- ggplot(mCOAt, aes(wavelength, average,colour=crosstype)) +
  facet_wrap(~crosstype)+
  geom_line(data=mCOAt, lty=1) +
  guides(colour='none')+
  scale_color_brewer(palette = 'Set1') +
  labs(x = expression("wavelength /nm"), 
       y = expression("cross-section /"*nm^2),
       colour = expression(type)) +
  ggtitle("Orientation-averaged cross-sections")

p1

## ----cleanup----
# unlink("*.dat")

