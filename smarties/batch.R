library(glue)
library(dplyr)
library(purrr)

# read in the template
template <- glue_collapse(readLines('smarties/_template.m'), sep = "\n")

parameters <- expand.grid(a = seq(10, 100, by=5), c = seq(10, 100, by=5),
                          material = c("Au", "Ag", "Si"), 
                          medium = c("vacuum", "water")) |> 
  filter(a != c) |> # use Mie for this
  mutate(shape = ifelse(a > c, "oblate", "prolate"),
         n = ifelse(medium == "water", 1.33, 1.0))

nrow(parameters)

