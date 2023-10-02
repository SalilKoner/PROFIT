# R packages

if (!require(pacman, quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, nlme, fda, face, fdapace, refund, RLRsim, mgcv, rlang)


mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
}