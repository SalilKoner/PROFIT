# R packages
require(nlme)
require(tidyverse)
require(fda)
require(doSNOW)
require(face)
require(fdapace)
require(refund)
require(RLRsim)
require(mgcv)
require(rlang)


mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
}