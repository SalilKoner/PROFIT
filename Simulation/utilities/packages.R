# R packages

if (!require(pacman, quietly = TRUE)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse, nlme, fda, face, fdapace, refund, RLRsim, mgcv, rlang, doSNOW)


mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
}