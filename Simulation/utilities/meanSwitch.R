# fn.mu <- fn.mu.al
# fn.mu <- fn.mu.anl
if (smooth_mean){
  cat("Mean function is infintely differentiable\n")
  fn.mu <- fn.mu.int
} else{
  cat("Mean function is upto twice differentiable\n")
  fn.mu <- fn.mu.c2
}


# delta <- 0
