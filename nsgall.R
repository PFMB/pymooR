rm(list = ls())
set.seed(2)
library(purrr)
library(magick)
source("./utils.R")

# need to choose conflicting objectives
obj1 <- F

if (obj1) {
  
  # pareto-front has cardinality > 1
  fn <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 # f1
    y[2] <- (x[1] + 2)^2 + (x[2] - 1)^2 # f2
    setNames(y, c("f1","f2")) # solutions need to start with little f
  }
} else {
  
  # pareto-front has cardinality = 1
  fn <- function(x) {
    y <- numeric(2)
    y[1] <- x[1] # f1
    y[2] <- 1 + x[2] - x[1]^2 # f2
    setNames(y, c("f1","f2")) # solutions need to start with little f
  }
  
}

# initialization
x1 <- runif(8, min = -5, max = 5)
#x2 <- runif(8, min = -5, max = 5)
x2 <- sample(-5:5, 8) # integer problem
X <- cbind(x1,x2) # force arguments to start with little x
f <- t(apply(X, 1, fn)) 
P <- cbind(X, f) # population
plot_sol(P)

ranked_P <- rank_cd_sort(P)

for (i in 1:20) {
  
  P <- create_off_springs(P = P, ranked_P = ranked_P, fn = fn)
  
  # stage 1 (first front)
  first_stage <- stage1_front(P[,3:4, drop = F])
  F_1 <- (1:first_stage$M)[as.logical(first_stage$rank)] # front one
  
  # stage 2 (further fronts)
  # all_fronts: the index in the list is the rank
  all_fronts <- stage2_front(F_1, first_stage$S, first_stage$n, first_stage$rank)
  cd_per_front <- unlist(map(all_fronts, ~crowd_dist(.x, P[, 3:4, drop = F])))
  
  selected_idx <- select_best_half(all_fronts, cd_per_front)
  P <- P[selected_idx,]
  
  #P <- P[sample(1:8, 8), ] need to reshuffle?
  
  ranked_P <- rank_cd_sort(P)
  #Sys.sleep(0.3)
  plot_sol(P, i, save = T)
  #plot_sol(P)
  
}

make_gif(2)