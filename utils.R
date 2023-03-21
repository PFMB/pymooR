
cross_sbx <- function(sorted_P, l_b = -5, u_b = 5, eta = 15, r = FALSE, prob_bin = NULL,
                      cross = NULL, rand_alpha = NULL, b_cross = NULL) {
  
  # pymoo implementation operators: crossover/sbx.py
  
  n_P <- nrow(sorted_P)
  
  if(is.null(prob_bin)) {
    # for reproduction
    prob_bin <- rep(0.5, n_P)
    prob_bin[runif(n_P) > 0.5] = 0
  }
  
  if(is.null(cross)) {
    # for reproduction
    cross <- runif(n_P) < 0.5 # cross-over probability
  }
  
  too_close <- abs(sorted_P[,1] - sorted_P[,2]) <= 1e-12
  cross[too_close] <- FALSE
  
  # TODO: No crossover when lower and upper bound are identical
  
  y1 <- apply(sorted_P, 1, min)[cross]
  y2 <- apply(sorted_P, 1, max)[cross]
  
  upper_b <- rep(u_b, n_P)[cross]
  lower_b <- rep(l_b, n_P)[cross]
  
  prob_bin <- prob_bin[cross]
  
  if(is.null(rand_alpha)) {
    # for reproduction
    rand_alpha <- runif(sum(cross))
  }
  
  calc_betaq <- function(beta, eta, rand_alpha) {
    
    alpha <- 2 - beta^(-(eta + 1))
    
    mask <- rand_alpha <= (1/alpha)
    mask_not <- rand_alpha > (1/alpha)
    
    betaq <- numeric(length(rand_alpha))
    
    tmp <- (rand_alpha*alpha)^(1/(eta + 1))
    betaq[mask] <- tmp[mask]
    
    tmp <- (1/(2 - rand_alpha*alpha))^(1/(eta + 1))
    betaq[mask_not] <- tmp[mask_not]
    
    return(betaq)
  }
  
  delta <- y2 - y1
  beta <- 1 + (2 * (y1 - lower_b) / delta)
  betaq <- calc_betaq(beta, eta, rand_alpha)
  c1 <- 0.5 * ((y1 + y2) - betaq * delta)
  
  beta <- 1 + (2 * (upper_b - y2))
  betaq <- calc_betaq(beta, eta, rand_alpha)
  c2 <- 0.5 * ((y1 + y2) + betaq * delta)
  
  if(is.null(b_cross)) {
    # for reproduction
    b_cross <- runif(sum(cross)) < prob_bin
  }
  
  tmp <- c1[b_cross]
  c1[b_cross] <- c2[b_cross]
  c2[b_cross] <- tmp
  
  # FIXME for a general number of columns
  Q <- sorted_P
  Q[cross, 1] <- c1
  Q[cross, 2] <- c2
  
  repair_clamp <- function(x, xl, xu) {
    
    xu <- rep(xu, length(x))
    xl <- rep(xl, length(x))
    
    I <- x < xl
    x[I] <- xl[I]
    
    I <- x > xu
    x[I] <- xu[I]
    x
  }
  
  # FIXME: apply to a general number of variables
  Q[,1] <- repair_clamp(Q[,1], l_b, u_b)
  Q[,2] <- repair_clamp(Q[,2], l_b, u_b)
  
  if (r) return(round(Q))
  
  return(Q)
}

polynnom_mutation <- function(Q, eta = 20, x_lower = NULL, x_upper = NULL, 
                              r = FALSE, mut_poly = NULL, rand_poly = NULL) {
  
  stopifnot(!is.null(x_lower) && !is.null(x_upper))
  
  Q <- c(Q)
  Qp <- Q
  
  if(is.null(mut_poly)) {
    # for reproduction
    mut_poly <- runif(length(Q)) > 0.5
  }
  
  # if(!any(mut_poly)) {
  #   # avoid that all are false
  #   mut_poly[sample(1:length(Q), 1)] <- T
  # }
  
  eta <- rep(eta, sum(mut_poly))
  Q <- Q[mut_poly]
  x_upper <- rep(x_upper, sum(mut_poly))
  x_lower <- rep(x_lower, sum(mut_poly))
  
  delta1 <- (Q - x_lower) / (x_upper - x_lower)
  delta2 <- (x_upper - Q) / (x_upper - x_lower)
  
  mut_pow <- 1 / (eta + 1)
  
  if(is.null(rand_poly)) {
    # for reproduction
    rand_poly <- runif(sum(mut_poly))
  }
  
  mask <- rand_poly <= 0.5
  not_mask <- !mask
  
  deltaq <- numeric(sum(mut_poly))
  
  xy <- 1 - delta1
  val <- 2 * rand_poly + (1 - 2 * rand_poly) * (xy^(eta + 1))
  d <- val^mut_pow - 1
  deltaq[mask] <- d[mask]
  
  xy <- 1 - delta2
  val <- 2 * (1 - rand_poly) + 2 * (rand_poly - 0.5) * (xy^(eta + 1))
  d <- 1 - val^mut_pow
  deltaq[not_mask] <- d[not_mask]
  
  Z <- Q + deltaq * (x_upper - x_lower) # _Y
  
  Z[Z < x_lower] <- x_lower[Z < x_lower] 
  Z[Z > x_upper] <- x_upper[Z > x_upper] 
  
  Qp[mut_poly] <- Z
  
  if (r) return(round(Qp))
  
  # TODO: Set to bounds if outside
  
  Qp
}

check_duplicate <- function(P, P_x_off) {
  
  # write tests to see that the results are permutation invariant
  
  d <- data.frame(rbind(P, P_x_off))
  l <- nrow(d)
  duplicated(d, fromLast = F)[(nrow(P) + 1):nrow(d)]
}

rank_cd_sort <- function(P) {
  
  # P: matrix with arguments in the left columns and solutions in the right columns
  # col name of solutions need to start with "f"
  # col name of arguments need to start with "x"
  
  P_y <- P[,grepl("^f", colnames(P)), drop = T]
  stopifnot(ncol(P) > ncol(P_y) && nrow(P_y) > 1)
  
  # stage 1
  first_stage <- stage1_front(P_y)
  F_1 <- (1:first_stage$M)[as.logical(first_stage$rank)] # front one
  
  # stage 2
  # all_fronts: the index in the list is the rank
  all_fronts <- stage2_front(F_1, first_stage$S, first_stage$n, 
                             first_stage$rank) # further fronts
  
  # crowding distance for each front
  cd_per_front <- unlist(map(all_fronts, ~crowd_dist(.x, P_y)))
  
  # rank before mating
  P[as.numeric(names(cd_per_front)),]
  
}


create_off_springs <- function(P, ranked_P, r, l2, fn, ...) {
  
  P_x <- P[,grepl("^x", colnames(P)), drop = F]
  
  # filter duplicates in the population
  stopifnot(ncol(ranked_P) == ncol(P) && nrow(P) == nrow(ranked_P))
  stopifnot(ncol(P) > ncol(P_x) && nrow(P_x) > 1)
  
  offs <- matrix(NA, dimnames = dimnames(P_x), ncol = 2)
  n_off <- n_remain <- nrow(P)
  n_parents <- 2
  off_spring_tries <- 1
  
  while (n_remain > 0) {
    
    P_x_off <- sapply(seq_len(ncol(P_x)), sbx_mutate, ranked_P = ranked_P, n_off = n_off)
    
    # check if some off-springs were not mutated (are equal to their parents)
    dupli_check <- check_duplicate(P_x, P_x_off)
    P_x_off <- P_x_off[!dupli_check, , drop = F]
    
    if (nrow(P_x_off) > n_remain) {
      
      # generated more off-springs than needed
      P_x_off <- P_x_off[1:n_remain, , drop = F]
      offs <- na.omit(rbind(offs, P_x_off)) # off-spring population
      
    } else {
      
      # generated too few off-springs
      offs <- na.omit(rbind(offs, P_x_off)) # off-spring population
    }
    
    offs <- offs[!duplicated(offs),, drop = F] # should be very unlikely
    
    n_remain <- nrow(P_x) - nrow(offs)
    n_off <- ceiling(n_remain/ n_parents) * n_parents
    cat("Off-spring try:", off_spring_tries, "\n")
    off_spring_tries <- off_spring_tries + 1
    
  }
  
  f <- t(apply(offs, 1, fn))
  colnames(f) <- c("f1","f2") # force solutions to start with little f
  
  P_off <- cbind(offs, f) # off-spring population
  return(rbind(ranked_P, P_off))
  
}

sbx_mutate <- function(col, ranked_P, parent_idx = NULL, mut_sbx = NULL, 
                       skip_Q = FALSE, n_off, ...) {
  
  ## per variable/column, perform sbx and mutate
  # n_off: number offsprings needed
  args <- list(...)
  
  # pick random parents for mating
  if (is.null(parent_idx)) {
    # reproducing
    parent_idx <- sample(1:nrow(ranked_P), replace = F)[1:n_off]
  }
  
  # first and sec parent
  mating_P <- matrix(ranked_P[parent_idx, col], ncol = 2, byrow = T)
  
  # check if you need an integer returned (round)
  is_integer <- all(as.integer(mating_P) == mating_P)
  
  if (skip_Q) {
    Q <- mating_P
  } else {
    Q <- cross_sbx(mating_P, r = is_integer, prob_bin = args$prob_bin,
                   cross = args$cross, rand_alpha = args$rand_alpha,
                   b_cross = args$b_cross) # off-springs
  }

  Xp <- polynnom_mutation(Q, x_lower = -5, x_upper = 5, r = is_integer)
  
  # likelihood for the mutation of an individual
  #mut_sbx <- args$mut_sbx
  
  if(is.null(mut_sbx)) {
    # reproduction
    mut_sbx <- runif(length(Xp)) <= 0.9
  }
  
  Q[mut_sbx] <- Xp[mut_sbx]
  
  return(Q)
}

plot_sol <- function(P, i = 1, save = F) {
  
  # plot the function values for two objectives
  
  if (save) png(paste0("./plots/p_", sprintf("%02d", i),".png"))
  solution_idx <- which(startsWith(colnames(P),"f"))
  f1 <- P[,solution_idx[1]]
  f2 <- P[,solution_idx[2]]
  f1_min_max <- c(min(f1), max(f1))
  f2_min_max <- c(min(f2), max(f2))
  M <- nrow(P)
  plot(f1, f2, xlim = c(-6,6), ylim = c(-30,5))
  text(f1, f2, paste(1:M), pos=2, cex = 0.8)
  legend("topright", paste0("Itearation ", i))
  if (save) dev.off()
}

make_gif <- function(fps) {
  
  imgs <- list.files("plots", full.names = TRUE)
  img_list <- lapply(imgs, image_read)
  
  ## join the images together
  img_joined <- image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = fps)
  
  ## view animated image
  img_animated
  
  ## save to disk
  image_write(image = img_animated,
              path = "plots/evolution.gif")
}

dominates <- function(p, q, P) {
  
  # p: row of the individual in P
  # q: row of the individual in P
  # P: matrix of objective values
  ## minimization ##
  # p dom q iff f_m of p <= f_m of q for all m and at least smaller in one m
  
  all(P[p,] <= P[q,]) && any(P[p,] < P[q,])
  
}

stage1_front <- function(P) {
  
  ## stage 1 of non-dominated sorting
  
  # avoid to get arguments in here as well
  stopifnot(all(grepl("^f", colnames(P))))
  
  M <- nrow(P)
  ranks <- numeric(length = M)
  S_list <- vector("list", length = M)
  n_vec <- numeric(length = M) 
  for (p in 1:M){
    n <- 0
    S <- numeric() # indices that dominate p
    
    for (q in 1:M){
      if(dominates(p, q, P)) {
        S <- c(S, q)
      } else if (dominates(q, p, P)){
        n <- n + 1
      }
    }
    if (n == 0) {
      ranks[p] <- 1
    }
    S_list[[p]] <- S
    n_vec[[p]] <- n
  }
  list(S = S_list, n = n_vec, rank = ranks, M = M)
}

stage2_front <- function(front, S_list, n_vec, ranks) {
  
  ## stage 2
  i <- 1
  fronts <- list(front) # front 1
  
  while (length(fronts[[i]]) > 0) {
    front <- numeric()
    for (p in fronts[[i]]) {
      for (q in S_list[[p]]) {
        n_vec[q] = n_vec[q] - 1
        if (n_vec[q] == 0) {
          ranks[q] = i + 1 
          front <- c(front, q)
        }
      }
    }
    i <- i + 1
    if (length(front) == 0) break
    fronts <- c(fronts, list(front))
  }
  fronts
}

crowd_dist <- function(front, obj) {
  
  # crowding distance is computed per objective/solution
  
  r <- length(front)
  
  if (r <= 2) {
    d <- rep(Inf, r)
    names(d) <- front
    return(d)
  }
  
  d <- numeric(length = r)
  names(d) <- front
  front_obj <- obj[front,, drop = F]
  no_obj <- ncol(front_obj)
  
  for (m in 1:no_obj) {
    
    obj_m <- front_obj[,m]
    names(obj_m) <- front
    
    I <- order(obj_m)
    obj_m_sorted <- obj_m[I]
    norm <- obj_m_sorted[r] - obj_m_sorted[1]
    
    dist <- c(obj_m_sorted, Inf) - c(Inf, obj_m_sorted)
    dist_to_last <- dist[1:(length(dist) - 1)]/norm
    dist_to_next <- dist[2:length(dist)]/norm
    
    J <- order(I)
    
    cd <- dist_to_next[J] + dist_to_last[J]
    cd[is.infinite(cd)] <- Inf # -Inf + Inf gives NaN so make -Inf also Inf 
    d <- d + cd # add across columns/objectives
  }
  
  d/no_obj
}

select_best_half <- function(all_fronts, cd_per_front) {
  
  ## survival stage: choose better solutions for the next generation
  # fill empty slots by rank, if tied with crowding distance
  
  front_len <- sapply(all_fronts, length)
  ranks <- rep(letters[1:length(front_len)], front_len)
  d <- data.frame(ranks = ranks, nm = names(cd_per_front), cd = cd_per_front)
  d <- d[order(d$ranks, d$cd),]
  as.numeric(d[1:(nrow(d) / 2), ]$nm)
}
