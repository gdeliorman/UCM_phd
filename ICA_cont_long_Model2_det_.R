ICA_contcont_long_ri_det <- function(p = numeric(),
                                     AR1 = TRUE,
                                     AR1_rho = numeric(),
                                     VT0S0,
                                     VT1S1,
                                     VT0T0,
                                     VT1T1,
                                     VS0S0,
                                     VS1S1,
                                     VT0T1 = seq(-1, 1, by = 0.1),
                                     VT0S1 = seq(-1, 1, by = 0.1),
                                     VT1S0 = seq(-1, 1, by = 0.1),
                                     VS0S1 = seq(-1, 1, by = 0.1),
                                     DT0S0,
                                     DT1S1,
                                     DT0T0,
                                     DT1T1,
                                     DS0S0,
                                     DS1S1,
                                     DT0T1 = seq(-1, 1, by = 0.1),
                                     DT0S1 = seq(-1, 1, by = 0.1),
                                     DT1S0 = seq(-1, 1, by = 0.1),
                                     DS0S1 = seq(-1, 1, by = 0.1)) {
  
  validate_scalar <- function(x, name, positive = FALSE) {
    if (length(x) != 1 || !is.finite(x)) {
      stop("`", name, "` must be a single finite numeric value.")
    }
    
    if (positive && x <= 0) {
      stop("`", name, "` must be strictly positive.")
    }
  }
  
  validate_correlation_vector <- function(x, name) {
    if (length(x) < 1 || any(!is.finite(x))) {
      stop("`", name, "` must contain at least one finite value.")
    }
    
    if (any(x < -1 | x > 1)) {
      stop("`", name, "` must contain values between -1 and 1.")
    }
  }
  
  validate_covariance_scalar <- function(cov_xy, var_x, var_y, name) {
    upper_bound <- sqrt(var_x * var_y)
    
    if (!is.finite(cov_xy) || abs(cov_xy) > upper_bound) {
      stop("`", name, "` must satisfy |cov| <= sqrt(var_x * var_y).")
    }
  }
  
  build_covariance_matrix <- function(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,
                                      T0T0, T1T1, S0S0, S1S1) {
    Sigma_c <- diag(c(T0T0, T1T1, S0S0, S1S1))
    
    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * sqrt(T0T0 * T1T1)
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * sqrt(T0T0 * S1S1)
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * sqrt(T1T1 * S0S0)
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * sqrt(S0S0 * S1S1)
    
    Sigma_c
  }
  
  is_positive_definite <- function(mat) {
    chol_result <- try(chol(mat), silent = TRUE)
    !inherits(chol_result, "try-error")
  }
  
  summarize_delta_components <- function(valid_candidates, T0T0, T1T1, S0S0, S1S1,
                                         T0S0_cov, T1S1_cov) {
    if (nrow(valid_candidates) == 0) {
      return(data.frame(
        c11 = numeric(0),
        c12 = numeric(0),
        c22 = numeric(0)
      ))
    }
    
    cov_t0t1 <- valid_candidates$T0T1 * sqrt(T0T0 * T1T1)
    cov_t0s1 <- valid_candidates$T0S1 * sqrt(T0T0 * S1S1)
    cov_t1s0 <- valid_candidates$T1S0 * sqrt(T1T1 * S0S0)
    cov_s0s1 <- valid_candidates$S0S1 * sqrt(S0S0 * S1S1)
    
    c11 <- T0T0 + T1T1 - (2 * cov_t0t1)
    c12 <- T0S0_cov + T1S1_cov - cov_t0s1 - cov_t1s0
    c22 <- S0S0 + S1S1 - (2 * cov_s0s1)
    det_2x2 <- (c11 * c22) - (c12^2)
    
    keep <- is.finite(c11) & is.finite(c12) & is.finite(c22) &
      c11 > 0 & c22 > 0 & det_2x2 > 0
    
    if (!any(keep)) {
      return(data.frame(
        c11 = numeric(0),
        c12 = numeric(0),
        c22 = numeric(0)
      ))
    }
    
    data.frame(
      c11 = c11[keep],
      c12 = c12[keep],
      c22 = c22[keep]
    )
  }
  
  enumerate_valid_matrices <- function(T0S0_cov, T1S1_cov, T0T0, T1T1, S0S0, S1S1,
                                       T0T1, T0S1, T1S0, S0S1) {
    combins <- expand.grid(
      T0T1 = T0T1,
      T0S0 = T0S0_cov,
      T0S1 = T0S1,
      T1S0 = T1S0,
      T1S1 = T1S1_cov,
      S0S1 = S0S1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (nrow(combins) == 0) {
      return(list(total = 0, valid = combins[0, , drop = FALSE]))
    }
    
    keep <- logical(nrow(combins))
    
    for (i in seq_len(nrow(combins))) {
      Sigma_c <- build_covariance_matrix(
        T0T1 = combins$T0T1[i],
        T0S0 = combins$T0S0[i],
        T0S1 = combins$T0S1[i],
        T1S0 = combins$T1S0[i],
        T1S1 = combins$T1S1[i],
        S0S1 = combins$S0S1[i],
        T0T0 = T0T0,
        T1T1 = T1T1,
        S0S0 = S0S0,
        S1S1 = S1S1
      )
      
      if (is_positive_definite(Sigma_c)) {
        keep[i] <- TRUE
      }
    }
    
    list(total = nrow(combins), valid = combins[keep, , drop = FALSE])
  }
  
  build_ar1_matrix <- function(n, rho) {
    rho ^ abs(outer(seq_len(n), seq_len(n), "-"))
  }
  
  logdet_positive_definite <- function(mat) {
    chol_mat <- chol(mat)
    2 * sum(log(diag(chol_mat)))
  }
  
  p <- as.numeric(p)
  AR1_rho <- as.numeric(AR1_rho)
  
  if (!isTRUE(AR1)) {
    stop("Currently only the AR(1) structure is implemented, so `AR1` must be TRUE.")
  }
  
  if (length(p) != 1 || !is.finite(p) || p < 1 || p != round(p)) {
    stop("`p` must be a single positive integer.")
  }
  
  if (length(AR1_rho) != 1 || !is.finite(AR1_rho) || AR1_rho <= -1 || AR1_rho >= 1) {
    stop("`AR1_rho` must be a single finite value in (-1, 1).")
  }
  
  scalar_values <- list(
    VT0S0 = VT0S0,
    VT1S1 = VT1S1,
    VT0T0 = VT0T0,
    VT1T1 = VT1T1,
    VS0S0 = VS0S0,
    VS1S1 = VS1S1,
    DT0S0 = DT0S0,
    DT1S1 = DT1S1,
    DT0T0 = DT0T0,
    DT1T1 = DT1T1,
    DS0S0 = DS0S0,
    DS1S1 = DS1S1
  )
  
  variance_names <- c("VT0T0", "VT1T1", "VS0S0", "VS1S1", "DT0T0", "DT1T1", "DS0S0", "DS1S1")
  
  for (name in names(scalar_values)) {
    validate_scalar(scalar_values[[name]], name, positive = name %in% variance_names)
  }
  
  validate_covariance_scalar(VT0S0, VT0T0, VS0S0, "VT0S0")
  validate_covariance_scalar(VT1S1, VT1T1, VS1S1, "VT1S1")
  validate_covariance_scalar(DT0S0, DT0T0, DS0S0, "DT0S0")
  validate_covariance_scalar(DT1S1, DT1T1, DS1S1, "DT1S1")
  
  correlation_vectors <- list(
    VT0T1 = VT0T1,
    VT0S1 = VT0S1,
    VT1S0 = VT1S0,
    VS0S1 = VS0S1,
    DT0T1 = DT0T1,
    DT0S1 = DT0S1,
    DT1S0 = DT1S0,
    DS0S1 = DS0S1
  )
  
  for (name in names(correlation_vectors)) {
    validate_correlation_vector(correlation_vectors[[name]], name)
  }
  
  v_candidates <- enumerate_valid_matrices(
    T0S0_cov = VT0S0,
    T1S1_cov = VT1S1,
    T0T0 = VT0T0,
    T1T1 = VT1T1,
    S0S0 = VS0S0,
    S1S1 = VS1S1,
    T0T1 = VT0T1,
    T0S1 = VT0S1,
    T1S0 = VT1S0,
    S0S1 = VS0S1
  )
  
  d_candidates <- enumerate_valid_matrices(
    T0S0_cov = DT0S0,
    T1S1_cov = DT1S1,
    T0T0 = DT0T0,
    T1T1 = DT1T1,
    S0S0 = DS0S0,
    S1S1 = DS1S1,
    T0T1 = DT0T1,
    T0S1 = DT0S1,
    T1S0 = DT1S0,
    S0S1 = DS0S1
  )
  
  num_pos_def_v <- nrow(v_candidates$valid)
  num_pos_def_d <- nrow(d_candidates$valid)
  
  v_components <- summarize_delta_components(
    valid_candidates = v_candidates$valid,
    T0T0 = VT0T0,
    T1T1 = VT1T1,
    S0S0 = VS0S0,
    S1S1 = VS1S1,
    T0S0_cov = VT0S0,
    T1S1_cov = VT1S1
  )
  
  d_components <- summarize_delta_components(
    valid_candidates = d_candidates$valid,
    T0T0 = DT0T0,
    T1T1 = DT1T1,
    S0S0 = DS0S0,
    S1S1 = DS1S1,
    T0S0_cov = DT0S0,
    T1S1_cov = DT1S1
  )
  
  max_pairs <- nrow(v_components) * nrow(d_components)
  r2_lambda <- numeric(max_pairs)
  out_index <- 0
  
  if (max_pairs > 0) {
    r_mat <- build_ar1_matrix(p, AR1_rho)
    one_mat <- matrix(1, nrow = p, ncol = p)
    tt_idx <- seq_len(p)
    ss_idx <- (p + 1):(2 * p)
    
    v_blocks <- vector("list", nrow(v_components))
    for (i in seq_len(nrow(v_components))) {
      sigma_u <- matrix(
        c(v_components$c11[i], v_components$c12[i],
          v_components$c12[i], v_components$c22[i]),
        nrow = 2,
        ncol = 2
      )
      v_blocks[[i]] <- kronecker(sigma_u, r_mat)
    }
    
    d_blocks <- vector("list", nrow(d_components))
    for (j in seq_len(nrow(d_components))) {
      sigma_d <- matrix(
        c(d_components$c11[j], d_components$c12[j],
          d_components$c12[j], d_components$c22[j]),
        nrow = 2,
        ncol = 2
      )
      d_blocks[[j]] <- kronecker(sigma_d, one_mat)
    }
    
    for (i in seq_len(nrow(v_components))) {
      for (j in seq_len(nrow(d_components))) {
        sigma_delta <- v_blocks[[i]] + d_blocks[[j]]
        
        if (!is_positive_definite(sigma_delta)) {
          next
        }
        
        sigma_tt <- sigma_delta[tt_idx, tt_idx, drop = FALSE]
        sigma_ss <- sigma_delta[ss_idx, ss_idx, drop = FALSE]
        
        if (!is_positive_definite(sigma_tt) || !is_positive_definite(sigma_ss)) {
          next
        }
        
        log_ratio <- logdet_positive_definite(sigma_delta) -
          logdet_positive_definite(sigma_tt) -
          logdet_positive_definite(sigma_ss)
        
        pair_value <- 1 - exp(log_ratio)
        
        if (!is.finite(pair_value)) {
          next
        }
        
        out_index <- out_index + 1
        r2_lambda[out_index] <- pair_value
      }
    }
  }
  
  if (out_index == 0) {
    r2_lambda <- numeric(0)
  } else {
    r2_lambda <- r2_lambda[seq_len(out_index)]
  }
  
  fit <- list(
    Num.Pos.Def.V = num_pos_def_v,
    Num.Pos.Def.D = num_pos_def_d,
    Num.Pos.Def.Pairs = length(r2_lambda),
    R2_Lambda = r2_lambda,
    Call = match.call()
  )
  
  class(fit) <- "ICA_contcont_long_ri_det"
  fit
}



ex11<- ICA_contcont_long_ri_det(
  p = 6,
  AR1 = TRUE,
  AR1_rho = 0.6223,
  VT0S0 = 0.03401,
  VT1S1 = 0.03640,
  VT0T0 = 0.03804,
  VT1T1 = 0.04099,
  VS0S0 = 0.03485,
  VS1S1 = 0.03653,
  VT0T1 = seq(-1, 1, by = 0.1),
  VT0S1 = seq(-1, 1, by = 0.1),
  VT1S0 = seq(-1, 1, by = 0.1),
  VS0S1 = seq(-1, 1, by = 0.1),
  DT0S0 = 0.03497,
  DT1S1 = 0.03412,
  DT0T0 = 0.04034,
  DT1T1 = 0.03470,
  DS0S0 = 0.03492,
  DS1S1 = 0.03800,
  DT0T1 = seq(-1, 1, by = 0.1),
  DT0S1 = seq(-1, 1, by = 0.1),
  DT1S0 = seq(-1, 1, by = 0.1),
  DS0S1 = seq(-1, 1, by = 0.1) )


ex11$Num.Pos.Def.V
ex11$Num.Pos.Def.D
ex11$Num.Pos.Def.Pairs
min(ex11$R2_Lambda)
max(ex11$R2_Lambda)
sd(ex11$R2_Lambda)
hist(ex11$R2_Lambda)
median(ex11$R2_Lambda)
mean(ex11$R2_Lambda)
hist(ex11$R2_Lambda, xlim=c(0,1))

ex2<- ICA_contcont_long_ri_det(
  p = 6,
  AR1 = TRUE,
  AR1_rho = 0.6223,
  VT0S0 = 0.03401,
  VT1S1 = 0.03640,
  VT0T0 = 0.03804,
  VT1T1 = 0.04099,
  VS0S0 = 0.03485,
  VS1S1 = 0.03653,
  VT0T1 = 0,
  VT0S1 = 0,
  VT1S0 = 0,
  VS0S1 = 0,
  DT0S0 = 0.03497,
  DT1S1 = 0.03412,
  DT0T0 = 0.04034,
  DT1T1 = 0.03470,
  DS0S0 = 0.03492,
  DS1S1 = 0.03800,
  DT0T1 = 0,
  DT0S1 = 0,
  DT1S0 = 0,
  DS0S1 = 0 )

ex2$R2_Lambda

