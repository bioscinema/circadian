simulate_mix_periods <- function(
    m = 1000,
    # ---- time settings ----
    times = NULL,                        
    time_mode = c("uniform", "nonuniform"),
    # uniform: start + by * (0:(length.out-1))
    uniform_start = 0, uniform_by = 1, uniform_length = 24,
    nonuniform_start = 0, nonuniform_end = 24, nonuniform_n = 24,
    # ---- signal settings ----
    periods = c(24, 17, 11, 35, NA),     
    frac = rep(1/5, 5),                  
    A_range = c(0.8, 2.0),
    c_mean_sd = c(10, 2),
    sigma = 1.0,
    seed = 2025
){
  stopifnot(length(periods) == length(frac))
  stopifnot(abs(sum(frac) - 1) < 1e-8)
  
  set.seed(seed)
  
  # ----- build time vector -----
  if (is.null(times)) {
    time_mode <- match.arg(time_mode)
    if (time_mode == "uniform") {
      times <- uniform_start + uniform_by * (0:(uniform_length - 1))
    } else {
      times <- sort(runif(nonuniform_n, min = nonuniform_start, max = nonuniform_end))
    }
  } else {
    if (!is.numeric(times)) stop("'times' must be numeric (in hours).")
    times <- as.numeric(times)
  }
  n <- length(times)
  
  # ----- assign gene types by periods/frac -----
  k <- as.vector(rmultinom(1, size = m, prob = frac))
  lab_each <- ifelse(is.na(periods), "nor", paste0("p", periods))
  type <- rep(lab_each, times = k)
  type <- sample(type)  
  
  # ----- parameters per gene -----
  A   <- ifelse(type == "nor", 0, runif(m, A_range[1], A_range[2]))
  c0  <- rnorm(m, mean = c_mean_sd[1], sd = c_mean_sd[2])
  phi <- runif(m, -pi, pi)
  
  period_map <- setNames(periods, lab_each)
  period_true <- unname(period_map[type])
  
  # ----- generate expression matrix (genes x samples) -----
  Y <- t(vapply(seq_len(m), function(j){
    if (is.na(period_true[j])) {
      mu <- c0[j] + 0*times
    } else {
      w  <- 2*pi/period_true[j]
      mu <- c0[j] + A[j] * sin(w * times + phi[j])
    }
    mu + rnorm(n, 0, sigma)
  }, numeric(n)))
  
  rownames(Y) <- paste0("gene_", seq_len(m))
  colnames(Y) <- paste0("T", seq_len(n))
  
  list(
    expr = Y,
    time = times,
    type = type,                 
    phi = phi,
    period_true = period_true,  
    meta = list(
      m = m,
      time_mode = if (is.null(match.call()$times)) match.arg(time_mode) else "provided",
      periods = periods,
      frac = frac,
      A_range = A_range,
      c_mean_sd = c_mean_sd,
      sigma = sigma,
      seed = seed
    )
  )
}
