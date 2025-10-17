#' Simulate Mixed-Period Circadian Gene Expression Data
#'
#' @description
#' Generate synthetic gene expression matrices containing a mixture of rhythmic and non-rhythmic genes
#' with user-defined oscillation periods, amplitudes, and noise. The function allows both uniformly
#' and non-uniformly spaced sampling times, providing flexible simulation settings for benchmarking
#' circadian analysis methods.
#'
#' @param m Integer. Number of genes (rows) to simulate. Default is 1000.
#' @param times Optional numeric vector of sampling times (in hours). If `NULL`, the time points
#'   are generated according to `time_mode`.
#' @param time_mode Character. Either `"uniform"` or `"nonuniform"`.  
#'   - `"uniform"`: equally spaced times specified by `uniform_start`, `uniform_by`, `uniform_length`.  
#'   - `"nonuniform"`: random sampling within [`nonuniform_start`, `nonuniform_end`] with `nonuniform_n` samples.
#' @param uniform_start,uniform_by,uniform_length Numeric. Parameters controlling the evenly spaced time grid
#'   when `time_mode = "uniform"`.
#' @param nonuniform_start,nonuniform_end Numeric. Range for randomly sampled times when
#'   `time_mode = "nonuniform"`.
#' @param nonuniform_n Integer. Number of samples when generating non-uniform times.
#' @param periods Numeric vector. A set of oscillation periods (in hours) assigned to rhythmic genes.
#'   Can include `NA` to represent non-rhythmic genes.
#' @param frac Numeric vector of the same length as `periods`, giving the fraction of genes
#'   assigned to each period type. Must sum to 1.
#' @param A_range Numeric length-2 vector. Range of amplitudes for rhythmic genes. Default is `c(0.8, 2.0)`.
#' @param c_mean_sd Numeric length-2 vector. Mean and standard deviation of the baseline expression level.
#'   Default is `c(10, 2)`.
#' @param sigma Numeric. Standard deviation of Gaussian noise added to each expression value.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @details
#' Each gene is simulated according to the following model:
#' \deqn{
#' Y_{ij} = c_{0i} + A_i \sin\left(\frac{2\pi t_j}{P_i} + \phi_i\right) + \epsilon_{ij},
#' }
#' where \(t_j\) is the sampling time, \(A_i\) is amplitude, \(c_{0i}\) is baseline expression,
#' \(P_i\) is period (possibly `NA` for non-rhythmic genes), and \(\epsilon_{ij}\sim N(0, \sigma^2)\).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{expr}{Matrix of simulated expression values (genes × samples).}
#'   \item{time}{Vector of sampling times.}
#'   \item{type}{Character vector labeling each gene type (`"p24"`, `"p17"`, `"nor"`, etc.).}
#'   \item{phi}{Vector of phase shifts (radians).}
#'   \item{period_true}{Vector of true periods assigned to each gene.}
#'   \item{meta}{List containing all simulation parameters.}
#' }
#'
#' @examples
#' # Example 1: Default settings with uniform time grid
#' sim <- simulate_mix_periods()
#' dim(sim$expr)
#'
#' # Example 2: Non-uniform time sampling
#' sim2 <- simulate_mix_periods(time_mode = "nonuniform", sigma = 0.5)
#' plot(sim2$time, sim2$expr[1, ], type = "l", main = sim2$type[1])
#'
#' @export
sim.mixed_freq <- function(
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
