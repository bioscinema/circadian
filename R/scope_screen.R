#' Spectrum Analysis with Optional Detrending and Power Correction
#'
#' @description
#' Compute the power spectrum of a time series with optional detrending
#' (via differencing or moving-average high-pass filtering), flexible
#' spectral estimation methods (FFT or Lomb–Scargle), and energy correction
#' for filter and window effects.
#'
#' @param t Numeric vector of sampling times.
#' @param y Numeric vector of observations (same length as `t`).
#' @param detrend Character. Detrending mode:
#'   - `"diff"`: apply \eqn{d}-order differencing;
#'   - `"ma"`: apply moving-average (MA) high-pass filtering;
#'   - `"none"`: no detrending.
#' @param d_order Integer. Differencing order when `detrend = "diff"`.
#' @param k_ma Integer. Window size for moving-average filter when `detrend = "ma"`.
#'   If `NULL`, the function automatically chooses an appropriate window length.
#' @param method Character. Spectrum estimation method:
#'   - `"auto"`: automatically choose `"fft"` for uniform sampling,
#'     otherwise `"lomb"`;
#'   - `"fft"`: Fast Fourier Transform;
#'   - `"lomb"`: Lomb–Scargle periodogram (requires the **lomb** package).
#' @param window Character. Window function for FFT. Either `"hann"` or `"rect"`.
#' @param correct Logical. Whether to apply power correction by dividing by
#'   \eqn{|H(\omega)|^2} (detrend filter gain) and coherent window gain.
#' @param interp_dt Numeric. Time-step used for interpolation when using FFT
#'   with nonuniform sampling. Defaults to the median time interval.
#' @param plot Logical. Whether to return a `ggplot2` spectrum plot.
#'
#' @details
#' The function first removes missing values, then optionally detrends
#' the signal. When `method = "auto"`, it checks sampling uniformity based
#' on the median interval deviation.  
#'
#' For FFT-based spectra, it supports Hann and rectangular windows, and applies
#' coherent gain correction when `correct = TRUE`. For nonuniform sampling,
#' the signal is interpolated to a uniform grid before FFT computation.
#'
#' The Lomb–Scargle option uses the \pkg{lomb} package to estimate the spectrum
#' directly from irregular time points. The correction factor \eqn{|H(\omega)|^2}
#' adjusts the power spectrum to approximate the "true" spectral power after
#' detrending.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{spectrum}{A `data.frame` containing frequencies (`freq`),
#'   raw power (`P_raw`), and corrected power (`P_corr`).}
#'   \item{plot}{A `ggplot2` object (if `plot = TRUE`).}
#'   \item{details}{A list of additional analysis details,
#'   including detrend method, filter parameters, and window corrections.}
#' }
#'
#' @examples
#' # Simulate a noisy sinusoidal signal
#' t <- seq(0, 48, by = 1)
#' y <- sin(2 * pi * t / 24) + rnorm(length(t), sd = 0.3)
#'
#' # Compute spectrum using FFT with Hann window
#' res <- spectrum_with_detrend(t, y, detrend = "diff", method = "fft", window = "hann")
#' res$plot
#'
#' # Use Lomb–Scargle method for nonuniform sampling
#' t2 <- sort(runif(30, 0, 48))
#' y2 <- sin(2 * pi * t2 / 24) + rnorm(length(t2), 0, 0.3)
#' res2 <- spectrum_with_detrend(t2, y2, method = "lomb")
#'
#' @export
spectrum_with_detrend <- function(t, y,
                                  detrend = c("diff","ma","none"),
                                  d_order = 1, k_ma = NULL,
                                  method  = c("auto","fft","lomb"),
                                  window  = c("hann","rect"),
                                  correct = TRUE,
                                  interp_dt = NULL,
                                  plot = TRUE){
  stopifnot(is.numeric(t), is.numeric(y), length(t)==length(y))
  detrend <- match.arg(detrend)
  method  <- match.arg(method)
  window  <- match.arg(window)
  
  # ---- prep & NA guard ----
  ok <- is.finite(t) & is.finite(y)
  t <- as.numeric(t[ok]); y <- as.numeric(y[ok])
  o <- order(t); t <- t[o]; y <- y[o]
  n <- length(y); if (n < 8) stop("Too few points.")
  
  # ---- helper: detrend filters H(omega) ----
  H_mag2 <- function(omega){
    if (detrend == "none") return(rep(1, length(omega)))
    if (detrend == "diff"){
      return((2*sin(omega/2))^(2*d_order))
    }
    if (detrend == "ma"){
      K <- k_ma
      num <- sin(K*omega/2)
      den <- sin(omega/2)
      den[abs(den)<.Machine$double.eps] <- .Machine$double.eps
      g <- (num/(K*den))
      return(1 + g^2 - 2*g*cos(omega*(K-1)/2))
    }
  }
  
  # ---- apply detrend in time domain ----
  y0 <- y - mean(y)
  if (detrend == "diff"){
    for (k in seq_len(d_order)){
      y0 <- diff(y0)
      t  <- t[-1]
    }
  } else if (detrend == "ma"){
    if (is.null(k_ma)) {
      dt_med <- median(diff(t))
      k_ma <- max(3, floor(length(y)/8))
    }
    filt <- rep(1/k_ma, k_ma)
    padL <- floor(k_ma/2); padR <- k_ma-1-padL
    pad  <- function(v) c(rep(v[1], padL), v, rep(v[length(v)], padR))
    yy   <- pad(y0)
    ma   <- as.numeric(stats::filter(yy, filt, sides = 2))
    ma   <- ma[(padL+1):(padL+length(y0))]
    y0   <- y0 - ma
  }
  
  # ---- choose method ----
  dt <- diff(t); dt_med <- median(dt)
  uniform <- max(abs(dt - dt_med)) <= 1e-8*max(1, dt_med)
  if (method == "auto") method <- if (uniform) "fft" else "lomb"
  
  # ---- window (FFT only) ----
  w <- NULL; Cw <- 1
  if (method == "fft"){
    M <- length(y0)
    if (window == "hann"){
      w <- 0.5 - 0.5*cos(2*pi*(0:(M-1))/(M-1))
    } else { w <- rep(1, M) }
    Cw <- (sum(w)^2) / sum(w^2)
  }
  
  # ---- spectrum calculation ----
  spec_df <- NULL
  details <- list()
  if (method == "fft"){
    if (!uniform){
      if (is.null(interp_dt) || !is.finite(interp_dt) || interp_dt <= 0)
        interp_dt <- dt_med
      grid <- seq(min(t), max(t), by = interp_dt)
      y0   <- stats::approx(x=t, y=y0, xout=grid, rule=2)$y
      t    <- grid
    }
    M <- length(y0)
    y_win <- y0 * w
    F <- stats::fft(y_win)
    K <- floor(M/2)
    fs <- 1/median(diff(t))
    freq <- (1:(2*K)) * fs / (2*M)
    P_all <- Re(F*Conj(F)) / (M^2)
    Ppos  <- P_all[2:(K+1)]
    omega <- 2*pi*freq/fs
    H2    <- H_mag2(omega)
    
    if (correct){
      P_corr <- Ppos / pmax(H2, .Machine$double.eps) / Cw
    } else {
      P_corr <- Ppos
    }
    
    spec_df <- data.frame(freq = freq, P_raw = Ppos, P_corr = P_corr)
    details$fs <- fs; details$M <- M; details$window <- window
    details$coherent_gain <- Cw; details$H2 <- H2
    
  } else { # Lomb–Scargle
    if (!requireNamespace("lomb", quietly = TRUE))
      stop('method="lomb" requires package lomb')
    Tspan <- diff(range(t))
    dt_min <- min(diff(t))
    f_lo <- 1/72
    f_hi <- 1/(2*dt_min)
    ls <- lomb::lsp(data.frame(t=t, y=y0), from=f_lo, to=f_hi,
                    type="frequency", ofac=1, plot=FALSE)
    freq <- as.numeric(ls$scanned)
    Ppos <- as.numeric(ls$power)
    dt_ref <- median(diff(t))
    omega  <- 2*pi*freq*dt_ref
    H2     <- H_mag2(omega)
    
    if (correct){
      P_corr <- Ppos / pmax(H2, .Machine$double.eps)
    } else {
      P_corr <- Ppos
    }
    spec_df <- data.frame(freq=freq, P_raw=Ppos, P_corr=P_corr)
    details$dt_ref <- dt_ref; details$H2 <- H2
  }
  
  # ---- return results ----
  if(plot){
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("Please install ggplot2")
    library(ggplot2)
    plt <- ggplot(spec_df, aes(x=freq, y=P_corr)) +
      geom_line(linewidth=0.7) +
      labs(x = "Frequency (1 / time unit)",
           y = if (correct) "Corrected power (≈ true power)"
           else "Power (arbitrary units)",
           title = sprintf("Spectrum (%s; detrend=%s)",
                           toupper(method), detrend)) +
      theme_minimal(base_size = 13)
    
    list(spectrum = spec_df, plot = plt,
         details = c(details, list(detrend=detrend, d_order=d_order, k_ma=k_ma,
                                   method=method, correct=correct)))
  } else {
    list(spectrum = spec_df, 
         details = c(details, list(detrend=detrend, d_order=d_order, k_ma=k_ma,
                                   method=method, correct=correct)))
  }
}
