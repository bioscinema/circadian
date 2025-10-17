# ================== Spectrum with optional detrend & correction ==================
# y: numeric vector of observations
# t: numeric vector of times (same length)
# detrend: "none" | "diff" | "ma"    (K阶差分 or K点滑动均值高通)
# d_order: 差分阶数（detrend="diff"时用）
# k_ma:    MA窗口长度（detrend="ma"时用）
# method: "auto" | "fft" | "lomb"
#         - auto: 若采样近似等距→FFT；否则→Lomb–Scargle（需 lomb 包）
# window: "rect" | "hann" （只用于 FFT）
# correct: 是否按 |H(ω)|^2 和窗口相干增益做能量校正（用于近似“真功率”）
# interp_dt: 不等距+FFT时的等距插值步长（默认用中位间隔）
# return: list(spectrum=data.frame, plot=ggplot, details=list(...))

spectrum_with_detrend <- function(t, y,
                                  detrend = c("diff","ma","none"),
                                  d_order = 1, k_ma = NULL,
                                  method  = c("auto","fft","lomb"),
                                  window  = c("hann","rect"),
                                  correct = TRUE,
                                  interp_dt = NULL,
                                  plot = T){
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
    # omega: rad/sample for FFT grid; for irregular we’ll map via 2*pi*f*dt_ref
    if (detrend == "none") return(rep(1, length(omega)))
    if (detrend == "diff"){
      # |1 - e^{-iω}|^{2d} = (2 sin(ω/2))^{2d}
      return( (2*sin(omega/2))^(2*d_order) )
    }
    if (detrend == "ma"){
      K <- k_ma
      num <- sin(K*omega/2)
      den <- sin(omega/2)
      den[abs(den)<.Machine$double.eps] <- .Machine$double.eps
      g <- (num/(K*den))
      # |1 - g e^{-iω(K-1)/2}|^2 = 1 + g^2 - 2g cos(ω(K-1)/2)
      return( 1 + g^2 - 2*g*cos(omega*(K-1)/2) )
    }
  }
  
  # ---- apply detrend in time domain (length-safe) ----
  y0 <- y - mean(y)
  if (detrend == "diff"){
    for (k in seq_len(d_order)){
      y0 <- diff(y0)
      t  <- t[-1]
    }
  } else if (detrend == "ma"){
    if (is.null(k_ma)) {
      # 取覆盖低频的窗口，尽量只留高于 ~π/K 的成分
      # 用时间间隔的中位数近似等距长度
      dt_med <- median(diff(t))
      # 目标：K ~ 覆盖 ~ W/4 这样的慢变化；此处用数据长度的 1/8 作为默认
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
  # check uniformity
  dt <- diff(t); dt_med <- median(dt); uniform <- max(abs(dt - dt_med)) <= 1e-8*max(1, dt_med)
  if (method == "auto") method <- if (uniform) "fft" else "lomb"
  
  # ---- window (FFT only) ----
  w <- NULL; Cw <- 1
  if (method == "fft"){
    M <- length(y0)
    if (window == "hann"){
      w <- 0.5 - 0.5*cos(2*pi*(0:(M-1))/(M-1))
    } else { # rect
      w <- rep(1, M)
    }
    # coherent gain Cw = (sum w)^2 / sum w^2
    Cw <- (sum(w)^2) / sum(w^2)
  }
  
  # ---- spectrum calc ----
  spec_df <- NULL
  details <- list()
  if (method == "fft"){
    # FFT expects uniform sampling; if not, interpolate
    if (!uniform){
      if (is.null(interp_dt) || !is.finite(interp_dt) || interp_dt <= 0) interp_dt <- dt_med
      grid <- seq(min(t), max(t), by = interp_dt)
      y0   <- stats::approx(x=t, y=y0, xout=grid, rule=2)$y
      t    <- grid
    }
    M <- length(y0)
    y_win <- y0 * w
    F <- stats::fft(y_win)
    # single-sided spectrum (exclude DC & Nyquist handling)
    K <- floor(M/2)
    # frequencies in Hz (= 1/time unit)
    fs <- 1/median(diff(t))
    freq <- (1:(2*K)) * fs / (2*M)
    # periodogram (Parseval-like): P = |F|^2 / (M^2)
    P_all <- Re(F*Conj(F)) / (M^2)
    Ppos  <- P_all[2:(K+1)]
    
    # map |H(ω)|^2 with ω = 2π f / fs * (radian/sample)
    omega <- 2*pi*freq/fs
    H2    <- H_mag2(omega)
    
    if (correct){
      # divide out filter gain and window coherent gain
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
    # 频率范围：从 1/Tspan 到 1/(2*dt_min)
    Tspan <- diff(range(t))
    dt_min <- min(diff(t))
    f_lo <- 1/72
    f_hi <- 1/(2*dt_min)
    ls <- lomb::lsp(data.frame(t=t, y=y0), from=f_lo, to=f_hi, type="frequency", ofac=1, plot=FALSE)
    freq <- as.numeric(ls$scanned)
    Ppos <- as.numeric(ls$power)
    # 近似把ω用 ω = 2π f * dt_ref (dt_ref=median diff as ref)
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
  
  if(plot){
    # ---- spectrum plot (ggplot) ----
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
  } else{
    list(spectrum = spec_df, 
         details = c(details, list(detrend=detrend, d_order=d_order, k_ma=k_ma,
                                   method=method, correct=correct)))
  }

}

