#' screen_circadian_24h
#' 
#' High-throughput 24-hour circadian screening with options for Cosinor/GLS,
#' neighborhood profile, local period-contrast, circular-shift permutation,
#' and multiple testing correction. Returns per-gene results and (optionally)
#' per-gene plots: fitted overlay and neighborhood profile.
#'
#' @param Y numeric matrix (genes x samples). Rows = genes/features, cols = samples.
#' @param time numeric vector length n (in hours). Sampling time for each sample.
#' @param id optional factor/character length n: replicate/batch/subject ID
#'           (used to scope AR(1) correlation and circular-shift permutation).
#' @param covariates optional data.frame (n x p) of additional covariates.
#' @param method "cosinor" or "gls". cosinor uses OLS (stats::lm), gls uses nlme::gls.
#' @param cor_struct correlation structure for GLS: "none" or "ar1". Ignored for cosinor.
#' @param trend_df integer >=0. If >0, include slow trend via splines::ns(time, df=trend_df).
#' @param periods numeric vector of candidate periods (hours) for neighborhood profile.
#' @param peak_tol numeric. Require |peak_period - 24| <= peak_tol to pass specificity S1.
#' @param contrast "bic" or "aic". Use information criterion to compare 24h vs best alt.
#' @param delta_ic_thresh numeric (default -4). Require IC(24) - IC(alt) <= thresh to pass S2.
#' @param use_permutation logical. If TRUE, compute circular-shift permutation p-values.
#' @param B integer. Number of permutations.
#' @param seed optional integer seed for reproducibility of permutations.
#' @param p_adjust multiple testing method: one of p.adjust methods ("BH","BY","bonferroni",
#'        "holm", etc.), "qvalue" (requires qvalue package), or "none".
#' @param alpha numeric in (0,1]. Significance level for binary call.
#' @param save_plots logical. If TRUE, save per-gene plots to `plot_dir`.
#' @param plot_dir character dir to save plots when save_plots=TRUE.
#' @param return_plots logical. If TRUE, return ggplot objects in the result.
#' @param verbose logical. Print progress.
#'
#' @return a list with elements:
#'   - results: data.frame with per-gene statistics and binary decision `is_24h`.
#'   - overlay_plots: named list of ggplot (if return_plots=TRUE).
#'   - profile_plots: named list of ggplot (if return_plots=TRUE).
#'   - params: a list of key settings used.
#'
#' @importFrom stats lm model.matrix coef na.omit p.adjust quantile predict residuals AIC BIC
#' @importFrom nlme gls corAR1
#' @importFrom splines ns
#' @importFrom utils tail
#' @export
circadian_screen <- function(
    Y,
    time,
    id = NULL,
    covariates = NULL,
    method = c("gls","cosinor"),
    cor_struct = c("ar1","none"),
    trend_df = 3,
    periods = seq(20, 28, by = 0.25),
    peak_tol = 0.5,
    contrast = c("bic","aic"),
    delta_ic_thresh = -4,
    use_permutation = FALSE,
    B = 499,
    seed = NULL,
    p_adjust = c("BH","BY","bonferroni","holm","hochberg","hommel","fdr","qvalue","none"),
    alpha = 0.05,
    save_plots = FALSE,
    plot_dir = NULL,
    return_plots = FALSE,
    verbose = TRUE
){
  # ---- Dependencies ----
  method      <- match.arg(method)
  cor_struct  <- match.arg(cor_struct)
  contrast    <- match.arg(contrast)
  p_adjust    <- match.arg(p_adjust)
  if (method == "gls" && !requireNamespace("nlme", quietly = TRUE)) stop("Package 'nlme' is required for GLS.")
  if (!requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for plotting.")
  if (p_adjust == "qvalue" && !requireNamespace("qvalue", quietly = TRUE)) stop("p_adjust='qvalue' requires the 'qvalue' package.")
  
  # ---- Basic checks ----
  if (is.null(dim(Y))) stop("Y must be a matrix with genes x samples")
  if (!is.numeric(time)) stop("time must be numeric (hours)")
  n <- length(time)
  if (ncol(Y) != n) stop("ncol(Y) must equal length(time)")
  if (!is.null(id) && length(id) != n) stop("length(id) must match length(time)")
  if (!is.null(covariates) && nrow(as.data.frame(covariates)) != n) stop("nrow(covariates) must match length(time)")
  if (trend_df < 0) stop("trend_df must be >= 0")
  
  # ---- Set up design helpers ----
  df <- data.frame(
    time = as.numeric(time),
    id   = if (is.null(id)) factor(rep("1", n)) else factor(id)
  )
  if (!is.null(covariates)) {
    cov_df <- as.data.frame(covariates)
    # sanitize names
    names(cov_df) <- make.names(names(cov_df))
    df <- cbind(df, cov_df)
  }
  
  # Cosinor bases at a given angular frequency
  cos_sin_terms <- function(time, w){
    data.frame(cos_wt = cos(w * time), sin_wt = sin(w * time))
  }
  
  # Build formula strings dynamically
  base_terms <- c()
  if (!is.null(covariates)) base_terms <- c(base_terms, paste0(names(cov_df), collapse = " + "))
  if (trend_df > 0)       base_terms <- c(base_terms, sprintf("splines::ns(time, df=%d)", trend_df))
  base_rhs <- if (length(base_terms)) paste(base_terms, collapse = " + ") else "1"
  
  # Stats utilities
  partial_R2_from_RSS <- function(RSS0, RSS1) max(0, min(1, (RSS0 - RSS1) / RSS0))
  
  # Fit function for a given period
  fit_period <- function(y, P, method, cor_struct, df, base_rhs){
    w <- 2*pi/P
    XY <- cbind(df, cos_sin_terms(df$time, w))
    # Build formulas
    f1 <- as.formula(paste0("y ~ ", base_rhs, " + cos_wt + sin_wt"))
    f0 <- as.formula(paste0("y ~ ", base_rhs))
    
    if (method == "cosinor") {
      d1 <- cbind(XY, y = y)
      d0 <- cbind(df, y = y)
      fit1 <- stats::lm(f1, data = d1)
      fit0 <- stats::lm(f0, data = d0)
      # Extract quantities
      RSS1 <- sum(residuals(fit1)^2, na.rm = TRUE)
      RSS0 <- sum(residuals(fit0)^2, na.rm = TRUE)
      ll1  <- as.numeric(logLik(fit1))
      ll0  <- as.numeric(logLik(fit0))
      ic_fun <- switch(contrast, bic = stats::BIC, aic = stats::AIC)
      ic1 <- ic_fun(fit1)
      ic0 <- ic_fun(fit0)
      coefs <- tryCatch(stats::coef(fit1), error = function(e) rep(NA_real_, 0))
      amp   <- if (all(c("cos_wt","sin_wt") %in% names(coefs))) sqrt(coefs["cos_wt"]^2 + coefs["sin_wt"]^2) else NA_real_
      phi   <- if (all(c("cos_wt","sin_wt") %in% names(coefs))) atan2(-coefs["sin_wt"], coefs["cos_wt"]) else NA_real_
      list(fit1 = fit1, fit0 = fit0, RSS1 = RSS1, RSS0 = RSS0, ll1 = ll1, ll0 = ll0, ic1 = ic1, ic0 = ic0,
           amp = amp, phi = phi, w = w)
    } else { # GLS
      if (cor_struct == "ar1") {
        corr <- nlme::corAR1(form = ~ time | id)
      } else {
        corr <- NULL
      }
      d1 <- cbind(XY, y = y)
      d0 <- cbind(df, y = y)
      fit1 <- nlme::gls(f1, data = d1, correlation = corr, method = "REML")
      fit0 <- nlme::gls(f0, data = d0, correlation = corr, method = "REML")
      # Use logLik and IC
      ll1 <- as.numeric(logLik(fit1))
      ll0 <- as.numeric(logLik(fit0))
      ic_fun <- switch(contrast, bic = stats::BIC, aic = stats::AIC)
      ic1 <- ic_fun(fit1)
      ic0 <- ic_fun(fit0)
      # RSS proxy via deviance (not exact under GLS); rely on logLik for comparisons
      coefs <- tryCatch(stats::coef(fit1), error = function(e) rep(NA_real_, 0))
      amp   <- if (all(c("cos_wt","sin_wt") %in% names(coefs))) sqrt(coefs["cos_wt"]^2 + coefs["sin_wt"]^2) else NA_real_
      phi   <- if (all(c("cos_wt","sin_wt") %in% names(coefs))) atan2(-coefs["sin_wt"], coefs["cos_wt"]) else NA_real_
      list(fit1 = fit1, fit0 = fit0, RSS1 = NA_real_, RSS0 = NA_real_, ll1 = ll1, ll0 = ll0, ic1 = ic1, ic0 = ic0,
           amp = amp, phi = phi, w = w)
    }
  }
  
  # Statistic used for permutation under either method: Likelihood-ratio (LR)
  lrt_stat <- function(ll1, ll0) 2 * (ll1 - ll0)
  
  # Circular-shift permutation within each id by rotating y along time order
  circ_shift_y <- function(y, time, id){
    if (is.null(id)) id <- factor(rep("1", length(y)))
    y_perm <- numeric(length(y))
    for (lev in levels(id)){
      idx <- which(id == lev)
      if (length(idx) <= 1){ y_perm[idx] <- y[idx]; next }
      ord <- idx[order(time[idx], na.last = NA)]
      k <- length(ord)
      s <- if (k > 1) sample.int(k, 1) - 1 else 0  # shift 0..k-1
      y_perm[ord] <- y[ord][ (seq_len(k) - 1 + s) %% k + 1 ]
    }
    y_perm
  }
  
  # Containers
  G <- nrow(Y)
  gene_names <- rownames(Y); if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(G))
  results <- vector("list", G)
  overlay_plots <- if (return_plots) vector("list", G) else NULL
  profile_plots <- if (return_plots) vector("list", G) else NULL
  
  if (save_plots) {
    if (is.null(plot_dir)) stop("Please provide plot_dir when save_plots=TRUE.")
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # Progress helper
  .msg <- function(...) if (verbose) cat(sprintf(...))
  
  # Main loop over genes
  for (g in seq_len(G)){
    y <- as.numeric(Y[g, ])
    ok <- is.finite(y) & is.finite(df$time)
    if (!all(ok)){
      y  <- y[ok]
      df_g <- df[ok, , drop = FALSE]
    } else {
      df_g <- df
    }
    if (length(y) < 6){
      results[[g]] <- data.frame(
        gene = gene_names[g], n = length(y), method = method,
        p_raw = NA_real_, p_perm = NA_real_, p_adj = NA_real_,
        amp = NA_real_, phi = NA_real_, peak_period = NA_real_,
        peak_within_tol = NA, alt_period = NA_real_, delta_ic = NA_real_,
        partial_R2 = NA_real_, lr_stat = NA_real_, is_24h = NA
      )
      next
    }
    
    # 1) Fit at 24h
    fit24 <- fit_period(y, P = 24, method = method, cor_struct = cor_struct, df = df_g, base_rhs = base_rhs)
    # F/LRT p-value with parametric anova when possible; fall back to chi-sq from LR.
    p_raw <- tryCatch({
      if (method == "cosinor") {
        a <- stats::anova(fit24$fit0, fit24$fit1)
        as.numeric(a$`Pr(>F)`[2])
      } else {
        a <- nlme::anova(fit24$fit0, fit24$fit1)
        as.numeric(a$`p-value`[2])
      }
    }, error = function(e){
      # Fallback: LR ~ ChiSq(df=2)
      stats::pchisq(lrt_stat(fit24$ll1, fit24$ll0), df = 2, lower.tail = FALSE)
    })
    
    # partial R^2 (only reliable for cosinor/OLS). For GLS we report NA.
    R2p <- if (!is.na(fit24$RSS0) && !is.na(fit24$RSS1)) partial_R2_from_RSS(fit24$RSS0, fit24$RSS1) else NA_real_
    
    amp <- fit24$amp; phi <- fit24$phi
    
    # 2) Neighborhood profile
    prof_ll <- numeric(length(periods)); names(prof_ll) <- as.character(periods)
    prof_ic <- numeric(length(periods)); names(prof_ic) <- as.character(periods)
    
    for (k in seq_along(periods)){
      Pk <- periods[k]
      fk <- try(fit_period(y, P = Pk, method = method, cor_struct = cor_struct, df = df_g, base_rhs = base_rhs), silent = TRUE)
      if (inherits(fk, "try-error")){
        prof_ll[k] <- NA_real_; prof_ic[k] <- NA_real_
      } else {
        prof_ll[k] <- fk$ll1
        prof_ic[k] <- switch(contrast, bic = stats::BIC(fk$fit1), aic = stats::AIC(fk$fit1))
      }
    }
    
    # Peak period (maximize logLik; tie-breaker = minimize IC)
    peak_idx <- which.max(prof_ll)
    peakP    <- periods[peak_idx]
    # Best alternative excluding 24 (by IC)
    alt_mask <- abs(periods - 24) > .Machine$double.eps
    alt_idx  <- which.min(prof_ic[alt_mask])
    altP     <- periods[alt_mask][alt_idx]
    
    # ΔIC = IC(24) - IC(best_alt)
    ic24 <- prof_ic[which(periods == 24)]
    icalt <- prof_ic[alt_mask][alt_idx]
    delta_ic <- unname(ic24 - icalt)
    
    # Specificity flags
    within_tol <- is.finite(peakP) && abs(peakP - 24) <= peak_tol
    pass_S1 <- isTRUE(within_tol)
    pass_S2 <- is.finite(delta_ic) && (delta_ic <= delta_ic_thresh)
    
    # 3) Optional permutation p-value (LR statistic)
    p_perm <- NA_real_
    lr_obs <- lrt_stat(fit24$ll1, fit24$ll0)
    if (use_permutation){
      ge <- 0L
      for (b in seq_len(B)){
        yb <- circ_shift_y(y, df_g$time, df_g$id)
        fb <- fit_period(yb, P = 24, method = method, cor_struct = cor_struct, df = df_g, base_rhs = base_rhs)
        lr_b <- lrt_stat(fb$ll1, fb$ll0)
        if (is.finite(lr_b) && lr_b >= lr_obs) ge <- ge + 1L
      }
      p_perm <- (ge + 1) / (B + 1)
    }
    
    results[[g]] <- data.frame(
      gene = gene_names[g], n = length(y), method = method,
      p_raw = p_raw, p_perm = p_perm, p_adj = NA_real_,
      amp = amp, phi = phi, peak_period = peakP,
      peak_within_tol = pass_S1, alt_period = altP, delta_ic = delta_ic,
      partial_R2 = R2p, lr_stat = lr_obs, stringsAsFactors = FALSE
    )
    
    # 4) Plots (optional)
    if (save_plots || return_plots){
      # Overlay plot on 24h fit
      # Build fitted line: predict at unique times sorted
      w24 <- 2*pi/24
      XY24 <- cbind(df_g, cos_wt = cos(w24 * df_g$time), sin_wt = sin(w24 * df_g$time))
      if (method == "cosinor") {
        yhat <- stats::predict(fit24$fit1, newdata = cbind(XY24))
      } else {
        yhat <- stats::predict(fit24$fit1, newdata = cbind(XY24))
      }
      ord <- order(df_g$time)
      dd <- data.frame(time = df_g$time, y = y, id = df_g$id, yhat = yhat)
      p_overlay <- ggplot2::ggplot(dd, ggplot2::aes(x = time, y = y)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::geom_line(ggplot2::aes(y = yhat), linewidth = 0.8) +
        ggplot2::labs(title = paste0(gene_names[g], " (fit @ 24h)"), x = "Time (h)", y = "Expression") +
        ggplot2::theme_minimal(base_size = 12)
      
      # Profile plot (logLik and scaled power)
      prof_df <- data.frame(period = periods, logLik = prof_ll, IC = prof_ic)
      # Scaled power 0-1 for display
      rng <- range(prof_ll, na.rm = TRUE)
      prof_df$power01 <- if (is.finite(rng[1]) && diff(rng) > 0) (prof_ll - rng[1]) / diff(rng) else NA_real_
      p_prof <- ggplot2::ggplot(prof_df, ggplot2::aes(x = period, y = logLik)) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = 24, linetype = "dashed") +
        ggplot2::labs(title = paste0(gene_names[g], " neighborhood profile"),
                      x = "Period (h)", y = "logLik (fit with cos/sin @ P)") +
        ggplot2::theme_minimal(base_size = 12)
      
      if (save_plots){
        fn1 <- file.path(plot_dir, paste0(gene_names[g], "__overlay_24h.png"))
        fn2 <- file.path(plot_dir, paste0(gene_names[g], "__profile.png"))
        ggplot2::ggsave(fn1, p_overlay, width = 6.5, height = 4.5, dpi = 300)
        ggplot2::ggsave(fn2, p_prof,    width = 6.5, height = 4.5, dpi = 300)
      }
      if (return_plots){
        overlay_plots[[g]] <- p_overlay
        profile_plots[[g]] <- p_prof
      }
    }
    
    if (verbose && (g %% max(1, floor(G/10)) == 0)) .msg("Processed %d / %d genes\n", g, G)
  }
  
  res_df <- do.call(rbind, results)
  
  # 5) Multiple testing adjustment
  # Choose the primary p-value: use permutation if requested, else raw
  p_primary <- if (use_permutation) res_df$p_perm else res_df$p_raw
  
  if (p_adjust == "none"){
    res_df$p_adj <- NA_real_
  } else if (p_adjust == "qvalue"){
    qobj <- qvalue::qvalue(p = p_primary)
    res_df$p_adj <- qobj$qvalues
  } else {
    res_df$p_adj <- stats::p.adjust(p_primary, method = p_adjust)
  }
  
  # 6) Binary decision: significant & specificity S1/S2
  sig_flag <- if (p_adjust == "none") p_primary <= alpha else res_df$p_adj <= alpha
  spec_flag <- isTRUE(all(!is.na(res_df$peak_within_tol))) # vectorized below
  res_df$is_24h <- (sig_flag) & (res_df$peak_within_tol) & (res_df$delta_ic <= delta_ic_thresh)
  
  # Return
  out <- list(
    results = res_df,
    overlay_plots = if (return_plots) setNames(overlay_plots, gene_names) else NULL,
    profile_plots = if (return_plots) setNames(profile_plots, gene_names) else NULL,
    params = list(method = method, cor_struct = cor_struct, trend_df = trend_df,
                  periods = periods, peak_tol = peak_tol, contrast = contrast,
                  delta_ic_thresh = delta_ic_thresh, use_permutation = use_permutation,
                  B = B, p_adjust = p_adjust, alpha = alpha)
  )
  class(out) <- c("circadian24_screen", class(out))
  out
}