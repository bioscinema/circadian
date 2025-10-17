## ============= Helpers =============

.wrap_angle <- function(x) {
  # wrap to (-pi, pi]
  a <- (x + pi) %% (2*pi)
  a[a == 0] <- 2*pi
  a - pi
}

.unwrap_to_ref <- function(phi, ref) {
  # unwrap phi to be closest to ref (elementwise)
  k <- round((ref - phi) / (2*pi))
  phi + 2*pi*k
}

.make_design <- function(t, omega) {
  cbind(Intercept = 1,
        S = sin(omega * t),
        C = cos(omega * t))
}

.vcov_sandwich <- function(X, resid, type = c("HC3","HC0")) {
  type <- match.arg(type)
  n <- nrow(X); p <- ncol(X)
  XtX_inv <- solve(crossprod(X))
  # leverage for HC3
  h <- rowSums((X %*% XtX_inv) * X)
  u <- as.numeric(resid)
  if (type == "HC0") {
    meat <- t(X) %*% diag(u^2, n, n) %*% X
  } else { # HC3
    w <- u / (1 - h)
    meat <- t(X) %*% diag(w^2, n, n) %*% X
  }
  XtX_inv %*% meat %*% XtX_inv
}

## delta variance for phi = atan2(C, B), given vcov of (B, C)
.var_phi_delta <- function(B, C, vcov_BC) {
  A2 <- B^2 + C^2
  if (A2 <= .Machine$double.eps) return(Inf)
  g <- c(-C, B) / A2
  as.numeric(t(g) %*% vcov_BC %*% g)
}

## ============= 1) Main fitting function =============

fit_circadian <- function(y, t,
                          period = 24, omega = 2*pi/period,
                          vcov_type = c("HC3","HC0","iid")) {
  stopifnot(length(y) == length(t))
  vcov_type <- match.arg(vcov_type)
  X <- .make_design(t, omega)
  fit <- lm.fit(x = X, y = y)
  beta_hat <- coef(fit)              # c, B, C
  names(beta_hat) <- c("c","B","C")
  resid <- fit$residuals
  n <- length(y); p <- ncol(X)
  rss <- sum(resid^2)
  sigma2_iid <- rss / (n - p)
  
  # vcov for (c,B,C)
  if (vcov_type == "iid") {
    vcov_beta <- sigma2_iid * solve(crossprod(X))
  } else {
    vcov_beta <- .vcov_sandwich(X, resid, type = vcov_type)
  }
  vcov_BC <- vcov_beta[2:3, 2:3, drop = FALSE]
  
  B <- beta_hat["B"]; C <- beta_hat["C"]
  A_hat <- sqrt(B^2 + C^2)
  phi_hat <- atan2(C, B)            # (-pi,pi]
  c_hat <- beta_hat["c"]
  
  se_phi <- sqrt(.var_phi_delta(B, C, vcov_BC))
  se_A <- if (A_hat <= .Machine$double.eps) Inf else {
    # delta for A = sqrt(B^2 + C^2)
    gA <- c(B, C) / A_hat
    sqrt(as.numeric(t(gA) %*% vcov_BC %*% gA))
  }
  
  structure(list(
    call = match.call(),
    period = period, omega = omega,
    coef = c(c = c_hat, A = A_hat, phi = phi_hat, B = B, C = C),
    se = c(se_c = sqrt(vcov_beta[1,1]),
           se_A = se_A,
           se_phi = se_phi,
           se_B = sqrt(vcov_BC[1,1]),
           se_C = sqrt(vcov_BC[2,2])),
    vcov_beta = vcov_beta,
    vcov_BC = vcov_BC,
    sigma2_iid = sigma2_iid,
    resid = resid, X = X, y = y, t = t,
    n = n, p = p, rss = rss,
    vcov_type = vcov_type,
    method = "fixed-frequency sinusoid (24h)"
  ), class = "circadian_fit")
}

## ============= 2) Confidence intervals =============

confint_circadian <- function(object,
                              parm = c("phi","A","c"),
                              level = 0.95,
                              method = c("delta","lrt"),
                              grid_size = 1440) {
  stopifnot(inherits(object, "circadian_fit"))
  parm <- match.arg(parm)
  method <- match.arg(method)
  z <- qnorm(0.5 + level/2)
  
  if (parm %in% c("A","c") && method == "lrt") {
    stop("Currently LRT only implemented for phi; use delta for A/c.")
  }
  
  if (parm == "phi") {
    phihat <- object$coef["phi"]
    if (method == "delta") {
      se <- object$se["se_phi"]
      ci <- phihat + c(-1,1) * z * se
      return(.wrap_angle(ci))
    } else {
      # profile-LRT over phi
      omega <- object$omega
      y <- object$y; t <- object$t
      n <- object$n
      sigma2 <- object$sigma2_iid
      
      phi_grid <- seq(phihat - pi, phihat + pi, length.out = grid_size)
      rss_phi <- numeric(length(phi_grid))
      for (i in seq_along(phi_grid)) {
        Sphi <- sin(omega * t + phi_grid[i])
        Xphi <- cbind(1, Sphi)
        bphi <- lm.fit(Xphi, y)
        rss_phi[i] <- sum(bphi$residuals^2)
      }
      rss_min <- min(rss_phi)
      cutoff <- qchisq(level, df = 1) * sigma2
      ok <- (rss_phi - rss_min) <= cutoff
      if (!any(ok)) {
        warning("LRT grid failed; falling back to delta SE.")
        se <- object$se["se_phi"]
        return(.wrap_angle(phihat + c(-1,1)*z*se))
      }
      idx <- which(ok)
      j0 <- idx[which.min(abs(phi_grid[idx] - phihat))]
      L <- R <- j0
      while (L > 1 && ok[L-1]) L <- L - 1
      while (R < length(ok) && ok[R+1]) R <- R + 1
      ci <- c(phi_grid[L], phi_grid[R])
      return(.wrap_angle(ci))
    }
  }
  
  if (parm == "A") {
    Ahat <- object$coef["A"]; seA <- object$se["se_A"]
    return(Ahat + c(-1,1) * z * seA)
  } else if (parm == "c") {
    chat <- object$coef["c"]; sec <- object$se["se_c"]
    return(chat + c(-1,1) * z * sec)
  }
}

## ============= 3) Contrast tests =============

phase_contrast_test <- function(fits, L, rhs = NULL,
                                alternative = c("two.sided","less","greater")) {
  if (inherits(fits, "circadian_fit")) fits <- list(fits)
  m <- length(fits)
  alternative <- match.arg(alternative)
  if (is.null(rhs)) rhs <- rep(0, nrow(L))
  stopifnot(ncol(L) == m, length(rhs) == nrow(L))
  
  phi <- sapply(fits, function(f) unname(f$coef["phi"]))
  u <- c(mean(cos(phi)), mean(sin(phi)))
  phi_center <- atan2(u[2], u[1])
  phi_u <- .unwrap_to_ref(phi, ref = phi_center)
  
  V <- diag(sapply(fits, function(f) unname(f$se["se_phi"]^2)), m, m)
  r <- nrow(L)
  diff <- as.numeric(L %*% phi_u - rhs)
  LVLt_inv <- solve(L %*% V %*% t(L))
  if (r == 1) {
    se <- sqrt(as.numeric(L %*% V %*% t(L)))
    z <- diff / se
    p <- switch(alternative,
                "two.sided" = 2*pnorm(-abs(z)),
                "less"      = pnorm(z),
                "greater"   = 1 - pnorm(z))
    out <- list(statistic = z, p.value = p, df = 1,
                estimate = as.numeric(L %*% phi_u),
                se = se, alternative = alternative,
                method = "Wald z-test for phase contrast (r=1)",
                contrast = L, rhs = rhs, phis_unwrapped = phi_u)
  } else {
    W <- as.numeric(t(diff) %*% LVLt_inv %*% diff)
    p <- pchisq(W, df = r, lower.tail = FALSE)
    out <- list(statistic = W, p.value = p, df = r,
                estimate = as.numeric(L %*% phi_u),
                vcov_contrast = L %*% V %*% t(L),
                method = "Wald chi-square test for phase contrasts (r>1)",
                contrast = L, rhs = rhs, phis_unwrapped = phi_u)
  }
  class(out) <- "phase_contrast_test"
  out
}

test_phase <- function(fit, phi0 = 0,
                       alternative = c("two.sided","less","greater")) {
  alternative <- match.arg(alternative)
  phi_hat <- unname(fit$coef["phi"])
  phi0_u <- .unwrap_to_ref(phi0, ref = phi_hat)
  se <- unname(fit$se["se_phi"])
  z <- (phi_hat - phi0_u) / se
  p <- switch(alternative,
              "two.sided" = 2*pnorm(-abs(z)),
              "less"      = pnorm(z),
              "greater"   = 1 - pnorm(z))
  list(statistic = z, p.value = p, df = 1,
       estimate = phi_hat, se = se, phi0 = phi0,
       alternative = alternative,
       method = "Wald z-test for single-gene phase")
}

test_phase_diff <- function(fit1, fit2, delta0 = 0,
                            alternative = c("two.sided","less","greater"),
                            method = c("difference","crossprod")) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  
  if (method == "difference") {
    phi1 <- unname(fit1$coef["phi"]); se1 <- unname(fit1$se["se_phi"])
    phi2 <- unname(fit2$coef["phi"]); se2 <- unname(fit2$se["se_phi"])
    target <- phi2 - phi1
    target <- .unwrap_to_ref(target, ref = delta0)
    se <- sqrt(se1^2 + se2^2)
    z <- (target - delta0) / se
    p <- switch(alternative,
                "two.sided" = 2*pnorm(-abs(z)),
                "less"      = pnorm(z),
                "greater"   = 1 - pnorm(z))
    return(list(statistic = z, p.value = p, df = 1,
                estimate = target, se = se, delta0 = delta0,
                alternative = alternative,
                method = "Wald z-test for phase difference"))
  } else {
    B1 <- unname(fit1$coef["B"]); C1 <- unname(fit1$coef["C"])
    B2 <- unname(fit2$coef["B"]); C2 <- unname(fit2$coef["C"])
    A1 <- unname(fit1$coef["A"]); A2 <- unname(fit2$coef["A"])
    U <- (B1*C2 - C1*B2) - A1*A2*sin(delta0)
    V1 <- fit1$vcov_BC; V2 <- fit2$vcov_BC
    g <- c( C2, -B2, -C1,  B1)
    VU <- as.numeric(t(g) %*% bdiag(V1, V2) %*% g)
    z <- U / sqrt(VU)
    p <- switch(alternative,
                "two.sided" = 2*pnorm(-abs(z)),
                "less"      = pnorm(z),
                "greater"   = 1 - pnorm(z))
    return(list(statistic = z, p.value = p, df = 1,
                estimate = U, se = sqrt(VU), delta0 = delta0,
                alternative = alternative,
                method = "Wald z-test via cross-product invariant"))
  }
}

## ============= 4) Prediction function =============

predict_circadian <- function(object, new_t,
                              interval = c("none","confidence","prediction"),
                              level = 0.95) {
  stopifnot(inherits(object, "circadian_fit"))
  interval <- match.arg(interval)
  Xnew <- .make_design(new_t, object$omega)
  beta_hat <- c(object$coef["c"], object$coef["B"], object$coef["C"])
  yhat <- as.numeric(Xnew %*% beta_hat)
  out <- list(fit = yhat)
  
  if (interval != "none") {
    V <- object$vcov_beta
    se_mean <- sqrt(rowSums((Xnew %*% V) * Xnew))
    z <- qnorm(0.5 + level/2)
    ci <- cbind(yhat - z*se_mean, yhat + z*se_mean)
    colnames(ci) <- c("lwr","upr")
    if (interval == "confidence") {
      out$conf_int <- ci
    } else {
      se_pred <- sqrt(se_mean^2 + object$sigma2_iid)
      pi <- cbind(yhat - z*se_pred, yhat + z*se_pred)
      colnames(pi) <- c("lwr","upr")
      out$pred_int <- pi
    }
  }
  out
}
