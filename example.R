## setting 1: 50% 24-hour and 50% non-periodic
## conclusion: works for uniform sampling, fails for non-uniform sampling

source("circadian_screen.R")
source("DGP.R")
source("circadian_analysis.R")
source("freq_domain_analysis.R")

#######################################
### case 1: uniform sampling time period
########################################
### data generating (to capture period > 24, need longer observation window)

sim1 <- simulate_mix_periods(
  m = 1000,
  # ---- time settings ----
  times = NULL,                        
  time_mode = "uniform",
  uniform_start = 0, uniform_by = 1, uniform_length = 48,
  nonuniform_start = 0, nonuniform_end = 24, nonuniform_n = 24,
  periods = c(36, 24, 12, NA),     
  frac = c(0.2, 0.2, 0.2, 0.4),                  
  A_range = c(0.8, 2.0),
  c_mean_sd = c(10, 2),
  sigma = 1.0,
  seed = 2025)

### scope screening

scope.fit <- rep(0,1000);names(scope.fit) <- row.names(sim1$expr)
for (i in 1:1000) {
  out1 <- spectrum_with_detrend(
    t = sim1$time,
    y = sim1$expr[i,],
    detrend="diff", 
    d_order=1,
    method="lomb", 
    window="hann", 
    correct=TRUE,
    plot = F
  )
  scope.fit[i] <- out1$spectrum$freq[which.max(out1$spectrum$P_corr)]
  if(i%%100 == 0){print(i)}
}

coverage0 <- sum(scope.fit < 0.04 & sim1$period_true == 36,na.rm = T)/sum(sim1$period_true == 36,na.rm = T)
TDR0 <- sum(scope.fit < 0.04 & sim1$period_true == 36,na.rm = T)/sum(scope.fit < 0.04)
FDR0 <- sum(scope.fit < 0.04 & sim1$period_true == 24,na.rm = T)/sum(scope.fit < 0.04)

to_be_removed <- names(scope.fit[scope.fit<0.04])
rownames.vec <- rownames(sim1$expr)
sim1$expr <- sim1$expr[!rownames.vec %in% to_be_removed, ]
sim1$phi <- sim1$phi[!rownames.vec %in% to_be_removed]
sim1$period_true <- sim1$period_true[!rownames.vec %in% to_be_removed]

### period screening
fit <- circadian_screen(
  Y = sim1$expr, time = sim1$time, covariates = NULL,
  method = "cosinor",          
  cor_struct = "none",
  trend_df = 0,                 
  periods = seq(23, 25, by = 0.025),    
  peak_tol = 1.2,               
  contrast = "bic",
  delta_ic_thresh = Inf,         
  use_permutation = FALSE,
  p_adjust = "none",            
  alpha = 0.05,
  save_plots = FALSE, return_plots = TRUE, verbose = TRUE
)
head(fit$results)

true_index <- which(sim1$period_true == 24)
selected_index <- which(fit$results$is_24h == T)
TDR <- mean(selected_index %in% true_index)
Coverage <- sum(selected_index %in% true_index)/length(true_index)
plot(fit$results$phi~sim1$phi)

summary(sim1$period_true[fit$results$is_24h == T])
hist(sim1$period_true[fit$results$is_24h == T])

### phase estimating
Y_filtered <- t(sim1$expr[selected_index,])
time_filtered <- sim1$time

fit.filtered <- apply(Y_filtered, 2, function(y) fit_circadian(y, t = time_filtered))
phi_estimates <- sapply(fit.filtered, function(f) f$coef["phi.C"])
phi_true <- sim1$phi[selected_index]

plot(phi_estimates~phi_true)

