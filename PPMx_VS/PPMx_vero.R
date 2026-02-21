# ==============================================================================
# PPMx con Selezione Variabile Cluster-Specifica (logit-Normal)
# Implementazione fedele a Quintana et al. (2015)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Marginal Predictive Log-Density (Normal-Inverse-Gamma)
# ------------------------------------------------------------------------------
# Calcola la verosimiglianza marginale integrando via i parametri locali.
# Utilizzata per determinare la similarità tra osservazioni.
calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  njj <- length(z_j)
  
  if (njj == 0) {
    # Caso cluster vuoto: predittiva basata solo sulla misura di base P0
    df_loc     <- 2 * alpha0
    mu_loc     <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  } else {
    # Aggiornamento Bayesiano dei parametri per il cluster esistente
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    ss_j <- sum(z_j^2) - (sum_z^2) / njj
    beta_loc <- beta0 + 0.5 * ss_j + 0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    
    df_loc     <- 2 * alpha_loc
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) / (alpha_loc * kappa_loc)
  }
  
  sig_loc <- sqrt(sig_sq_loc)
  # Ritorna la densità log-t (Student)
  return(dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc))
}

# ------------------------------------------------------------------------------
# 2. Algoritmo MCMC principale con Variable Selection
# ------------------------------------------------------------------------------

PPMx_VS_Quintana_LogitNormal <- function(
    y, x,
    n_iter, burn_in, thin,
    hyperpar_P0_y, hyperpar_P0_x,
    hyperpar_Dir, hyperpar_VS,
    mle_xi_list, init
) {
  
  n      <- length(y)
  n_cov  <- ncol(x)
  saved_iters <- seq(burn_in + 1, n_iter, by = thin)
  G <- length(saved_iters)
  
  # --- Inizializzazione ---
  clu_curr <- as.integer(init$clu)
  k_curr   <- length(unique(clu_curr))
  
  # z_{jl}: Variabile latente logit-Normal per la selezione (ex eta)
  z_jl_curr   <- matrix(0, k_curr, n_cov)      
  mu_z_curr   <- rep(hyperpar_VS$mu_eta0, n_cov)
  
  sigma_z  <- hyperpar_VS$sigma_eta
  sigma_z0 <- hyperpar_VS$sigma_eta0
  
  # gamma_{jl}: Indicatore binario di inclusione covariata l nel cluster j
  gamma_jl_curr <- matrix(
    rbinom(k_curr * n_cov, 1, plogis(z_jl_curr)),
    k_curr, n_cov
  )
  
  # Output storage
  k_out      <- rep(NA, G)
  clu_out    <- matrix(NA, G, n)
  gamma_out  <- vector("list", G)
  mu_z_out   <- matrix(NA, G, n_cov)
  
  # --- Inizio Loop MCMC ---
  for (g in 1:n_iter) {
    
    # STEP 1: Allocazione Cluster (Neal Algorithm 8)
    for (i in 1:n) {
      y_i <- y[i]
      x_i <- x[i, ]
      old_c <- clu_curr[i]
      clu_curr[i] <- -1
      
      # Rimozione cluster vuoti e pulizia matrici VS
      if (!any(clu_curr == old_c)) {
        gamma_jl_curr <- gamma_jl_curr[-old_c, , drop = FALSE]
        z_jl_curr     <- z_jl_curr[-old_c, , drop = FALSE]
        clu_curr[clu_curr > old_c] <- clu_curr[clu_curr > old_c] - 1
      }
      
      k_curr <- length(unique(clu_curr[clu_curr > 0]))
      log_w  <- numeric(k_curr + 1)
      
      # Calcolo pesi per cluster esistenti
      for (j in 1:k_curr) {
        idx_j <- which(clu_curr == j)
        log_w[j] <- log(length(idx_j)) + calcola_log_pred_dens_t(y_i, y[idx_j], hyperpar_P0_y)
        
        # Aggiunta contributo covariate informative (gamma == 1)
        for (l in 1:n_cov) {
          if (gamma_jl_curr[j, l] == 1) {
            log_w[j] <- log_w[j] + 
              calcola_log_pred_dens_t(x_i[l], x[idx_j, l], hyperpar_P0_x[[l]]) -
              calcola_log_densita_mle(x_i[l], l, mle_xi_list)
          }
        }
      }
      
      # Peso per nuovo cluster (misura di base)
      log_w[k_curr + 1] <- log(hyperpar_Dir$M) + calcola_log_pred_dens_t(y_i, numeric(0), hyperpar_P0_y)
      
      # Sampling dell'allocazione
      w <- exp(log_w - max(log_w))
      new_c <- sample.int(k_curr + 1, 1, prob = w)
      clu_curr[i] <- new_c
      
      if (new_c > k_curr) {
        # Espansione parametri VS per il nuovo cluster
        gamma_jl_curr <- rbind(gamma_jl_curr, rbinom(n_cov, 1, plogis(mu_z_curr)))
        z_jl_curr     <- rbind(z_jl_curr, rnorm(n_cov, mu_z_curr, sigma_z))
        k_curr <- k_curr + 1
      }
    }
    
    # STEP 2: Aggiornamento gamma_{jl} (Selezione Variabili)
    # Implementa la "Calibrated Similarity": Eq. 6 Quintana et al.
    for (j in 1:k_curr) {
      idx_j <- which(clu_curr == j)
      for (l in 1:n_cov) {
        x_jl <- x[idx_j, l]
        
        # log_g_tilde: Rapporto tra marginal likelihood del cluster e marginal likelihood baseline
        log_g_tilde <- sum(sapply(1:length(x_jl), function(m) 
          calcola_log_pred_dens_t(x_jl[m], x_jl[1:m][-m], hyperpar_P0_x[[l]]))) -
          sum(sapply(x_jl, function(xi) calcola_log_densita_mle(xi, l, mle_xi_list)))
        
        # Posterior odds per gamma (Bernoulli)
        log_odds <- z_jl_curr[j, l] + log_g_tilde
        gamma_jl_curr[j, l] <- rbinom(1, 1, plogis(log_odds))
      }
    }
    
    # STEP 3: Aggiornamento z_{jl} (Metropolis-Hastings)
    # Aggiorna la variabile latente logit-Normal che guida gamma
    for (j in 1:k_curr) {
      for (l in 1:n_cov) {
        z_old <- z_jl_curr[j, l]
        z_new <- rnorm(1, z_old, hyperpar_VS$mh_sd_eta)
        
        log_acc <- (dbinom(gamma_jl_curr[j, l], 1, plogis(z_new), log = TRUE) + 
                      dnorm(z_new, mu_z_curr[l], sigma_z, log = TRUE)) -
          (dbinom(gamma_jl_curr[j, l], 1, plogis(z_old), log = TRUE) + 
             dnorm(z_old, mu_z_curr[l], sigma_z, log = TRUE))
        
        if (log(runif(1)) < log_acc) z_jl_curr[j, l] <- z_new
      }
    }
    
    # STEP 4: Aggiornamento mu_{z,l} (Gibbs Sampling)
    # Aggiorna la media globale della probabilità di inclusione per covariata l
    for (l in 1:n_cov) {
      var_post <- 1 / (k_curr / sigma_z^2 + 1 / sigma_z0^2)
      mean_post <- var_post * (sum(z_jl_curr[, l]) / sigma_z^2 + hyperpar_VS$mu_eta0 / sigma_z0^2)
      mu_z_curr[l] <- rnorm(1, mean_post, sqrt(var_post))
    }
    
    # Salvataggio iterazioni
    if (g %in% saved_iters) {
      idx <- which(saved_iters == g)
      k_out[idx]      <- k_curr
      clu_out[idx, ]  <- clu_curr
      gamma_out[[idx]] <- gamma_jl_curr
      mu_z_out[idx, ] <- mu_z_curr
    }
    if(g %% 50 == 0) cat("Iterazione:", g, " | Cluster attivi:", k_curr, "\n")
  }
  
  return(list(k_out = k_out, clu_out = clu_out, gamma_out = gamma_out, mu_z_out = mu_z_out))
}