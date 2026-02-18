# ============================================================
# PPMx con selezione cluster-specifica (logit-Normal)
# Implementazione fedele a Quintana et al. (2015)
# ============================================================

# ------------------------------------------------------------
# Marginal predictive log-density (Normal-Inverse-Gamma)
# ------------------------------------------------------------
calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  njj <- length(z_j)
  
  if (njj == 0) {
    df_loc     <- 2 * alpha0
    mu_loc     <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  } else {
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    ss_j <- sum(z_j^2) - (sum_z^2) / njj
    beta_loc <- beta0 +
      0.5 * ss_j +
      0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    df_loc     <- 2 * alpha_loc
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) /
      (alpha_loc * kappa_loc)
  }
  
  sig_loc <- sqrt(sig_sq_loc)
  return(
    dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) -
      log(sig_loc)
  )
}

if(TRUE){
# ------------------------------------------------------------
# PPMx con selezione variabili logit-Normal (paper-faithful)
# ------------------------------------------------------------
PPMx_VS_Quintana_LogitNormal <- function(
    y, x,
    n_iter, burn_in, thin,
    hyperpar_P0_y,
    hyperpar_P0_x,
    hyperpar_Dir,
    hyperpar_VS,
    mle_xi_list,
    init
) {
  
  n      <- length(y)
  n_cov  <- ncol(x)
  saved_iters <- seq(burn_in + 1, n_iter, by = thin)
  G <- length(saved_iters)
  
  # ----------------------------------------------------------
  # Inizializzazione
  # ----------------------------------------------------------
  clu_curr <- as.integer(init$clu)
  k_curr   <- length(unique(clu_curr))
  
  # logit-Normal gerarchico (eq. 5)
  eta_jl_curr <- matrix(0, k_curr, n_cov)      # η_{jℓ}
  mu_eta_curr <- rep(hyperpar_VS$mu_eta0, n_cov)
  
  sigma_eta  <- hyperpar_VS$sigma_eta
  sigma_eta0 <- hyperpar_VS$sigma_eta0
  
  gamma_jl_curr <- matrix(
    rbinom(k_curr * n_cov, 1, plogis(eta_jl_curr)),
    k_curr, n_cov
  )
  
  # Output
  k_out    <- rep(NA, G)
  clu_out  <- matrix(NA, G, n)
  gamma_out <- vector("list", G)
  mu_eta_out <- matrix(NA, G, n_cov)
  
  # ----------------------------------------------------------
  # MCMC
  # ----------------------------------------------------------
  for (g in 1:n_iter) {
    
    # ========================================================
    # 1. Allocazione cluster (Neal Algorithm 8)
    # ========================================================
    for (i in 1:n) {
      y_i <- y[i]
      x_i <- x[i, ]
      
      old_c <- clu_curr[i]
      clu_curr[i] <- -1
      
      # rimuove cluster vuoti
      if (!any(clu_curr == old_c)) {
        gamma_jl_curr <- gamma_jl_curr[-old_c, , drop = FALSE]
        eta_jl_curr   <- eta_jl_curr[-old_c, , drop = FALSE]
        clu_curr[clu_curr > old_c] <- clu_curr[clu_curr > old_c] - 1
      }
      
      k_curr <- length(unique(clu_curr[clu_curr > 0]))
      log_w  <- numeric(k_curr + 1)
      
      for (j in 1:k_curr) {
        idx_j <- which(clu_curr == j)
        
        log_w[j] <- log(length(idx_j)) +
          calcola_log_pred_dens_t(y_i, y[idx_j], hyperpar_P0_y)
        
        for (l in 1:n_cov) {
          if (gamma_jl_curr[j, l] == 1) {
            log_w[j] <- log_w[j] +
              calcola_log_pred_dens_t(
                x_i[l],
                x[idx_j, l],
                hyperpar_P0_x[[l]]
              ) -
              calcola_log_densita_mle(x_i[l], l, mle_xi_list)
          }
        }
      }
      
      # nuovo cluster
      log_w[k_curr + 1] <-
        log(hyperpar_Dir$M) +
        calcola_log_pred_dens_t(y_i, numeric(0), hyperpar_P0_y)
      
      w <- exp(log_w - max(log_w))
      new_c <- sample.int(k_curr + 1, 1, prob = w)
      clu_curr[i] <- new_c
      
      if (new_c > k_curr) {
        gamma_jl_curr <- rbind(
          gamma_jl_curr,
          rbinom(n_cov, 1, plogis(mu_eta_curr))
        )
        eta_jl_curr <- rbind(
          eta_jl_curr,
          rnorm(n_cov, mu_eta_curr, sigma_eta)
        )
        k_curr <- k_curr + 1
      }
    }
    
    # ========================================================
    # 2. Aggiornamento gamma_{jℓ} (eq. 4–6)
    # ========================================================
    for (j in 1:k_curr) {
      idx_j <- which(clu_curr == j)
      for (l in 1:n_cov) {
        x_jl <- x[idx_j, l]
        
        log_g_tilde <-
          sum(sapply(
            1:length(x_jl),
            function(m)
              calcola_log_pred_dens_t(
                x_jl[m],
                x_jl[1:m][-m],
                hyperpar_P0_x[[l]]
              )
          )) -
          sum(sapply(
            x_jl,
            function(xi)
              calcola_log_densita_mle(xi, l, mle_xi_list)
          ))
        
        logit_p <- eta_jl_curr[j, l]
        log_odds <- logit_p + log_g_tilde
        
        gamma_jl_curr[j, l] <- rbinom(1, 1, plogis(log_odds))
      }
    }
    
    # ========================================================
    # 3. Aggiornamento η_{jℓ} (Metropolis–Hastings)
    # ========================================================
    for (j in 1:k_curr) {
      for (l in 1:n_cov) {
        eta_old <- eta_jl_curr[j, l]
        eta_new <- rnorm(1, eta_old, hyperpar_VS$mh_sd_eta)
        
        log_prior_old <- dnorm(
          eta_old, mu_eta_curr[l], sigma_eta, log = TRUE
        )
        log_prior_new <- dnorm(
          eta_new, mu_eta_curr[l], sigma_eta, log = TRUE
        )
        
        log_like_old <- dbinom(
          gamma_jl_curr[j, l], 1, plogis(eta_old), log = TRUE
        )
        log_like_new <- dbinom(
          gamma_jl_curr[j, l], 1, plogis(eta_new), log = TRUE
        )
        
        log_acc <- (log_like_new + log_prior_new) -
          (log_like_old + log_prior_old)
        
        if (log(runif(1)) < log_acc) {
          eta_jl_curr[j, l] <- eta_new
        }
      }
    }
    
    # ========================================================
    # 4. Aggiornamento μ_{η,ℓ} (Gibbs – eq. 5)
    # ========================================================
    for (l in 1:n_cov) {
      var_post <- 1 / (k_curr / sigma_eta^2 + 1 / sigma_eta0^2)
      mean_post <- var_post * (
        sum(eta_jl_curr[, l]) / sigma_eta^2 +
          hyperpar_VS$mu_eta0 / sigma_eta0^2
      )
      mu_eta_curr[l] <- rnorm(1, mean_post, sqrt(var_post))
    }
    
    # ========================================================
    # Salvataggio
    # ========================================================
    if (g %in% saved_iters) {
      idx <- which(saved_iters == g)
      k_out[idx]      <- k_curr
      clu_out[idx, ]  <- clu_curr
      gamma_out[[idx]] <- gamma_jl_curr
      mu_eta_out[idx, ] <- mu_eta_curr
    }
    if(g %% 10 == 0) cat("Iterazione:", g, " | Cluster attivi:", k_curr, "\n")
  }
  
  return(list(
    k_out = k_out,
    clu_out = clu_out,
    gamma_out = gamma_out,
    mu_eta_out = mu_eta_out
  ))
}
}
# ============================================================
# PPMx - VERSIONE SEMPLIFICATA (Senza Variable Selection)
# Per test di funzionamento e performance
# ============================================================

if(FALSE){
PPMx_VS_Quintana_LogitNormal <- function(
    y, x,
    n_iter, burn_in, thin,
    hyperpar_P0_y,
    hyperpar_P0_x,
    hyperpar_Dir,
    hyperpar_VS,
    mle_xi_list,
    init
) {
  
  n      <- length(y)
  n_cov  <- ncol(x)
  saved_iters <- seq(burn_in + 1, n_iter, by = thin)
  G <- length(saved_iters)
  
  # ----------------------------------------------------------
  # Inizializzazione
  # ----------------------------------------------------------
  clu_curr <- as.integer(init$clu)
  k_curr   <- length(unique(clu_curr))
  
  # Inizializziamo gamma a 1 (tutte le variabili incluse)
  # Questa matrice non verrà più aggiornata
  gamma_jl_curr <- matrix(1, k_curr, n_cov) 
  
  # --- COMMENTATO: Parametri latenti VS ---
  # eta_jl_curr <- matrix(0, k_curr, n_cov)
  # mu_eta_curr <- rep(hyperpar_VS$mu_eta0, n_cov)
  # sigma_eta  <- hyperpar_VS$sigma_eta
  # sigma_eta0 <- hyperpar_VS$sigma_eta0
  
  # Output
  k_out    <- rep(NA, G)
  clu_out  <- matrix(NA, G, n)
  gamma_out <- vector("list", G)
  mu_eta_out <- matrix(NA, G, n_cov)
  
  # ----------------------------------------------------------
  # MCMC
  # ----------------------------------------------------------
  for (g in 1:n_iter) {
    
    # ========================================================
    # 1. Allocazione cluster (Neal Algorithm 8)
    # ========================================================
    for (i in 1:n) {
      y_i <- y[i]
      x_i <- x[i, ]
      
      old_c <- clu_curr[i]
      clu_curr[i] <- -1
      
      if (!any(clu_curr == old_c)) {
        gamma_jl_curr <- gamma_jl_curr[-old_c, , drop = FALSE]
        # eta_jl_curr   <- eta_jl_curr[-old_c, , drop = FALSE] # COMMENTATO
        clu_curr[clu_curr > old_c] <- clu_curr[clu_curr > old_c] - 1
      }
      
      k_curr <- length(unique(clu_curr[clu_curr > 0]))
      log_w  <- numeric(k_curr + 1)
      
      for (j in 1:k_curr) {
        idx_j <- which(clu_curr == j)
        
        log_w[j] <- log(length(idx_j)) +
          calcola_log_pred_dens_t(y_i, y[idx_j], hyperpar_P0_y)
        
        # In questa versione semplificata gamma_jl_curr è sempre 1
        for (l in 1:n_cov) {
          if (gamma_jl_curr[j, l] == 1) {
            log_w[j] <- log_w[j] +
              calcola_log_pred_dens_t(
                x_i[l],
                x[idx_j, l],
                hyperpar_P0_x[[l]]
              ) -
              calcola_log_densita_mle(x_i[l], l, mle_xi_list)
          }
        }
      }
      
      # nuovo cluster
      log_w[k_curr + 1] <-
        log(hyperpar_Dir$M) +
        calcola_log_pred_dens_t(y_i, numeric(0), hyperpar_P0_y)
      
      w <- exp(log_w - max(log_w))
      new_c <- sample.int(k_curr + 1, 1, prob = w)
      clu_curr[i] <- new_c
      
      if (new_c > k_curr) {
        # Aggiungiamo una riga di 1 per il nuovo cluster
        gamma_jl_curr <- rbind(gamma_jl_curr, rep(1, n_cov))
        # eta_jl_curr <- rbind(eta_jl_curr, rnorm(n_cov, 0, 1)) # COMMENTATO
        k_curr <- k_curr + 1
      }
    }
    
    # ========================================================
    # 2. Aggiornamento gamma_{jℓ} -> DISATTIVATO
    # ========================================================
    # (Questa era la parte più lenta del codice originale)
    # for (j in 1:k_curr) { ... } 
    
    # ========================================================
    # 3. Aggiornamento η_{jℓ} -> DISATTIVATO
    # ========================================================
    # for (j in 1:k_curr) { ... }
    
    # ========================================================
    # 4. Aggiornamento μ_{η,ℓ} -> DISATTIVATO
    # ========================================================
    # for (l in 1:n_cov) { ... }
    
    # ========================================================
    # Salvataggio
    # ========================================================
    if (g %in% saved_iters) {
      idx <- which(saved_iters == g)
      k_out[idx]      <- k_curr
      clu_out[idx, ]  <- clu_curr
      gamma_out[[idx]] <- gamma_jl_curr
      # mu_eta_out[idx, ] <- mu_eta_curr # COMMENTATO
    }
    
    # Feedback visivo per il debug
    if(g %% 10 == 0) cat("Iterazione:", g, "Cluster attivi:", k_curr, "\n")
  }
  
  return(list(
    k_out = k_out,
    clu_out = clu_out,
    gamma_out = gamma_out
    # mu_eta_out = mu_eta_out # COMMENTATO
  ))
}
}


# ============================================================
#ALGORITMO CON VARIABLE SELECTION EFFICIENTE
# Per test di funzionamento e performance
# ============================================================
if(FALSE){
  # ------------------------------------------------------------
  # PPMx con selezione variabili logit-Normal (paper-faithful)
  # ------------------------------------------------------------
  PPMx_VS_Quintana_LogitNormal <- function(
    y, x,
    n_iter, burn_in, thin,
    hyperpar_P0_y,
    hyperpar_P0_x,
    hyperpar_Dir,
    hyperpar_VS,
    mle_xi_list,
    init
  ) {
    
    n      <- length(y)
    n_cov  <- ncol(x)
    saved_iters <- seq(burn_in + 1, n_iter, by = thin)
    G <- length(saved_iters)
    
    # ----------------------------------------------------------
    # Inizializzazione
    # ----------------------------------------------------------
    clu_curr <- as.integer(init$clu)
    k_curr   <- length(unique(clu_curr))
    
    # logit-Normal gerarchico (eq. 5)
    eta_jl_curr <- matrix(0, k_curr, n_cov)      # η_{jℓ}
    mu_eta_curr <- rep(hyperpar_VS$mu_eta0, n_cov)
    
    sigma_eta  <- hyperpar_VS$sigma_eta
    sigma_eta0 <- hyperpar_VS$sigma_eta0
    
    gamma_jl_curr <- matrix(
      rbinom(k_curr * n_cov, 1, plogis(eta_jl_curr)),
      k_curr, n_cov
    )
    
    # Output
    k_out    <- rep(NA, G)
    clu_out  <- matrix(NA, G, n)
    gamma_out <- vector("list", G)
    mu_eta_out <- matrix(NA, G, n_cov)
    
    # ----------------------------------------------------------
    # MCMC
    # ----------------------------------------------------------
    for (g in 1:n_iter) {
      
      # ========================================================
      # 1. Allocazione cluster (Neal Algorithm 8)
      # ========================================================
      for (i in 1:n) {
        y_i <- y[i]
        x_i <- x[i, ]
        
        old_c <- clu_curr[i]
        clu_curr[i] <- -1
        
        # rimuove cluster vuoti
        if (!any(clu_curr == old_c)) {
          gamma_jl_curr <- gamma_jl_curr[-old_c, , drop = FALSE]
          eta_jl_curr   <- eta_jl_curr[-old_c, , drop = FALSE]
          clu_curr[clu_curr > old_c] <- clu_curr[clu_curr > old_c] - 1
        }
        
        k_curr <- length(unique(clu_curr[clu_curr > 0]))
        log_w  <- numeric(k_curr + 1)
        
        for (j in 1:k_curr) {
          idx_j <- which(clu_curr == j)
          
          log_w[j] <- log(length(idx_j)) +
            calcola_log_pred_dens_t(y_i, y[idx_j], hyperpar_P0_y)
          
          for (l in 1:n_cov) {
            if (gamma_jl_curr[j, l] == 1) {
              log_w[j] <- log_w[j] +
                calcola_log_pred_dens_t(
                  x_i[l],
                  x[idx_j, l],
                  hyperpar_P0_x[[l]]
                ) -
                calcola_log_densita_mle(x_i[l], l, mle_xi_list)
            }
          }
        }
        
        # nuovo cluster
        log_w[k_curr + 1] <-
          log(hyperpar_Dir$M) +
          calcola_log_pred_dens_t(y_i, numeric(0), hyperpar_P0_y)
        
        w <- exp(log_w - max(log_w))
        new_c <- sample.int(k_curr + 1, 1, prob = w)
        clu_curr[i] <- new_c
        
        if (new_c > k_curr) {
          gamma_jl_curr <- rbind(
            gamma_jl_curr,
            rbinom(n_cov, 1, plogis(mu_eta_curr))
          )
          eta_jl_curr <- rbind(
            eta_jl_curr,
            rnorm(n_cov, mu_eta_curr, sigma_eta)
          )
          k_curr <- k_curr + 1
        }
      }
      
      # ========================================================
      # 2. Aggiornamento gamma_{jℓ} (OTTIMIZZATO - Fedele al PDF)
      # ========================================================
      # Implementa l'Eq. 6 del PDF e l'Eq. 4 di Quintana et al.
      for (j in 1:k_curr) {
        idx_j <- which(clu_curr == j)
        njj <- length(idx_j)
        
        for (l in 1:n_cov) {
          x_jl <- x[idx_j, l]
          hp <- hyperpar_P0_x[[l]]
          
          # --- 1. Calcolo log g_tilde (Calibrated Similarity) ---
          # Formula della Marginal Likelihood NIG (Normal-Inverse-Gamma)
          sum_x <- sum(x_jl)
          ss_x  <- sum(x_jl^2)
          
          # Parametri posterior (come da sezione 2.1 del PDF)
          kn <- hp$kappa0 + njj
          an <- (hp$nu0 + njj) / 2
          bn <- (hp$nu0 * hp$sig_square0 + ss_x + hp$kappa0 * hp$mu0^2 - (hp$kappa0 * hp$mu0 + sum_x)^2 / kn) / 2
          
          # Log-Marginal Likelihood (g_l) [cite: 589]
          log_marg_P0 <- (hp$nu0/2) * log(hp$nu0 * hp$sig_square0 / 2) - 
            (njj/2) * log(pi) + 0.5 * log(hp$kappa0 / kn) + 
            lgamma(an) - lgamma(hp$nu0/2) - an * log(bn)
          
          # Log-Marginal Baseline (m_l) - Prodotto delle marginali indipendenti
          # Come richiesto dalla "Calibrated Similarity" (Eq. 6 PDF) [cite: 356, 426]
          log_marg_baseline <- sum(sapply(x_jl, function(xi) 
            calcola_log_densita_mle(xi, l, mle_xi_list)))
          
          log_g_tilde_jl <- log_marg_P0 - log_marg_baseline
          
          # --- 2. Sampling gamma_{jl} ---
          # Eq. 8 del PDF: gamma | eta ~ Bernoulli(logit^-1(eta + log_g_tilde))
          logit_p <- eta_jl_curr[j, l] + log_g_tilde_jl
          
          # Protezione numerica per evitare NaN
          if (is.na(logit_p)) logit_p <- -100 
          
          gamma_jl_curr[j, l] <- rbinom(1, 1, plogis(logit_p))
        }
      }
      
      # ========================================================
      # 3. Aggiornamento η_{jℓ} (Metropolis–Hastings)
      # ========================================================
      for (j in 1:k_curr) {
        for (l in 1:n_cov) {
          eta_old <- eta_jl_curr[j, l]
          eta_new <- rnorm(1, eta_old, hyperpar_VS$mh_sd_eta)
          
          log_prior_old <- dnorm(
            eta_old, mu_eta_curr[l], sigma_eta, log = TRUE
          )
          log_prior_new <- dnorm(
            eta_new, mu_eta_curr[l], sigma_eta, log = TRUE
          )
          
          log_like_old <- dbinom(
            gamma_jl_curr[j, l], 1, plogis(eta_old), log = TRUE
          )
          log_like_new <- dbinom(
            gamma_jl_curr[j, l], 1, plogis(eta_new), log = TRUE
          )
          
          log_acc <- (log_like_new + log_prior_new) -
            (log_like_old + log_prior_old)
          
          if (log(runif(1)) < log_acc) {
            eta_jl_curr[j, l] <- eta_new
          }
        }
      }
      
      # ========================================================
      # 4. Aggiornamento μ_{η,ℓ} (Gibbs – eq. 5)
      # ========================================================
      for (l in 1:n_cov) {
        var_post <- 1 / (k_curr / sigma_eta^2 + 1 / sigma_eta0^2)
        mean_post <- var_post * (
          sum(eta_jl_curr[, l]) / sigma_eta^2 +
            hyperpar_VS$mu_eta0 / sigma_eta0^2
        )
        mu_eta_curr[l] <- rnorm(1, mean_post, sqrt(var_post))
      }
      
      # ========================================================
      # Salvataggio
      # ========================================================
      if (g %in% saved_iters) {
        idx <- which(saved_iters == g)
        k_out[idx]      <- k_curr
        clu_out[idx, ]  <- clu_curr
        gamma_out[[idx]] <- gamma_jl_curr
        mu_eta_out[idx, ] <- mu_eta_curr
      }
      if(g %% 10 == 0) cat("Iterazione:", g, " | Cluster attivi:", k_curr, "\n")
    }
    
    return(list(
      k_out = k_out,
      clu_out = clu_out,
      gamma_out = gamma_out,
      mu_eta_out = mu_eta_out
    ))
  }
}
