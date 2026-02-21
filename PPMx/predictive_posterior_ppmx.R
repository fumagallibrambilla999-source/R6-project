# -----------------------------------------------------------------
# FUNZIONE HELPER: Verosimiglianza Marginale (Log-Space)
# -----------------------------------------------------------------
# Necessaria per calcolare la similarità tra la nuova osservazione x_new 
# e i cluster esistenti nello spazio delle covariate.

calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  alpha0  <- nu0 / 2
  beta0   <- 0.5 * (nu0 * sig_sq0)
  njj     <- length(z_j)
  
  if (njj == 0) {
    # Parametri a priori se il cluster è vuoto
    df_loc  <- 2 * alpha0
    mu_loc  <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  } else {
    # Aggiornamento posterior per il cluster j
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    ss_j      <- sum(z_j^2) - (sum_z^2) / njj 
    beta_loc  <- beta0 + 0.5 * ss_j + 0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    
    df_loc  <- 2 * alpha_loc
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) / (alpha_loc * kappa_loc)
  }
  
  if (sig_sq_loc <= 0) { return(-Inf) }
  sig_loc <- sqrt(sig_sq_loc)
  # Ritorna la densità della t-Student in scala logaritmica
  return(dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc))
}

# -----------------------------------------------------------------
# CALCOLO DENSITÀ PREDITTIVA POSTERIORE (PPMx)
# -----------------------------------------------------------------

#' Calcola la Densità Predittiva Posteriore per un modello PPMx
#' 
#' @description 
#' Implementa l'integrazione Monte Carlo: p(y_new | x_new, D) ≈ 1/G * Σ p(y_new | x_new, θ_g).
#' La probabilità che y_new appartenga a un cluster dipende dalla vicinanza di x_new alle covariate del cluster.

calcola_densita_predittiva_ppmx <- function(x_new,
                                            grid_points_y, 
                                            risultati_mcmc, 
                                            x_dati,
                                            hyperpar_P0_y, 
                                            hyperpar_P0_x, 
                                            hyperpar_Dir) {
  
  G <- length(risultati_mcmc$k_out)
  if (G == 0) stop("Nessuna iterazione salvata nell'output MCMC.")
  
  n_covariate <- length(x_new)
  alpha_val   <- hyperpar_Dir$alpha
  dens_finale_accumulator <- numeric(length(grid_points_y))
  
  # --- 1. Densità a priori per Y ---
  # Rappresenta il contributo di un potenziale "nuovo cluster"
  log_dens_y_prior <- sapply(grid_points_y, function(y_val) {
    calcola_log_pred_dens_t(y_val, numeric(0), hyperpar_P0_y)
  })
  dens_y_prior <- exp(log_dens_y_prior)
  
  # --- 2. Loop sulle iterazioni MCMC ---
  for (g in 1:G) {
    k_g    <- risultati_mcmc$k_out[g]
    nj_g   <- risultati_mcmc$nj_out[[g]]
    muj_g  <- risultati_mcmc$muj_out[[g]]
    tauj_g <- risultati_mcmc$tauj_out[[g]]
    clu_g  <- risultati_mcmc$clu_out[g, ]
    
    # --- FASE A: Calcolo pesi di allocazione condizionati a x_new ---
    # Calcoliamo quanto x_new è simile ai cluster esistenti nell'iterazione g
    log_w_alloc <- numeric(k_g + 1)
    
    for (j in 1:k_g) {
      log_prior_j <- log(nj_g[j])
      indici_j    <- which(clu_g == j)
      
      # Similarità basata sul prodotto delle densità marginali delle covariate
      log_sim_x_j <- 0.0
      for (l in 1:n_covariate) {
        log_sim_x_j <- log_sim_x_j + 
          calcola_log_pred_dens_t(x_new[l], x_dati[indici_j, l], hyperpar_P0_x[[l]])
      }
      log_w_alloc[j] <- log_prior_j + log_sim_x_j
    }
    
    # Peso per l'allocazione a un nuovo cluster
    log_sim_x_new <- 0.0
    for (l in 1:n_covariate) {
      log_sim_x_new <- log_sim_x_new + calcola_log_pred_dens_t(x_new[l], numeric(0), hyperpar_P0_x[[l]])
    }
    log_w_alloc[k_g + 1] <- log(alpha_val) + log_sim_x_new
    
    # Trasformazione da log-space a probabilità (Normalizzazione)
    w_alloc <- exp(log_w_alloc - max(log_w_alloc))
    prob_alloc <- w_alloc / sum(w_alloc)
    
    # --- FASE B: Calcolo densità della risposta Y ---
    # Miscela di densità Normali pesata dalle probabilità di allocazione appena calcolate
    dens_y_mat <- matrix(NA, nrow = length(grid_points_y), ncol = k_g)
    for (j in 1:k_g) {
      dens_y_mat[, j] <- dnorm(grid_points_y, mean = muj_g[j], sd = 1 / sqrt(tauj_g[j]))
    }
    
    # Media pesata: p(y | x_new) = Σ [ P(c_i=j|x_new) * N(y|μ_j, τ_j) ]
    dens_g_esistenti <- dens_y_mat %*% prob_alloc[1:k_g]
    dens_g_nuovo     <- dens_y_prior * prob_alloc[k_g + 1]
    
    dens_finale_accumulator <- dens_finale_accumulator + (dens_g_esistenti + dens_g_nuovo)
  }
  
  # Media Monte Carlo finale
  return(as.vector(dens_finale_accumulator / G))
}