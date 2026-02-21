# -----------------------------------------------------------------
# FUNZIONE HELPER: Verosimiglianza Marginale (Log-Space)
# -----------------------------------------------------------------
# Calcola la densità predittiva per un singolo valore z_i dato un cluster z_j.
# Viene usata per determinare la "similarità" tra x_nuovo e i cluster esistenti.

calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  # Estrazione iperparametri base measure
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  njj    <- length(z_j)
  
  if (njj == 0) {
    # Caso: Nuovo cluster (densità predittiva a priori)
    df_loc  <- 2 * alpha0
    mu_loc  <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  } else {
    # Caso: Cluster esistente (densità predittiva a posteriori)
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
  
  # Protezione numerica
  if (sig_sq_loc <= 0) { return(-Inf) }
  sig_loc <- sqrt(sig_sq_loc)
  
  # Restituisce la log-densità t-Student
  return(dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc))
}

# -----------------------------------------------------------------
# FUNZIONE PER IL CALCOLO DELLA DENSITÀ PREDITTIVA (PPMx)
# -----------------------------------------------------------------
# Questa funzione risponde alla domanda: "Dato un nuovo set di covariate x_new,
# quale è la distribuzione probabile della risposta y?"
# -----------------------------------------------------------------

calcola_densita_predittiva_ppmx <- function(x_new,
                                            grid_points_y, 
                                            risultati_mcmc, 
                                            x_dati,
                                            hyperpar_P0_y, 
                                            hyperpar_P0_x, 
                                            hyperpar_Dir) {
  
  # --- 1. Controlli e Inizializzazione ---
  G <- length(risultati_mcmc$k_out)
  if (G == 0) stop("Nessuna iterazione salvata nell'output MCMC.")
  
  n_covariate <- length(x_new)
  if (n_covariate != ncol(x_dati)) stop("Dimensioni di x_new non coerenti con x_dati.")
  
  alpha_val <- hyperpar_Dir$alpha
  dens_finale_accumulator <- numeric(length(grid_points_y))
  
  # Calcolo preventivo della densità a priori per Y (per il termine "nuovo cluster")
  log_dens_y_prior <- sapply(grid_points_y, function(y_val) {
    calcola_log_pred_dens_t(y_val, numeric(0), hyperpar_P0_y)
  })
  dens_y_prior <- exp(log_dens_y_prior)
  
  # --- 2. Media Monte Carlo sulle iterazioni MCMC ---
  # Per ogni stato salvato della catena, valutiamo la densità predittiva
  for (g in 1:G) {
    k_g    <- risultati_mcmc$k_out[g]
    nj_g   <- risultati_mcmc$nj_out[[g]]
    muj_g  <- risultati_mcmc$muj_out[[g]]
    tauj_g <- risultati_mcmc$tauj_out[[g]]
    clu_g  <- risultati_mcmc$clu_out[g, ]
    
    # --- FASE A: Calcolo Pesi di Allocazione (Condizionati a x_new) ---
    # Calcoliamo la probabilità che la nuova osservazione appartenga a ciascun cluster
    # basandoci sulla "vicinanza" di x_new alle X già osservate.
    log_w_alloc <- numeric(k_g + 1)
    
    for (j in 1:k_g) {
      log_prior_j <- log(nj_g[j])
      indici_j    <- which(clu_g == j) 
      
      # Similarità: prodotto delle densità predittive per ogni covariata l
      log_sim_x_j <- 0.0
      for (l in 1:n_covariate) {
        log_sim_x_j <- log_sim_x_j + 
          calcola_log_pred_dens_t(x_new[l], x_dati[indici_j, l], hyperpar_P0_x[[l]])
      }
      log_w_alloc[j] <- log_prior_j + log_sim_x_j
    }
    
    # Peso per l'allocazione a un cluster completamente nuovo
    log_sim_x_new <- 0.0
    for (l in 1:n_covariate) {
      log_sim_x_new <- log_sim_x_new + calcola_log_pred_dens_t(x_new[l], numeric(0), hyperpar_P0_x[[l]])
    }
    log_w_alloc[k_g + 1] <- log(alpha_val) + log_sim_x_new
    
    # Normalizzazione pesi (Log-Sum-Exp trick implicito)
    w_alloc <- exp(log_w_alloc - max(log_w_alloc))
    prob_alloc <- w_alloc / sum(w_alloc)
    
    # --- FASE B: Calcolo Densità della Risposta Y ---
    # La densità predittiva per questa iterazione è una mistura di normali
    # pesata dalle probabilità di allocazione calcolate sopra.
    
    dens_g_esistenti <- numeric(length(grid_points_y))
    for (j in 1:k_g) {
      # Densità Normale basata sui parametri (mu, tau) del cluster j all'iterazione g
      dens_y_j <- dnorm(grid_points_y, mean = muj_g[j], sd = 1 / sqrt(tauj_g[j]))
      dens_g_esistenti <- dens_g_esistenti + (prob_alloc[j] * dens_y_j)
    }
    
    # Termine per il nuovo cluster (pesato)
    dens_g_nuovo <- dens_y_prior * prob_alloc[k_g + 1]
    
    # Accumulo del risultato dell'iterazione g
    dens_finale_accumulator <- dens_finale_accumulator + (dens_g_esistenti + dens_g_nuovo)
  }
  
  # --- 3. Output Finale ---
  # Restituisce la media di tutte le densità valutate (Integrazione Monte Carlo)
  return(dens_finale_accumulator / G)
}