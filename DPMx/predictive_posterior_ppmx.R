# -----------------------------------------------------------------
# FUNZIONE HELPER (Copiata da PPMx_kernal_gaussiano.R)
# -----------------------------------------------------------------
# È necessaria per calcolare la similarità di x_new
# -----------------------------------------------------------------

calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  njj <- length(z_j)
  
  if (njj == 0) {
    df_loc  <- 2 * alpha0
    mu_loc  <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  } else {
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    ss_j <- sum(z_j^2) - (sum_z^2) / njj 
    beta_loc  <- beta0 + 
      0.5 * ss_j + 
      0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    df_loc  <- 2 * alpha_loc
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) / (alpha_loc * kappa_loc)
  }
  
  if (sig_sq_loc <= 0) { return(-Inf) }
  sig_loc <- sqrt(sig_sq_loc)
  log_dens <- dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc)
  return(log_dens)
}

# -----------------------------------------------------------------
# FUNZIONE PER IL CALCOLO DELLA DENSITÀ PREDITTIVA (PPMx)
# -----------------------------------------------------------------
# Implementa la formula di integrazione Monte Carlo per PPMx:
# p(y_nuovo | x_nuovo, Dati) = (1/G) * sum_{g=1 a G} [ p(y_nuovo | x_nuovo, stato_g) ]
# -----------------------------------------------------------------

#' Calcola la Densità Predittiva Posteriore per un PPMx
#'
#' @param x_new Vettore di covariate (lungo P) per la *nuova* osservazione
#'              per cui si vuole la predizione.
#' @param grid_points_y Vettore di punti 'y' su cui valutare la densità.
#' @param risultati_mcmc Output della funzione `PPMx_kernel_gaussiano`.
#' @param x_dati Matrice (n x P) delle covariate *originali* usate per il training.
#' @param hyperpar_P0_y Lista iperparametri per Y.
#' @param hyperpar_P0_x Lista di liste di iperparametri per X.
#' @param hyperpar_Dir Lista contenente 'alpha'.
#'
#' @return Un vettore (lungo quanto grid_points_y) con la densità predittiva.

calcola_densita_predittiva_ppmx <- function(x_new,
                                            grid_points_y, 
                                            risultati_mcmc, 
                                            x_dati,
                                            hyperpar_P0_y, 
                                            hyperpar_P0_x, 
                                            hyperpar_Dir) {
  
  G <- length(risultati_mcmc$k_out)
  if (G == 0) stop("Nessuna iterazione salvata.")
  
  n_covariate <- length(x_new)
  if (n_covariate != ncol(x_dati)) stop("x_new ha un numero errato di covariate.")
  
  alpha_val <- hyperpar_Dir$alpha
  
  # Accumulatore finale per la densità (somma su G)
  dens_finale_accumulator <- numeric(length(grid_points_y))
  
  # --- Calcola densità predittiva a priori per Y (costante) ---
  # p(y_nuovo | P0_y) -> t-Student
  # Usiamo 'sapply' perché calcola_log_pred_dens_t non è vettorizzata su z_i
  log_dens_y_prior <- sapply(grid_points_y, function(y_val) {
    calcola_log_pred_dens_t(y_val, numeric(0), hyperpar_P0_y)
  })
  dens_y_prior <- exp(log_dens_y_prior)
  
  
  # --- Loop su tutte le G iterazioni MCMC salvate ---
  for (g in 1:G) {
    # Estrai stato all'iterazione g
    k_g   <- risultati_mcmc$k_out[g]
    nj_g  <- risultati_mcmc$nj_out[[g]]
    muj_g <- risultati_mcmc$muj_out[[g]]
    tauj_g<- risultati_mcmc$tauj_out[[g]]
    clu_g <- risultati_mcmc$clu_out[g, ]
    
    # --- 1. Calcola i pesi di allocazione per x_new ---
    # Dobbiamo calcolare P(c_nuovo = j | x_nuovo, stato_g) per j=1...k_g, k_g+1
    log_w_alloc <- numeric(k_g + 1)
    
    # Cluster Esistenti (j = 1...k_g)
    for (j in 1:k_g) {
      log_prior_j <- log(nj_g[j])
      
      # Similarità di x_new con il cluster j
      indici_j <- which(clu_g == j) # Indici delle osservazioni nel cluster j
      
      log_similarita_x_j <- 0.0
      for (l in 1:n_covariate) {
        log_similarita_x_j <- log_similarita_x_j + 
          calcola_log_pred_dens_t(x_new[l], x_dati[indici_j, l], hyperpar_P0_x[[l]])
      }
      log_w_alloc[j] <- log_prior_j + log_similarita_x_j
    }
    
    # Nuovo Cluster (j = k_g+1)
    log_prior_new <- log(alpha_val)
    log_similarita_x_new <- 0.0
    for (l in 1:n_covariate) {
      log_similarita_x_new <- log_similarita_x_new +
        calcola_log_pred_dens_t(x_new[l], numeric(0), hyperpar_P0_x[[l]])
    }
    log_w_alloc[k_g + 1] <- log_prior_new + log_similarita_x_new
    
    # Normalizza i pesi di allocazione
    max_log_w <- max(log_w_alloc)
    w_alloc <- exp(log_w_alloc - max_log_w)
    prob_alloc <- w_alloc / sum(w_alloc)
    
    # --- 2. Calcola la densità di Y pesata ---
    # p(y | x_nuovo, stato_g) = sum_j [ P(c=j | x_nuovo, g) * p(y | c=j, g) ]
    
    # Matrice (grid_size x k_g) delle densità per Y
    dens_y_mat <- matrix(NA, nrow = length(grid_points_y), ncol = k_g)
    
    for (j in 1:k_g) {
      # Densità di Y se appartiene al cluster j
      # (usiamo i parametri campionati, che è più efficiente)
      dens_y_mat[, j] <- dnorm(grid_points_y, 
                               mean = muj_g[j], 
                               sd = 1 / sqrt(tauj_g[j]))
    }
    
    # Densità pesata per i cluster esistenti
    dens_g_esistenti <- dens_y_mat %*% prob_alloc[1:k_g]
    
    # Densità pesata per il nuovo cluster
    dens_g_nuovo <- dens_y_prior * prob_alloc[k_g + 1]
    
    # Densità totale per l'iterazione g
    dens_g <- dens_g_esistenti + dens_g_nuovo
    
    # Aggiungi all'accumulatore finale
    dens_finale_accumulator <- dens_finale_accumulator + dens_g
  }
  
  # --- 3. Restituisci la media ---
  densita_media <- dens_finale_accumulator / G
  return(as.vector(densita_media))
}
