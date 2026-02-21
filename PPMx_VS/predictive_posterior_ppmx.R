# --- predictive_posterior_ppmx.R (Versione Marginalizzata) ---
source("PPMx_VS_gaussiano.R")

#' Calcola la densità predittiva per una nuova osservazione x_new
#' 
#' @param x_new Vettore delle covariate per il nuovo punto.
#' @param grid_y Griglia di valori su cui valutare la densità della risposta Y.
#' @param res Output della simulazione MCMC (contiene cluster, k, gamma, mu_z).
#' @param x_dati Matrice delle covariate originali.
#' @param y_dati Vettore della risposta originale.
#' @param h_P0_y Iperparametri della prior per la risposta Y.
#' @param h_P0_x Lista di iperparametri della prior per le covariate X.
#' @param h_Dir Parametri della distribuzione di Dirichlet (incluso alpha).
#' @param mle_list Lista dei parametri MLE (per il termine di normalizzazione).
#'
#' @return Un vettore della stessa lunghezza di grid_y con i valori di densità media.
calcola_densita_predittiva_ppmx_VS <- function(x_new, grid_y, res, x_dati, y_dati, 
                                               h_P0_y, h_P0_x, h_Dir, mle_list) {
  G <- length(res$k_out)         # Numero di iterazioni MCMC salvate
  dens_acc <- numeric(length(grid_y)) # Accumulatore per la densità media
  n_cov <- ncol(x_dati)
  
  # Ciclo sulle iterazioni della catena MCMC
  for (g in 1:G) {
    k_g <- res$k_out[g]          # Numero di cluster nell'iterazione g
    clu_g <- res$clu_out[g, ]    # Assegnazione dei cluster per ogni n
    gamma_g <- res$gamma_out[[g]] # Matrice Variable Selection (cluster x covariate)
    
    log_w <- numeric(k_g + 1)    # Log-pesi per l'allocazione del nuovo punto
    
    # --- 1. Calcolo pesi per i cluster ESISTENTI ---
    for (j in 1:k_g) {
      idx_j <- which(clu_g == j) # Indici delle osservazioni nel cluster j
      log_w[j] <- log(length(idx_j)) # Componente della "PPM" (cardinalità del cluster)
      
      # Influenza delle covariate sulla somiglianza (Similarity Function)
      for (l in 1:n_cov) {
        # Se la covariata l è "attiva" per il cluster j (gamma=1)
        if (gamma_g[j, l] == 1) {
          # Aggiungiamo il contributo della densità predittiva marginale (t di Student)
          # Sottraiamo la log-densità MLE per normalizzare rispetto al modello globale
          log_w[j] <- log_w[j] + (calcola_log_pred_dens_t(x_new[l], x_dati[idx_j, l], h_P0_x[[l]]) - 
                                    calcola_log_densita_mle(x_new[l], l, mle_list))
        }
      }
    }
    
    # --- 2. Calcolo peso per un NUOVO cluster ---
    # Si basa sul parametro di concentrazione alpha della Dirichlet Process
    log_w[k_g + 1] <- log(h_Dir$alpha)
    
    # Trasformazione log-sum-exp per ottenere probabilità normalizzate (evita l'underflow)
    w <- exp(log_w - max(log_w))
    prob_w <- w / sum(w)
    
    # --- 3. Calcolo densità su grid_y (Marginalizzazione) ---
    # Per ogni cluster esistente, calcoliamo la probabilità che y_new sia in grid_y
    for (j in 1:k_g) {
      idx_j <- which(clu_g == j)
      # Calcolo densità predittiva marginale (Y | Y_cluster_j)
      log_dens_y <- sapply(grid_y, function(val) calcola_log_pred_dens_t(val, y_dati[idx_j], h_P0_y))
      dens_acc <- dens_acc + exp(log_dens_y) * prob_w[j]
    }
    
    # Componente per il nuovo cluster (si basa solo sulla prior P0_y, dati vuoti)
    log_dens_y_new <- sapply(grid_y, function(val) calcola_log_pred_dens_t(val, numeric(0), h_P0_y))
    dens_acc <- dens_acc + exp(log_dens_y_new) * prob_w[k_g + 1]
  }
  
  # Restituiamo la densità media su tutte le iterazioni MCMC (Posterior Predictive)
  return(dens_acc / G)
}