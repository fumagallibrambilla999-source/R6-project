# --- predictive_posterior_ppmx.R (Versione Marginalizzata) ---
source("PPMx_VS_gaussiano.R")

calcola_densita_predittiva_ppmx_VS <- function(x_new, grid_y, res, x_dati, y_dati, 
                                               h_P0_y, h_P0_x, h_Dir, mle_list) {
  G <- length(res$k_out)
  dens_acc <- numeric(length(grid_y))
  n_cov <- ncol(x_dati)
  
  for (g in 1:G) {
    k_g <- res$k_out[g]
    clu_g <- res$clu_out[g, ]
    gamma_g <- res$gamma_out[[g]]
    mu_z_g <- res$mu_z_out[[g]]
    
    log_w <- numeric(k_g + 1)
    # 1. Pesi per cluster esistenti
    for (j in 1:k_g) {
      idx_j <- which(clu_g == j)
      log_w[j] <- log(length(idx_j))
      for (l in 1:n_cov) {
        if (gamma_g[j, l] == 1) {
          log_w[j] <- log_w[j] + (calcola_log_pred_dens_t(x_new[l], x_dati[idx_j, l], h_P0_x[[l]]) - 
                                    calcola_log_densita_mle(x_new[l], l, mle_list))
        }
      }
    }
    # 2. Peso per nuovo cluster (neutrale come in allocazione)
    log_w[k_g + 1] <- log(h_Dir$alpha)
    
    # Normalizzazione pesi
    w <- exp(log_w - max(log_w))
    prob_w <- w / sum(w)
    
    # 3. Calcolo densità su grid_y
    for (j in 1:k_g) {
      idx_j <- which(clu_g == j)
      log_dens_y <- sapply(grid_y, function(val) calcola_log_pred_dens_t(val, y_dati[idx_j], h_P0_y))
      dens_acc <- dens_acc + exp(log_dens_y) * prob_w[j]
    }
    # Componente nuovo cluster
    log_dens_y_new <- sapply(grid_y, function(val) calcola_log_pred_dens_t(val, numeric(0), h_P0_y))
    dens_acc <- dens_acc + exp(log_dens_y_new) * prob_w[k_g + 1]
  }
  
  return(dens_acc / G)
}