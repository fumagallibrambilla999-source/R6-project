# -----------------------------------------------------------------
# FUNZIONE PER IL CALCOLO DELLA DENSITÀ PREDITTIVA POSTERIORE
# -----------------------------------------------------------------
# Implementa la formula di integrazione Monte Carlo per modelli DPMM:
# p(y_new | data) ≈ (1/G) * sum_{g=1}^G [ p(y_new | state_g) ]
# -----------------------------------------------------------------

#' Calcola la Densità Predittiva Posteriore
#'
#' @param risultati_mcmc Lista contenente l'output della catena (k_out, nj_out, muj_out, tauj_out).
#' @param hyperpar_P0 Lista iperparametri: mu0, kappa0, nu0, sig_square0.
#' @param hyperpar_Dir Lista contenente alpha (parametro di concentrazione).
#' @param n_obs Numero di osservazioni originali (usato per i pesi n/(n+alpha)).
#' @param grid_points Punti su cui valutare la densità (es. seq(...)).
#'
#' @return Vettore della densità stimata sui punti della griglia.

calcola_densita_predittiva <- function(risultati_mcmc, 
                                       hyperpar_P0, 
                                       hyperpar_Dir, 
                                       n_obs, 
                                       grid_points) {
  
  # --- 1. Parametri Generali ---
  G <- length(risultati_mcmc$k_out) # Numero di campioni MCMC salvati
  if (G == 0) stop("Nessuna iterazione salvata trovata.")
  
  alpha_val <- hyperpar_Dir$alpha
  denom <- n_obs + alpha_val # Il denominatore della regola predittiva del DP
  
  # --- 2. Termine "Nuovo Cluster" (Integrazione della Base Measure P0) ---
  # p(y_new | G0) segue una distribuzione t-Student derivante dalla
  # combinazione Normal-Inverse-Gamma della misura di base.
  
  mu0    <- hyperpar_P0$mu0
  kappa0 <- hyperpar_P0$kappa0
  alpha0 <- hyperpar_P0$nu0 / 2
  beta0  <- 0.5 * hyperpar_P0$nu0 * hyperpar_P0$sig_square0
  
  df0    <- 2 * alpha0
  # Fattore di scala per la t: sqrt( beta*(kappa+1) / (alpha*kappa) )
  sig0   <- sqrt((beta0 * (kappa0 + 1)) / (alpha0 * kappa0))
  
  # Densità t-Student scalata
  dens_t_prior <- dt((grid_points - mu0) / sig0, df = df0) * (1 / sig0)
  
  # Questo termine è pesato per alpha/(n+alpha) e mediato su G
  # Essendo costante rispetto alle iterazioni, lo calcoliamo una volta sola.
  dens_term2_final <- (alpha_val / denom) * dens_t_prior
  
  # --- 3. Termine "Cluster Esistenti" (Media sulle Misture Posteriore) ---
  # Per ogni iterazione g, calcoliamo la densità della mistura di Normali definita dai cluster.
  
  dens_term1_accumulator <- numeric(length(grid_points))
  
  for (g in 1:G) {
    k_g    <- risultati_mcmc$k_out[g]
    nj_g   <- risultati_mcmc$nj_out[[g]]
    muj_g  <- risultati_mcmc$muj_out[[g]]
    tauj_g <- risultati_mcmc$tauj_out[[g]]
    
    # Accumulatore locale per la mistura dell'iterazione g
    dens_g_inner_sum <- numeric(length(grid_points))
    
    if (k_g > 0) {
      for (j in 1:k_g) {
        # Deviazione standard = 1/sqrt(precisione)
        sigma_j <- 1 / sqrt(tauj_g[j])
        # Peso del cluster j = n_j / (n + alpha)
        weight_j <- nj_g[j] / denom
        
        # Aggiungiamo il contributo della componente normale j alla mistura
        dens_g_inner_sum <- dens_g_inner_sum + weight_j * dnorm(grid_points, mean = muj_g[j], sd = sigma_j)
      }
    }
    
    # Sommiamo per poi fare la media Monte Carlo (1/G)
    dens_term1_accumulator <- dens_term1_accumulator + dens_g_inner_sum
  }
  
  dens_term1_final <- dens_term1_accumulator / G
  
  # --- 4. Sintesi Finale ---
  # La densità predittiva totale è la somma delle due componenti
  return(dens_term1_final + dens_term2_final)
}