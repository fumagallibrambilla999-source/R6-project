# --- mle_functions.R ---

#' Calcola le stime di Massima Verosimiglianza (MLE) per media e precisione
#' (parametri xi) per ciascuna covariata, assumendo una distribuzione Normale.
#'
#' @param x_dati Matrice (n x P) delle covariate (standardizzate).
#'
#' @return Una lista di liste. Ogni elemento 'l' contiene:
#'         - mu (media campionaria)
#'         - tau (precisione campionaria 1/sigma_sq_MLE)
calcola_mle_xi <- function(x_dati) {
  n <- nrow(x_dati)
  n_covariate <- ncol(x_dati)
  
  mle_list <- vector("list", n_covariate)
  
  for (l in 1:n_covariate) {
    x_l <- x_dati[, l]
    
    # 1. Media MLE
    mu_bar_l <- mean(x_l)
    
    # 2. Varianza MLE (non corretta, divisa per n)
    sigma_sq_bar_l <- sum((x_l - mu_bar_l)^2) / n
    
    # Prevenzione di divisione per zero
    if (sigma_sq_bar_l <= 0) {
      sigma_sq_bar_l <- 1e-6 # Valore molto piccolo
    }
    
    # 3. Precisione MLE
    tau_bar_l <- 1 / sigma_sq_bar_l
    
    mle_list[[l]] <- list(mu = mu_bar_l, tau = tau_bar_l)
  }
  
  names(mle_list) <- colnames(x_dati)
  return(mle_list)
}

#' Calcola la log-densità del modello ausiliario q(x_i | xi_bar)
#'
#' @param x_i Il valore della singola covariata da valutare.
#' @param l L'indice (colonna) della covariata.
#' @param mle_xi_list L'output di calcola_mle_xi.
#'
#' @return Il valore della log-densità Normale (scalare).
calcola_log_densita_mle <- function(x_i, l, mle_xi_list) {
  
  mle_params <- mle_xi_list[[l]]
  
  log_dens <- dnorm(x_i, 
                    mean = mle_params$mu, 
                    sd = 1 / sqrt(mle_params$tau), 
                    log = TRUE)
  return(log_dens)
}