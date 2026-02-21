# --- mle_functions.R ---

#' Calcola le stime di Massima Verosimiglianza (MLE) per media e precisione
#' (parametri xi) per ciascuna covariata, assumendo una distribuzione Normale.
#' 
#' @details In una distribuzione Normale, la MLE della varianza è la varianza campionaria 
#'          biassata (divisa per n e non n-1). La precisione è il reciproco della varianza.
#'
#' @param x_dati Matrice (n x P) delle covariate (standardizzate).
#'
#' @return Una lista di liste. Ogni elemento contiene mu (media) e tau (precisione).
calcola_mle_xi <- function(x_dati) {
  # n = numero di osservazioni, n_covariate = numero di predittori
  n <- nrow(x_dati)
  n_covariate <- ncol(x_dati)
  
  # Inizializziamo una lista per contenere i parametri stimati per ogni colonna
  mle_list <- vector("list", n_covariate)
  
  # Iteriamo su ogni covariata (colonna)
  for (l in 1:n_covariate) {
    x_l <- x_dati[, l]
    
    # 1. Media MLE: Coincide con la media aritmetica campionaria
    mu_bar_l <- mean(x_l)
    
    # 2. Varianza MLE: Calcolata dividendo per 'n'. 
    # Nota: var() in R usa 'n-1', quindi calcoliamo manualmente per coerenza con la MLE teorica.
    sigma_sq_bar_l <- sum((x_l - mu_bar_l)^2) / n
    
    # Prevenzione di divisione per zero:
    # Se la varianza è zero (es. colonna costante), la precisione diventerebbe infinita.
    if (sigma_sq_bar_l <= 0) {
      sigma_sq_bar_l <- 1e-6 # Small epsilon per stabilità numerica
    }
    
    # 3. Precisione MLE: Definita come tau = 1 / sigma^2
    tau_bar_l <- 1 / sigma_sq_bar_l
    
    # Salviamo i risultati per la specifica covariata 'l'
    mle_list[[l]] <- list(mu = mu_bar_l, tau = tau_bar_l)
  }
  
  # Assegniamo i nomi originali delle colonne alla lista per facilitare il richiamo
  names(mle_list) <- colnames(x_dati)
  return(mle_list)
}

#' Calcola la log-densità del modello ausiliario q(x_i | xi_bar)
#'
#' @description Questa funzione valuta quanto è probabile un valore x_i 
#'              rispetto ai parametri stimati precedentemente.
#'
#' @param x_i Il valore della singola covariata da valutare.
#' @param l L'indice (numero) o il nome della colonna della covariata.
#' @param mle_xi_list L'output della funzione calcola_mle_xi.
#'
#' @return Il valore della log-densità Normale (log-likelihood puntuale).
calcola_log_densita_mle <- function(x_i, l, mle_xi_list) {
  
  # Recuperiamo i parametri stimati (mu e tau) per la covariata l
  mle_params <- mle_xi_list[[l]]
  
  # Calcoliamo la densità della distribuzione Normale in scala logaritmica.
  # Poiché dnorm richiede la deviazione standard (sd), convertiamo la precisione:
  # sd = sqrt(varianza) = sqrt(1 / tau) = 1 / sqrt(tau)
  log_dens <- dnorm(x_i, 
                    mean = mle_params$mu, 
                    sd = 1 / sqrt(mle_params$tau), 
                    log = TRUE)
  
  return(log_dens)
}