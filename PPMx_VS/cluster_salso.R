# --- cluster_salso.R (Versione Ottimizzata per PPMx-VS) ---

#' Riassume l'output MCMC in un'unica partizione ottimale
#'
#' @description
#' Questa funzione risolve il problema del "label switching" trovando la partizione
#' che minimizza la perdita attesa (Expected Loss) rispetto alla distribuzione 
#' a posteriori del clustering. È essenziale per interpretare i risultati di 
#' modelli non parametrici come il PPMx.
#'
#' @param res La lista di output restituita dal campionatore PPMx.
#' @param loss_type Il tipo di funzione di perdita: 
#'        - "binder": Bilancia il rischio di accorpare cluster diversi e separare osservazioni simili.
#'        - "VI": Variation of Information, basata sulla teoria dell'informazione (spesso più granulare).
#'
#' @return Un vettore di interi che rappresenta l'assegnazione finale dei cluster.

riassumi_cluster_salso <- function(res, loss_type = "binder") {
  
  # 1. Verifica la presenza del pacchetto SALSO (Specialized Additive Loss Sampling Approximation)
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Il pacchetto 'salso' è necessario per sintetizzare i risultati. Installalo con install.packages('salso')")
  }
  
  # 2. Controllo coerenza input
  if (is.null(res$clu_out)) {
    stop("L'oggetto fornito non contiene la matrice delle allocazioni 'clu_out'.")
  }
  
  cat("Inizio sintesi della partizione ottimale (Loss:", loss_type, ")... \n")
  
  # 3. Calcolo della partizione ottimale
  # salso() analizza la matrice G x N e trova la configurazione 'c' 
  # che meglio rappresenta la Posterior Similarity Matrix (PSM).
  
  
  
  partizione <- salso::salso(
    x = res$clu_out, 
    loss = loss_type, 
    maxIter = 100 # Numero di iterazioni per l'algoritmo di ricerca locale
  )
  
  # 4. Diagnostica rapida
  n_cluster <- length(unique(partizione))
  cat("Fatto. Numero di cluster identificati nella partizione ottimale:", n_cluster, "\n")
  
  return(as.integer(partizione))
}