#' Riassumere l'output MCMC in un'unica partizione (Modello PPMx)
#'
#' Utilizza il pacchetto 'salso' per trovare una stima puntuale (point estimate)
#' ottimale del clustering. Questo passaggio è cruciale nei modelli Bayesiani 
#' non parametrici per risolvere l'incertezza dovuta al "label switching".
#'
#' @param risultati_mcmc La lista di output della funzione `PPMx_kernel_gaussiano`.
#' @param loss_type La funzione di perdita da minimizzare. 
#'                  'binder' (Binder's loss) penalizza equamente coppie erroneamente 
#'                  unite o divise. 'VI' (Variation of Information) si basa 
#'                  sull'entropia ed è spesso più sensibile a strutture fini.
#'
#' @return Un vettore numerico (lungo n) che rappresenta la partizione 
#'         di clustering "ottimale" secondo la metrica scelta.

riassumi_cluster_salso <- function(risultati_mcmc, loss_type = "binder") {
  
  # --- 1. Verifica Dipendenze ---
  # Il pacchetto SALSO implementa algoritmi di ricerca euristica per ottimizzare
  # funzioni di perdita su spazi di partizioni.
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Il pacchetto 'salso' è necessario. Installalo con: install.packages('salso')")
  }
  
  # --- 2. Estrazione Matrice delle Allocazioni ---
  # Recuperiamo la matrice clu_out: dimensioni [G x n]
  # G = numero di iterazioni salvate, n = numero di osservazioni.
  allocations_matrix <- risultati_mcmc$clu_out
  
  if (is.null(allocations_matrix) || nrow(allocations_matrix) == 0) {
    stop("Errore: la matrice 'clu_out' è nulla o vuota. Verificare l'output del MCMC.")
  }
  
  # --- 3. Ottimizzazione della Partizione ---
  # Obiettivo: trovare la partizione c* che minimizza la perdita attesa a posteriori:
  # argmin_{c*} E[ L(c, c*) | y ] = sum_{c} L(c, c*) p(c | y)
  
  cat("Esecuzione di SALSO (metodo:", loss_type, ")...\n")
  
  # Utilizziamo tryCatch per gestire eventuali parametri incompatibili tra versioni
  optimal_partition <- tryCatch({
    # Nota: Nelle versioni recenti l'argomento principale è 'x'
    salso::salso(
      x = allocations_matrix, 
      loss = loss_type
    )
  }, error = function(e) {
    # Fallback nel caso in cui la versione installata usi ancora 'draws'
    message("Avviso: tentata chiamata con parametro 'draws' a causa di un errore.")
    salso::salso(
      draws = allocations_matrix,
      loss = loss_type
    )
  })
  
  cat("Fatto. Partizione ottimale identificata.\n")
  
  # --- 4. Output ---
  # Restituisce un oggetto di classe 'salso.estimate', che contiene anche 
  # informazioni sulla convergenza dell'ottimizzatore.
  return(optimal_partition)
}