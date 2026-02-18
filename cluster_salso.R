#' Riassumere l'output MCMC in un'unica partizione
#'
#' Utilizza il pacchetto 'salso' per trovare una stima puntuale (point estimate)
#' ottimale del clustering, risolvendo il problema del "label switching".
#'
#' @param risultati_mcmc La lista di output della funzione `DPMM_kernel_gaussiano`.
#' @param loss_type La funzione di perdita da minimizzare. 'binder' (Binder's
#'                  loss) è un'ottima scelta standard. 'VI' (Variation of
#'                  Information) è un'altra opzione comune.
#'
#' @return Un vettore numerico (lungo n) che rappresenta la singola
#'         partizione di clustering "migliore".

riassumi_cluster_salso <- function(risultati_mcmc, loss_type = "binder") {
  
  # 1. Controlla che 'salso' sia caricato
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Il pacchetto 'salso' è necessario. Installalo con: install.packages('salso')")
  }
  
  # 2. Estrai la matrice delle allocazioni
  # Questa è la matrice G x n dove ogni riga è un campione di allocazione
  allocations_matrix <- risultati_mcmc$clu_out
  
  if (is.null(allocations_matrix) || nrow(allocations_matrix) == 0) {
    stop("La matrice 'clu_out' nei risultati è vuota o mancante.")
  }
  
  # 3. Esegui salso
  # La funzione salso::salso() (o salso::estimate_clustering)
  # calcola la partizione ottimale.
  # Trova il clustering c* che minimizza la perdita attesa a posteriori
  # E[L(c, c*) | dati], dove L è la funzione di perdita (es. 'binder').
  
  cat("Esecuzione di salso per trovare la partizione ottimale...\n")
  
  # L'argomento 'draws' si aspetta la matrice di allocazioni G x n
  optimal_partition <- salso::salso(
    draws = allocations_matrix,
    loss = loss_type
  )
  
  cat("Fatto. Partizione ottimale trovata.\n")
  
  # 4. Restituisci il vettore di clustering ottimale
  return(optimal_partition)
}