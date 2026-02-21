#' Riassumere l'output MCMC in un'unica partizione
#'
#' @param risultati_mcmc La lista di output della funzione `DPMM_kernel_gaussiano`.
#' @param loss_type La funzione di perdita da minimizzare ('binder' o 'VI').
#'
#' @return Un vettore di clustering ottimale.

riassumi_cluster_salso <- function(risultati_mcmc, loss_type = "binder") {
  
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Il pacchetto 'salso' è necessario.")
  }
  
  # Estraiamo la matrice: deve essere una matrice dove ogni riga è un'osservazione 
  # e ogni colonna è un'iterazione MCMC (o viceversa, salso di solito vuole n_iter x n_obs)
  allocations_matrix <- risultati_mcmc$clu_out
  
  if (is.null(allocations_matrix) || nrow(allocations_matrix) == 0) {
    stop("La matrice 'clu_out' è vuota.")
  }
  
  cat("Inizio elaborazione SALSO...\n")
  
  # Salso richiede spesso che l'input sia esplicitamente una matrice.
  # Proviamo la chiamata più robusta specifica per la versione recente:
  optimal_partition <- tryCatch({
    # Nota: alcune versioni di salso usano 'x' invece di 'draws'
    salso::salso(
      x = allocations_matrix, 
      loss = loss_type
    )
  }, error = function(e) {
    # Fallback: se 'x' fallisce, prova con 'draws'
    salso::salso(
      draws = allocations_matrix,
      loss = loss_type
    )
  })
  
  cat("Fatto. Partizione ottimale trovata.\n")
  return(optimal_partition)
}