# --- cluster_salso.R (Versione Semplificata) ---
riassumi_cluster_salso <- function(res, loss_type = "binder") {
  if (!requireNamespace("salso", quietly = TRUE)) stop("Install 'salso'")
  
  # Calcola la partizione ottimale basata sulla matrice di associazione
  # tra le osservazioni in tutte le iterazioni salvate
  partizione <- salso::salso(res$clu_out, loss = loss_type)
  
  return(partizione)
}