#' Riassumere l'output MCMC in un'unica partizione (Dati Reali)
#'
#' Utilizza il pacchetto 'salso' per trovare una stima puntuale (point estimate)
#' ottimale del clustering dal modello PPMx, risolvendo il problema del "label switching".
#' In contesti di dati reali, questo passaggio permette di passare da migliaia di 
#' possibili partizioni a una singola struttura decisionale interpretabile.
#'
#' @param risultati_mcmc La lista di output della funzione `PPMx_kernel_gaussiano`.
#' @param loss_type La funzione di perdita da minimizzare. 
#'                  'binder' (più comune, bilancia falsi positivi/negativi) 
#'                  'VI' (Variation of Information, basata sulla teoria dell'informazione).
#'
#' @return Un vettore numerico (lungo n) che rappresenta la singola
#'          partizione di clustering "migliore" (stima puntuale).

riassumi_cluster_salso <- function(risultati_mcmc, loss_type = "binder") {
  
  # --- 1. Controllo Dipendenze ---
  if (!requireNamespace("salso", quietly = TRUE)) {
    stop("Il pacchetto 'salso' è necessario. Installalo con: install.packages('salso')")
  }
  
  # --- 2. Estrazione della Matrice delle Allocazioni ---
  # clu_out è la matrice G (iterazioni) x n (osservazioni reali).
  # Ogni riga rappresenta un "voto" del modello sulla struttura dei dati.
  allocations_matrix <- risultati_mcmc$clu_out
  
  if (is.null(allocations_matrix) || nrow(allocations_matrix) == 0) {
    stop("La matrice 'clu_out' nei risultati è vuota o mancante. Verificare il run MCMC.")
  }
  
  # --- 3. Ottimizzazione della Partizione (SALSO) ---
  # L'algoritmo minimizza la perdita attesa a posteriori. 
  # Questo è il "Gold Standard" per riassumere modelli Bayesiani non parametrici.
  
  cat("Esecuzione di SALSO per trovare la partizione ottimale...\n")
  cat("Metodo selezionato:", loss_type, "\n")
  
  # Utilizziamo l'argomento 'x' come richiesto dalle versioni recenti del pacchetto
  optimal_partition <- salso::salso(
    x = allocations_matrix, 
    loss = loss_type
  )
  
  cat("Elaborazione completata. Partizione ottimale trovata.\n")
  
  # --- 4. Restituzione Risultato ---
  # Il vettore restituito può essere usato direttamente per colorare i plot dei dati reali
  return(optimal_partition)
}