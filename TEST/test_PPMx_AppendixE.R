# ==============================================================================
# TEST PPMx - VERSIONE BILANCIATA per Dataset Appendice E
# ==============================================================================
# OBIETTIVO: Risolvere il problema del collasso della catena (K=1) e della
# scarsa selezione (PIP basse) osservato nelle versioni precedenti.
#
# STRATEGIA DI BILANCIAMENTO:
# 1. Prior logit-Normal meno punitiva per favorire l'inclusione iniziale.
# 2. Incremento dello step size Metropolis-Hastings per migliorare l'esplorazione.
# 3. Riduzione della precisione a priori sulle X per dare più "peso" ai dati
#    rispetto alla struttura di base.
# ==============================================================================

# Caricamento moduli necessari
source("PPMx_vero.R")      # Motore del modello
source("mle_functions.R")  # Calcolo MLE per baseline marginale
source("genera_dataset_simulati.R")

# Caricamento dei 10 dataset sintetici generati precedentemente
load("Dataset_Simulati_AppendixE_WithNoise.RData")

# ==============================================================================
# CONFIGURAZIONE MCMC
# ==============================================================================
N_DATASET_DA_USARE <- 2    # Numero di simulazioni da testare
N_ITER   <- 5000           # Numero totale di iterazioni per catena
N_BURN   <- 1000           # Fase di riscaldamento (scartata)
N_THIN   <- 10             # Diradamento per ridurre l'autocorrelazione

stats_results <- list()    # Lista per archiviare le statistiche finali
par(mfrow = c(1, 1))

# ==============================================================================
# ANALISI ITERATIVA SUI DATASET
# ==============================================================================

for(d in 1:N_DATASET_DA_USARE) {
  
  cat("\n")
  cat("========================================================================\n")
  cat(" Dataset", d, "/", N_DATASET_DA_USARE, "(BALANCED VERSION)\n")
  cat("========================================================================\n")
  
  # --- 1. Preparazione Dati ---
  data_sim <- DATA[[d]]
  x <- data_sim$X
  y <- data_sim$Y
  true_cols  <- data_sim$selected.cols  # Indici 1, 2, 3, 4
  noise_cols <- data_sim$noise.cols     # Indici 5, 6, 7
  true_groups <- data_sim$true_clusters
  
  # Standardizzazione: Cruciale per rendere confrontabili le distanze tra covariate
  x <- scale(x)
  y <- as.vector(scale(y))
  
  if(is.null(colnames(x))) colnames(x) <- paste0("X", 1:ncol(x))
  
  # Pre-calcolo delle densità MLE per la normalizzazione della similarità
  mle_xi <- calcola_mle_xi(x)
  
  # Inizializzazione: partiamo con tutti i dati in un unico cluster
  init <- list(clu = rep(1, length(y)))
  n_cov <- ncol(x)
  n_obs <- length(y)
  
  # --- 2. Definizione Iperparametri Bilanciati ---
  
  # Prior sulla risposta Y: prior debole centrata sulla media campionaria
  hyper_y <- list(
    mu0 = mean(y),
    kappa0 = 0.1,
    nu0 = 2,
    sig_square0 = var(y)
  )
  
  # Prior sulle covariate X:
  # L'aumento di kappa0 e la riduzione di nu0 permettono ai dati di 
  # influenzare maggiormente la funzione di similarità.
  hyper_x <- lapply(1:n_cov, function(l) {
    list(
      mu0 = mean(x[, l]),
      kappa0 = 10,     
      nu0 = 3,         
      sig_square0 = var(x[, l])
    )
  })
  
  # Parametro di concentrazione Dirichlet (M): 
  # Un valore più basso facilita la suddivisione in più gruppi.
  hyper_dir <- list(
    M = 0.1, 
    alpha = 0.1
  )
  
  # Parametri Variable Selection (logit-Normal):
  # mu_eta0 = -1.2 alza la probabilità di inclusione a priori (~30%)
  # sigma_eta e sigma_eta0 aumentati per permettere eterogeneità tra cluster
  # mh_sd_eta aumentato per evitare che z_jl rimanga bloccata
  hp_ln <- list(
    mu_eta0 = -1.2,    
    sigma_eta0 = 1.0,  
    sigma_eta = 0.8,   
    mh_sd_eta = 0.4    
  )
  
  # --- 3. Esecuzione del Campionatore Gibbs (PPMx-VS) ---
  cat("\n Esecuzione MCMC...\n")
  res_ln <- PPMx_VS_Quintana_LogitNormal(
    y, x, N_ITER, N_BURN, N_THIN,
    hyper_y, hyper_x, hyper_dir, hp_ln,
    mle_xi, init
  )
  
  # --- 4. Calcolo PIP (Posterior Inclusion Probabilities) ---
  # La PIP rappresenta la probabilità che una variabile sia utile al clustering.
  n_save <- length(res_ln$gamma_out)
  pip_accumulo <- matrix(0, nrow = n_obs, ncol = n_cov)
  
  for(g in 1:n_save) {
    cluster_assegnati <- res_ln$clu_out[g, ]
    gamma_iterazione <- res_ln$gamma_out[[g]]
    # Sommiamo l'indicatore gamma relativo al cluster di appartenenza di ogni punto i
    pip_accumulo <- pip_accumulo + gamma_iterazione[cluster_assegnati, ]
  }
  
  pip_final <- pip_accumulo / n_save
  pip_avg <- colMeans(pip_final) # Media globale per ogni variabile
  
  avg_k <- mean(res_ln$k_out)
  
  # --- 5. Selezione Variabili (Soglia Adattiva) ---
  # Implementazione del "Median Probability Model" con meccanismi di sicurezza
  THRESHOLD <- 0.5
  selected_vars <- which(pip_avg > THRESHOLD)
  
  if(length(selected_vars) == 0) {
    THRESHOLD <- 0.3 # Rilassamento soglia se il segnale è debole
    selected_vars <- which(pip_avg > THRESHOLD)
  }
  
  if(length(selected_vars) == 0) {
    # Fallback estremo: prendiamo le migliori 4 (numero atteso di variabili vere)
    selected_vars <- order(pip_avg, decreasing = TRUE)[1:4]
    THRESHOLD <- pip_avg[selected_vars[4]]
  }
  
  # --- 6. Calcolo Metriche di Accuratezza ---
  tp <- length(intersect(selected_vars, true_cols))  # True Positives
  fp <- length(intersect(selected_vars, noise_cols)) # False Positives
  fn <- length(setdiff(true_cols, selected_vars))    # False Negatives
  tn <- length(setdiff(noise_cols, selected_vars))   # True Negatives
  
  precision <- if(length(selected_vars) > 0) tp / length(selected_vars) else 0
  recall <- tp / length(true_cols)
  f1_score <- if((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0
  
  # Archiviazione risultati
  stats_results[[d]] <- list(
    prec = precision, rec = recall, f1 = f1_score, pip = pip_avg,
    selected_vars = selected_vars, true_vars = true_cols,
    avg_k = avg_k, tp = tp, fp = fp, fn = fn, tn = tn, threshold = THRESHOLD
  )
  
  # --- 7. Visualizzazione Risultati (Barplot PIP) ---
  # Definizione colori: Verde=Signal/Correct, Grey=Noise/Correct, Orange/Red=Errors
  colors <- rep("grey80", n_cov) 
  colors[true_cols] <- "firebrick" # FN potenziale
  if(length(selected_vars) > 0) {
    colors[intersect(selected_vars, true_cols)] <- "darkgreen"
    colors[intersect(selected_vars, noise_cols)] <- "orange"
  }
  
  par(mar = c(5, 5, 5, 10)) # Margine largo per legenda
  bp <- barplot(pip_avg, main = paste0("FREQUENZA SELEZIONE VARIABILI (PIP)\nDataset ", d),
                col = colors, names.arg = paste0("X", 1:n_cov), las = 1,
                ylim = c(0, 1.15), ylab = "PIP", xlab = "Covariate", border = "white")
  
  abline(h = THRESHOLD, col = "blue", lty = 2, lwd = 2)
  text(x = max(bp), y = THRESHOLD + 0.03, "Soglia", col = "blue", pos = 4, cex = 0.8)
  
  # Etichette "SELECTED" per variabili sopra soglia
  for (j in 1:n_cov) {
    text(bp[j], pip_avg[j] + 0.02, labels = paste0(round(pip_avg[j]*100, 1), "%"), 
         cex = 0.7, pos = 3, font = 2)
    if (j %in% selected_vars) {
      text(bp[j], -0.1, labels = "SELECTED", cex = 0.6, font = 2, srt = 90, pos = 2)
    }
  }
  
  # Legenda
  par(xpd = TRUE)
  legend(x = max(bp) + 1.5, y = 1, title = "Classificazione",
         legend = c("TP (Vera Selezionata)", "FP (Noise Selezionata)", 
                    "FN (Vera Mancata)", "TN (Noise Esclusa)"),
         fill = c("darkgreen", "orange", "firebrick", "grey80"), bty = "n", cex = 0.8)
  par(xpd = FALSE)
  
  # Stampa report sintetico
  cat(sprintf("\nDataset %d - Precision: %.2f | Recall: %.2f | K medio: %.1f\n", d, precision, recall, avg_k))
}

# ==============================================================================
# ANALISI DELLA CLUSTERIZZAZIONE E GRAFICI FINALI
# ==============================================================================

# Analisi dell'ultimo dataset processato
clu_final <- res_ln$clu_out[nrow(res_ln$clu_out), ]



# Calcolo dell'Indice di Rand Corretto (ARI)
library(mclust) 
ari_score <- adjustedRandIndex(true_groups, clu_final)

cat("\nANALISI CLUSTERIZZAZIONE:\n")
cat("----------------------------------------\n")
cat("Indice di Rand Corretto (ARI):", round(ari_score, 3), "\n")
cat("----------------------------------------\n")

# Layout Diagnostico 2x2
par(mfrow = c(2, 2))

# A & B: Confronto Jittered tra Vero e Stimato
plot(true_groups + rnorm(n_obs, 0, 0.05), y, col = clu_final + 1, pch = 16,
     main = "Veri vs Y (Colori=Stimati)", xlab = "Cluster Reali", ylab = "Y")
plot(clu_final + rnorm(n_obs, 0, 0.05), y, col = true_groups + 1, pch = 16,
     main = "Stimati vs Y (Colori=Reali)", xlab = "Cluster Stimati", ylab = "Y")

# C: Evoluzione di K (Verifica convergenza stazionaria)
k_values <- res_ln$k_out
plot(k_values, type = "s", col = "darkblue", lwd = 1.5,
     main = "Evoluzione del numero di Cluster", xlab = "MCMC Iterations", ylab = "K", yaxt = "n")
axis(2, at = seq(min(k_values), max(k_values), by = 1), las = 1)
abline(h = length(unique(true_groups)), col = "red", lty = 2) # Target K=3

# D: Distribuzione a posteriori del numero di cluster
boxplot(res_ln$k_out, col = "lightblue", main = "Posterior Distribution of K")

# Reset layout e salvataggio
par(mfrow = c(1, 1))
save(stats_results, file = "Results_PPMx_BALANCED.RData")
cat("\nAnalisi completata. Risultati salvati in 'Results_PPMx_BALANCED.RData'\n")

