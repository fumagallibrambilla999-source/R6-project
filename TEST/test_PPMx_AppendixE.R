# ==============================================================================
# TEST PPMx - VERSIONE BILANCIATA per Dataset Appendice E
# ==============================================================================
# PROBLEMA con versione ultra-sparse:
# - K collassa a 1
# - Tutte le PIP ~ 0.02 (troppo basse)
# - MH acceptance rate troppo alto (0.79)
#
# SOLUZIONE:
# - Prior meno estrema
# - MH step size più grande
# - Prior meno informativa sulle covariate (per dare più peso alla similarity)
# ==============================================================================

source("PPMx_vero.R")
source("mle_functions.R")

load("Dataset_Simulati_AppendixE_WithNoise.RData")

# ==============================================================================
# CONFIGURAZIONE
# ========================================================++======================
N_DATASET_DA_USARE <- 2
N_CHAINS <- 2
N_ITER   <- 5000    # Aumentato per migliore convergenza
N_BURN   <- 1000     # Burn-in più lungo
N_THIN   <- 10

stats_results <- list()

par(mfrow = c(1, 1))

# ==============================================================================
# ANALISI
# ==============================================================================

for(d in 1:N_DATASET_DA_USARE) {
  
  cat("\n")
  cat("========================================================================\n")
  cat(" Dataset", d, "/", N_DATASET_DA_USARE, "(BALANCED VERSION)\n")
  cat("========================================================================\n")
  
  # --------------------------------------------------------------------------
  # Carica dati
  # --------------------------------------------------------------------------
  data_sim <- DATA[[d]]
  x <- data_sim$X
  y <- data_sim$Y
  true_cols <- data_sim$selected.cols
  noise_cols <- data_sim$noise.cols
  true_groups <- data_sim$true_clusters
  
  # Standardizzazione
  x <- scale(x)
  y <- as.vector(scale(y))
  
  if(is.null(colnames(x))) colnames(x) <- paste0("X", 1:ncol(x))
  
  mle_xi <- calcola_mle_xi(x)
  init <- list(clu = rep(1, length(y)))
  n_cov <- ncol(x)
  n_obs <- length(y)
  
  cat("\n Dataset info:\n")
  cat("   N osservazioni:      ", n_obs, "\n")
  cat("   N covariate:         ", n_cov, "\n")
  cat("   Covariate vere:      ", paste(true_cols, collapse = ", "), "\n")
  cat("   Covariate noise:     ", paste(noise_cols, collapse = ", "), "\n")
  cat("   N gruppi veri:       ", length(unique(true_groups)), "\n")
  
  # --------------------------------------------------------------------------
  # IPERPARAMETRI BILANCIATI
  # --------------------------------------------------------------------------
  
  # Risposta: prior debole
  hyper_y <- list(
    mu0 = mean(y),
    kappa0 = 0.1,
    nu0 = 2,
    sig_square0 = var(y)
  )
  
  # CAMBIAMENTO CHIAVE: Prior DEBOLE sulle covariate
  # Questo permette alla similarity function di avere più influenza
  hyper_x <- lapply(1:n_cov, function(l) {
    list(
      mu0 = mean(x[, l]),
      kappa0 = 10,      # aumentando kappa0 mi permette covariate eterogenee in un cluster
      nu0 = 3,         # RIDOTTO da 5
      sig_square0 = var(x[, l])
    )
  })
  
  # Dirichlet: favorisce pochi cluster ma non troppo
  hyper_dir <- list(
    M = 0.1,           # RIDOTTO da 0.5 → permette più cluster
    alpha = 0.1
  )
  
  # PRIOR LOGIT-NORMAL BILANCIATA
  hp_ln <- list(
    mu_eta0 = -1.2,    # AUMENTATO da -3.5 → meno sparso (~8% inclusione invece di ~3%)
    sigma_eta0 = 1.0,  # AUMENTATO da 0.3 → più flessibilità
    sigma_eta = 0.8,   # AUMENTATO da 0.3
    mh_sd_eta = 0.4    # AUMENTATO da 0.2 → esplora meglio (target acceptance ~0.3-0.4)
  )
  
  cat("\n Iperparametri:\n")
  cat("   mu_eta0 (sparsità):  ", hp_ln$mu_eta0, 
      " (", round(plogis(hp_ln$mu_eta0) * 100, 1), "% prior inclusion)\n", sep = "")
  cat("   sigma_eta (cluster var):", hp_ln$sigma_eta, "\n")
  cat("   M (concentrazione):  ", hyper_dir$M, "\n")
  cat("   kappa0 (info prior): ", hyper_x[[1]]$kappa0, "\n")
  
  # --------------------------------------------------------------------------
  # MCMC
  # --------------------------------------------------------------------------
  
  cat("\n Esecuzione MCMC...\n")
  cat("   Iterazioni:", N_ITER, "\n")
  cat("   Burn-in:   ", N_BURN, "\n")
  cat("   Thinning:  ", N_THIN, "\n\n")
  
  res_ln <- PPMx_VS_Quintana_LogitNormal(
    y, x, N_ITER, N_BURN, N_THIN,
    hyper_y, hyper_x, hyper_dir, hp_ln,
    mle_xi, init
  )
  
  # --------------------------------------------------------------------------
  # PIP
  # --------------------------------------------------------------------------
  
  n_save <- length(res_ln$gamma_out)
  pip_accumulo <- matrix(0, nrow = n_obs, ncol = n_cov)
  
  for(g in 1:n_save) {
    cluster_assegnati <- res_ln$clu_out[g, ]
    gamma_iterazione <- res_ln$gamma_out[[g]]
    pip_accumulo <- pip_accumulo + gamma_iterazione[cluster_assegnati, ]
  }
  
  pip_final <- pip_accumulo / n_save
  pip_avg <- colMeans(pip_final)
  
  avg_k <- mean(res_ln$k_out)
  sd_k <- sd(res_ln$k_out)
  
  cat("\n Diagnostica MCMC:\n")
  cat("   K medio:  ", round(avg_k, 2), " (SD =", round(sd_k, 2), ")\n")
  cat("   K range:  [", min(res_ln$k_out), ",", max(res_ln$k_out), "]\n")
  
  # --------------------------------------------------------------------------
  # Selezione (soglia adattiva)
  # --------------------------------------------------------------------------
  
  # Prova soglia 0.5 (median probability model)
  THRESHOLD <- 0.5
  selected_vars <- which(pip_avg > THRESHOLD)
  
  # Se troppo poche o troppe, aggiusta
  if(length(selected_vars) == 0) {
    THRESHOLD <- 0.3
    selected_vars <- which(pip_avg > THRESHOLD)
    cat("\n   NOTA: Soglia ridotta a", THRESHOLD, "\n")
  }
  
  if(length(selected_vars) == 0) {
    cat("\n   FALLBACK: Seleziono top-4 variabili\n")
    selected_vars <- order(pip_avg, decreasing = TRUE)[1:4]
    THRESHOLD <- pip_avg[selected_vars[4]]
  }
  
  # --------------------------------------------------------------------------
  # Metriche
  # --------------------------------------------------------------------------
  
  tp <- length(intersect(selected_vars, true_cols))
  fp <- length(intersect(selected_vars, noise_cols))
  fn <- length(setdiff(true_cols, selected_vars))
  tn <- length(setdiff(noise_cols, selected_vars))
  
  precision <- if(length(selected_vars) > 0) tp / length(selected_vars) else 0
  recall <- tp / length(true_cols)
  specificity <- tn / length(noise_cols)
  f1_score <- if((precision + recall) > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    0
  }
  
  # --------------------------------------------------------------------------
  # Salva
  # --------------------------------------------------------------------------
  
  stats_results[[d]] <- list(
    prec = precision,
    rec = recall,
    spec = specificity,
    f1 = f1_score,
    pip = pip_avg,
    selected_vars = selected_vars,
    true_vars = true_cols,
    noise_vars = noise_cols,
    avg_k = avg_k,
    sd_k = sd_k,
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    threshold = THRESHOLD
  )
  
  # --------------------------------------------------------------------------
  # Visualizzazione Migliorata: Histogram delle PIP e Significatività
  # --------------------------------------------------------------------------
  
  # Definizione colori basata sui risultati (Matrice di confusione)
  # - Darkgreen: Variabile vera e selezionata (TP)
  # - Orange: Variabile noise ma selezionata (FP)
  # - Firebrick: Variabile vera ma NON selezionata (FN)
  # - Grey80: Variabile noise correttamente esclusa (TN)
  
  colors <- rep("grey80", n_cov) # Default: TN
  colors[true_cols] <- "firebrick" # FN potenziale
  
  if(length(selected_vars) > 0) {
    tp_vars <- intersect(selected_vars, true_cols)
    fp_vars <- intersect(selected_vars, noise_cols)
    colors[tp_vars] <- "darkgreen"
    colors[fp_vars] <- "orange"
  }
  
  # Layout del grafico
  par(mar = c(5, 5, 5, 10)) # Più spazio a destra per la legenda
  
  # Creazione del barplot (Frequenza di selezione = PIP)
  bp <- barplot(
    pip_avg,
    main = paste0("FREQUENZA SELEZIONE VARIABILI (PIP)\nDataset ", d, " - Selezione Bilanciata"),
    col = colors,
    names.arg = paste0("X", 1:n_cov),
    las = 1,
    ylim = c(0, 1.15), # Spazio extra in alto per le etichette
    ylab = "Posterior Inclusion Probability (PIP)",
    xlab = "Covariate",
    cex.axis = 0.9,
    border = "white"
  )
  
  # Linee di soglia
  abline(h = THRESHOLD, col = "blue", lty = 2, lwd = 2)
  text(x = max(bp), y = THRESHOLD + 0.03, "Soglia", col = "blue", pos = 4, cex = 0.8)
  
  # Aggiunta etichette sopra le barre (Percentuale e Significatività)
  for (j in 1:n_cov) {
    # Percentuale PIP
    text(bp[j], pip_avg[j] + 0.02, 
         labels = paste0(round(pip_avg[j]*100, 1), "%"), 
         cex = 0.7, pos = 3, font = 2)
    
    # Etichetta "SIGNIFICATIVA"
    if (j %in% selected_vars) {
      rect(bp[j]-0.4, -0.05, bp[j]+0.4, -0.02, col = "black", border = NA)
      text(bp[j], -0.1, labels = "SELECTED", cex = 0.6, font = 2, srt = 90, pos = 2)
    }
  }
  
  # Legenda esterna migliorata
  par(xpd = TRUE)
  legend(
    x = max(bp) + 1.5, y = 1,
    title = "Classificazione",
    legend = c("TP (Vera Selezionata)", "FP (Noise Selezionata)", 
               "FN (Vera Mancata)", "TN (Noise Esclusa)"),
    fill = c("darkgreen", "orange", "firebrick", "grey80"),
    bty = "n",
    cex = 0.8
  )
  par(xpd = FALSE)
  
  # --------------------------------------------------------------------------
  # Report Testuale Consolidato
  # --------------------------------------------------------------------------
  
  cat("\n--- CONCLUSIONI DATASET", d, "---\n")
  cat(sprintf("Variabili Significative (>%.0f%%): %s\n", 
              THRESHOLD*100, paste0("X", selected_vars, collapse = ", ")))
  cat(sprintf("Accuratezza Selezione: Precision %.2f | Recall %.2f | F1 %.2f\n", 
              precision, recall, f1_score))
  
  # Tabella di riepilogo rapida
  res_table <- data.frame(
    Var = paste0("X", 1:n_cov),
    Vera = ifelse(1:n_cov %in% true_cols, "SI", "no"),
    Freq_Inclusione = paste0(round(pip_avg * 100, 1), "%"),
    Esito = status_vec <- sapply(1:n_cov, function(j) {
      if(j %in% true_cols && j %in% selected_vars) "✓ Correct (TP)"
      else if(!(j %in% true_cols) && j %in% selected_vars) "⚠ False Pos (FP)"
      else if(j %in% true_cols && !(j %in% selected_vars)) "✖ Missed (FN)"
      else "✓ Correct (TN)"
    })
  )
  print(res_table)
  cat("----------------------------------------\n\n")
}

# ==============================================================================
# REPORT FINALE
# ==============================================================================

cat("\n")
cat("========================================================================\n")
cat(" REPORT FINALE (BALANCED VERSION)\n")
cat("========================================================================\n")

prec_mean <- mean(sapply(stats_results, function(x) x$prec))
rec_mean <- mean(sapply(stats_results, function(x) x$rec))
spec_mean <- mean(sapply(stats_results, function(x) x$spec))
f1_mean <- mean(sapply(stats_results, function(x) x$f1))
k_mean <- mean(sapply(stats_results, function(x) x$avg_k))

tp_tot <- sum(sapply(stats_results, function(x) x$tp))
fp_tot <- sum(sapply(stats_results, function(x) x$fp))
fn_tot <- sum(sapply(stats_results, function(x) x$fn))
tn_tot <- sum(sapply(stats_results, function(x) x$tn))

cat("\n METRICHE MEDIE:\n")
cat(" ----------------------------------------\n")
cat(" Precision:   ", round(prec_mean, 3), "\n")
cat(" Recall:      ", round(rec_mean, 3), "\n")
cat(" Specificity: ", round(spec_mean, 3), "\n")
cat(" F1-Score:    ", round(f1_mean, 3), "\n")
cat(" K medio:     ", round(k_mean, 2), "(vero = 3)\n")
cat(" ----------------------------------------\n")

cat("\n CONTEGGI AGGREGATI:\n")
cat(" ----------------------------------------\n")
cat(" TP totali:", tp_tot, "(su", N_DATASET_DA_USARE * 4, ")\n")
cat(" FP totali:", fp_tot, "(su", N_DATASET_DA_USARE * 3, ")\n")
cat(" FN totali:", fn_tot, "\n")
cat(" TN totali:", tn_tot, "\n")
cat(" ----------------------------------------\n")

prec_agg <- tp_tot / (tp_tot + fp_tot)
rec_agg <- tp_tot / (tp_tot + fn_tot)
spec_agg <- tn_tot / (tn_tot + fp_tot)

cat("\n AGGREGATE:\n")
cat(" Precision:   ", round(prec_agg, 3), "\n")
cat(" Recall:      ", round(rec_agg, 3), "\n")
cat(" Specificity: ", round(spec_agg, 3), "\n")

# Tabella
cat("\n DETTAGLIO:\n")
cat(" ", paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf(" %3s | %4s | %4s | %4s | %4s | %2s | %2s | %4s\n",
            "D#", "Prec", "Rec", "Spec", "F1", "TP", "FP", "K"))
cat(" ", paste(rep("-", 70), collapse = ""), "\n")

for(d in 1:N_DATASET_DA_USARE) {
  s <- stats_results[[d]]
  cat(sprintf(" %3d | %.2f | %.2f | %.2f | %.2f | %2d | %2d | %.1f\n",
              d, s$prec, s$rec, s$spec, s$f1,
              s$tp, s$fp, s$avg_k))
}

cat(" ", paste(rep("-", 70), collapse = ""), "\n")

# PIP medie
cat("\n PIP MEDIE:\n")
pip_matrix <- sapply(stats_results, function(x) x$pip)
pip_means <- rowMeans(pip_matrix)
pip_sds <- apply(pip_matrix, 1, sd)

for(j in 1:7) {
  is_true <- j <= 4
  cat(sprintf(" X%d (%s): %.3f (±%.3f)\n",
              j,
              ifelse(is_true, "VERA ", "NOISE"),
              pip_means[j],
              pip_sds[j]))
}

# ==============================================================================
# ANALISI DELLA CLUSTERIZZAZIONE E GRAFICI
# ==============================================================================

# 1. Calcolo accuratezza Cluster (confronto partizione finale vs reale)
# ------------------------------------------------------------------------------
# Prendiamo l'ultima iterazione per la clusterizzazione stimata
clu_final <- res_ln$clu_out[nrow(res_ln$clu_out), ]

# Per calcolare "giusti/sbagliati" in modo formale si usa la matrice di contingenza
tab_cluster <- table(true_groups, clu_final)

# Calcolo semplificato di osservazioni "coerenti":
# (Quante coppie di punti sono messe nello stesso cluster sia nel vero che nello stimato)
library(mclust) # Se non ce l'hai: install.packages("mclust")
ari_score <- adjustedRandIndex(true_groups, clu_final)

cat("\n ANALISI CLUSTERIZZAZIONE:\n")
cat(" ----------------------------------------\n")
cat(" Indice di Rand Corretto (ARI):", round(ari_score, 3), "\n")
cat(" (Nota: ARI = 1 significa partizioni identiche)\n")
cat(" Matrice di Confusione Cluster (Righe=Veri, Colonne=Stimati):\n")
print(tab_cluster)
cat(" ----------------------------------------\n")

# 2. Visualizzazione Grafica
# ------------------------------------------------------------------------------
par(mfrow = c(2, 2)) # Griglia 2x2 per i grafici

# Grafico A: Cluster Reali vs Y (Colori = Cluster Stimati)
plot(true_groups + rnorm(n_obs, 0, 0.05), y, 
     col = clu_final + 1, pch = 16,
     main = "Y per Cluster REALI\n(Colori = Cluster Stimati)",
     xlab = "Cluster Reali (jittered)", ylab = "Y (scaled)")
grid()

# Grafico B: Cluster Stimati vs Y (Colori = Cluster Reali)
plot(clu_final + rnorm(n_obs, 0, 0.05), y, 
     col = true_groups + 1, pch = 16,
     main = "Y per Cluster STIMATI\n(Colori = Cluster Reali)",
     xlab = "Cluster Stimati (jittered)", ylab = "Y (scaled)")
grid()

# Grafico C: Numero di Cluster (K) per iterazione
k_values <- res_ln$k_out
plot(k_values, type = "s", col = "darkblue", lwd = 1.5,
     main = "Evoluzione del numero di Cluster",
     xlab = "Iterazioni salvate", ylab = "Numero di Cluster (K)",
     yaxt = "n") # Rimuovo l'asse Y automatico

# Forza l'asse Y a mostrare solo numeri interi dal minimo al massimo trovato
axis(2, at = seq(min(k_values), max(k_values), by = 1), las = 1)

# Aggiungo la linea del valore vero
abline(h = length(unique(true_groups)), col = "red", lty = 2, lwd = 2)

legend("topright", legend = c("K stimato (Step)", "K vero"), 
       col = c("darkblue", "red"), lty = c(1, 2), bty = "n", cex = 0.8)
grid(ny = NA) # Griglia solo verticale per non sporcare i livelli interi

# Grafico D: Boxplot della distribuzione di K
boxplot(res_ln$k_out, col = "lightblue", main = "Distribuzione di K", 
        ylab = "K")

# Reset layout
par(mfrow = c(1, 1))

cat("\n========================================================================\n\n")

save(stats_results, pip_matrix, file = "Results_PPMx_BALANCED.RData")
cat("Salvato: Results_PPMx_BALANCED.RData\n\n")