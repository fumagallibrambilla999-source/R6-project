# --- run_ppmx_logit_normal_FINAL.R ---

# ---------------------------------------------------------
# 1. CARICAMENTO SCRIPT E LIBRERIE
# ---------------------------------------------------------
source("load_data.R")      # Caricamento e preprocessing
source("PPMx_vero.R")      # Motore MCMC con Variable Selection
source("mle_functions.R")  # Funzioni ausiliarie per la verosimiglianza marginale

# salso: per trovare la partizione "ottimale" che minimizza la perdita (Variation of Information)
# pheatmap: per la visualizzazione delle matrici di adiacenza e inclusione
if (!requireNamespace("salso", quietly = TRUE)) install.packages("salso")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

# ---------------------------------------------------------
# 2. PREPARAZIONE INPUT E IPERPARAMETRI
# ---------------------------------------------------------
mle_xi_list <- calcola_mle_xi(x)
n_cov <- ncol(x)

# Iperparametri per la risposta Y (Prior Normal-Inverse-Gamma)
hyper_y <- list(mu0 = mean(y), kappa0 = 0.1, nu0 = 2, sig_square0 = var(y))

# Iperparametri per le covariate X (Liste indipendenti per ogni colonna)
hyper_x <- lapply(1:n_cov, function(l) {
  list(mu0 = mean(x[,l]), kappa0 = 10, nu0 = 3, sig_square0 = var(x[,l]))
})

# Parametri del processo PPMx: M controlla la coesione, alpha la nascita di nuovi cluster
hyper_dir <- list(M = 0.1, alpha = 0.1) 

# --- IPERPARAMETRI LOGIT-NORMAL (Variable Selection) ---
# Questi controllano la probabilità a priori che una variabile entri nel modello
hyper_vs_ln <- list(
  mu_eta0    = -1.5,   # Media negativa = prior "sparsa" (poche variabili incluse di base)
  sigma_eta0 = 1.0,    # Incertezza sulla media globale dell'inclusione
  sigma_eta  = 0.8,    # Variabilità dell'inclusione tra i diversi cluster
  mh_sd_eta  = 0.4     # Ampiezza del passo di Metropolis-Hastings (tuning del campionamento)
)

init <- list(clu = rep(1, length(y))) # Partenza con un singolo cluster globale

# ---------------------------------------------------------
# 3. ESECUZIONE MULTI-CHAIN (Robustezza Statistica)
# ---------------------------------------------------------
# Eseguiamo più catene per evitare di rimanere bloccati in ottimi locali
N_CHAINS <- 5 
N_ITER   <- 5000
N_BURN   <- 1000
N_THIN   <- 50 # Sfoltimento per ridurre l'autocorrelazione tra i campioni

chain_pips <- matrix(NA, nrow = N_CHAINS, ncol = n_cov)
all_res <- list()

cat("Avvio Multi-Chain MCMC (Modello Logit-Normal)...\n")
for(c in 1:N_CHAINS) {
  cat("  Esecuzione Catena", c, "...\n")
  res_c <- PPMx_VS_Quintana_LogitNormal(
    y, x, 
    n_iter = N_ITER, burn_in = N_BURN, thin = N_THIN, 
    hyperpar_P0_y = hyper_y, hyperpar_P0_x = hyper_x, 
    hyperpar_Dir = hyper_dir, hyperpar_VS = hyper_vs_ln, 
    mle_xi_list = mle_xi_list, init = init)
  
  all_res[[c]] <- res_c
  
  # Calcolo delle PIP (Posterior Inclusion Probabilities) medie per questa catena
  # Rappresentano la frequenza con cui ogni variabile ha gamma=1
  chain_pips[c,] <- colMeans(do.call(rbind, lapply(res_c$gamma_out, colMeans)))
}

# ---------------------------------------------------------
# 4. LOGICA DI CONSENSO
# ---------------------------------------------------------
# Media delle PIP tra tutte le catene per identificare le variabili stabili
pip_media <- colMeans(chain_pips)
pred_final <- which(pip_media > 0.5) # Regola del Modello Mediano (Barbieri & Berger)

cat("\n--- Risultati Variable Selection ---\n")
cat("Variabili selezionate:", paste(colnames(x)[pred_final], collapse=", "), "\n")

# ---------------------------------------------------------
# 5. DIAGNOSTICA E MATRICI DI SIMILARITÀ
# ---------------------------------------------------------
res_final <- all_res[[N_CHAINS]] # Analizziamo l'ultima catena come rappresentativa

# Stima della partizione ottimale (Point Estimate dei Cluster)
partizione_ottimale <- salso::salso(res_final$clu_out, loss = "VI")
psm <- salso::psm(res_final$clu_out) # Posterior Similarity Matrix

# 

par(mfrow = c(1, 2))
# Traceplot per verificare la stabilità del numero di cluster K
plot(res_final$k_out, type = "l", main = "Traceplot K", col="blue")
# Barplot delle PIP: evidenzia in rosso le variabili "vincenti"
colors <- rep("grey80", n_cov)
colors[pred_final] <- "firebrick"
barplot(pip_media, names.arg = colnames(x), las = 2, col = colors, 
        main = "PIP Consenso (Soglia 0.5)", ylim = c(0,1))
abline(h = 0.5, col = "blue", lty = 2)

# ---------------------------------------------------------
# 6. VISUALIZZAZIONE FISICA (Interpretazione dei Risultati)
# ---------------------------------------------------------

# A. BOXPLOT: Mostra come le covariate variano tra i cluster identificati
# Le variabili selezionate (titolo rosso) dovrebbero mostrare separazioni nette tra i boxplot
cat("\nGenerazione Boxplot per cluster...\n")
par(mfrow = c(ceiling(n_cov/3), 3), mar = c(4, 4, 3, 1))
for (l in 1:n_cov) {
  nome_var <- colnames(x)[l]
  colore_titolo <- if(l %in% pred_final) "firebrick" else "black"
  boxplot(x[, l] ~ partizione_ottimale, main = nome_var, 
          col = rainbow(length(unique(partizione_ottimale))), col.main = colore_titolo)
}

# B. SCATTERPLOT: Y vs X colorati per Cluster con rette di regressione locali
# Utile per vedere se i cluster catturano relazioni lineari diverse
cat("\nGenerazione Scatterplot Y vs Covariate...\n")
for (l in 1:n_cov) {
  plot(x[, l], y, col = rainbow(length(unique(partizione_ottimale)))[partizione_ottimale], 
       pch = 16, main = paste("Y vs", colnames(x)[l]))
  # Aggiunta regressione locale per ogni cluster significativo
  for(k in unique(partizione_ottimale)){
    idx_k <- which(partizione_ottimale == k)
    if(length(idx_k) > 5) abline(lm(y[idx_k] ~ x[idx_k, l]), col = rainbow(max(partizione_ottimale))[k])
  }
}

# ---------------------------------------------------------
# 7. HEATMAP DI SIMILARITÀ (SALSO) E PCA
# ---------------------------------------------------------
cat("\nGenerazione Heatmap di similarità e PCA...\n")

# A. HEATMAP DELLA MATRICE DI SIMILARITÀ POSTERIORE (PSM)
# La PSM contiene la probabilità (0-1) che l'osservazione i e j siano nello stesso cluster.
# Usiamo salso::psm per estrarla e pheatmap per visualizzarla.

# Ordiniamo le osservazioni per "partizione ottimale" così i blocchi sono visibili
ordine_ottimale <- order(partizione_ottimale)
psm_ordinata <- psm[ordine_ottimale, ordine_ottimale]

pheatmap::pheatmap(psm_ordinata, 
                   cluster_rows = FALSE,  # Ordine gestito manualmente dalla partizione
                   cluster_cols = FALSE, 
                   main = "Heatmap di Similarità Posteriore (Ordinata per Cluster)",
                   color = colorRampPalette(c("white", "#FEE08B", "#D73027"))(100),
                   show_rownames = FALSE, 
                   show_colnames = FALSE)



# B. ANALISI DELLE COMPONENTI PRINCIPALI (PCA)
# Proiettiamo le covariate in 2D per visualizzare la separazione dei cluster.

# Eseguiamo la PCA (i dati x sono già standardizzati)
pca_res <- prcomp(x, scale. = FALSE)
perc_var <- summary(pca_res)$importance[2, 1:2] * 100 # % varianza spiegata

# Creiamo un dataframe per il plot
df_pca <- data.frame(
  PC1 = pca_res$x[,1], 
  PC2 = pca_res$x[,2], 
  Cluster = as.factor(partizione_ottimale),
  Risposta_Y = y
)

# Plot PCA con ggplot2
library(ggplot2)
ggplot(df_pca, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(aes(size = Risposta_Y), alpha = 0.6) +
  scale_size_continuous(range = c(1, 5)) +
  theme_minimal() +
  labs(title = "Visualizzazione Cluster PPMx via PCA",
       subtitle = paste0("Varianza spiegata: PC1 (", round(perc_var[1], 1), 
                         "%) + PC2 (", round(perc_var[2], 1), "%)"),
       x = "Prima Componente Principale",
       y = "Seconda Componente Principale",
       size = "Valore Y (max90)") +
  scale_color_brewer(palette = "Set1")

