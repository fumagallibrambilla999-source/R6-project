# --- run_ppmx_logit_normal_FINAL.R ---

# ---------------------------------------------------------
# 1. Caricamento Script e Dati
# ---------------------------------------------------------
source("load_data.R")            
source("PPMx_vero.R") # File contenente PPMx_VS_Quintana_LogitNormal
source("mle_functions.R")
if (!requireNamespace("salso", quietly = TRUE)) install.packages("salso")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

# ---------------------------------------------------------
# 2. Preparazione Input e Iperparametri Ottimizzati
# ---------------------------------------------------------
mle_xi_list <- calcola_mle_xi(x)
n_cov <- ncol(x)

# Setup Iperparametri risposta Y
hyper_y <- list(mu0 = mean(y), kappa0 = 0.1, nu0 = 2, sig_square0 = var(y))

# Iperparametri covariate X
hyper_x <- lapply(1:n_cov, function(l) {
  list(mu0 = mean(x[,l]), kappa0 = 10, nu0 = 3, sig_square0 = var(x[,l]))
})

# Coesione del processo PPMx
hyper_dir <- list(M = 0.1, alpha = 0.1) 

# --- IPERPARAMETRI LOGIT-NORMAL (VS) ---
# Basati sull'ottimizzazione per bilanciare Precision/Recall
hyper_vs_ln <- list(
  mu_eta0    = -1.5,   # Prior media della propensione d'inclusione (logit scale)
  sigma_eta0 = 1.0,    # Deviazione standard della prior per mu_eta
  sigma_eta  = 0.8,    # Deviazione standard intra-cluster (rigidità della selezione)
  mh_sd_eta  = 0.4     # Step del Metropolis-Hastings per l'aggiornamento di eta
)

# Inizializzazione
init <- list(clu = rep(1, length(y)))

# ---------------------------------------------------------
# 3. Esecuzione Multi-Chain (Consenso)
# ---------------------------------------------------------
N_CHAINS <- 5  
N_ITER   <- 15000
N_BURN   <- 8000
N_THIN   <- 50

chain_pips <- matrix(NA, nrow = N_CHAINS, ncol = n_cov)
all_res <- list()

cat("Avvio Multi-Chain MCMC (Modello Logit-Normal)...\n")
for(c in 1:N_CHAINS) {
  cat("  Esecuzione Catena", c, "...\n")
  res_c <- PPMx_VS_Quintana_LogitNormal(
    y, x, 
    n_iter = N_ITER, 
    burn_in = N_BURN, 
    thin = N_THIN, 
    hyperpar_P0_y = hyper_y, 
    hyperpar_P0_x = hyper_x, 
    hyperpar_Dir = hyper_dir, 
    hyperpar_VS = hyper_vs_ln, 
    mle_xi_list = mle_xi_list, 
    init = init)
  
  all_res[[c]] <- res_c
  
  # Estrazione PIP della catena corrente
  # Calcolata mediando le matrici gamma prodotte nelle iterazioni post burn-in
  chain_pips[c,] <- colMeans(do.call(rbind, lapply(res_c$gamma_out, colMeans)))
}


# ---------------------------------------------------------
# 4. LOGICA DI CONSENSO E VARIABLE SELECTION
# ---------------------------------------------------------
pip_media <- colMeans(chain_pips)
pred_final <- which(pip_media > 0.5)

cat("\n--- Risultati Variable Selection (Consenso Logit-Normal) ---\n")
cat("Variabili selezionate:", paste(colnames(x)[pred_final], collapse=", "), "\n")

# ---------------------------------------------------------
# 5. DIAGNOSTICA E REPORTING (Su ultima catena valida)
# ---------------------------------------------------------
res_final <- all_res[[N_CHAINS]]
partizione_ottimale <- salso::salso(res_final$clu_out, loss = "VI")
psm <- salso::psm(res_final$clu_out)

par(mfrow = c(1, 2))

# A. Traceplot K
plot(res_final$k_out, type = "l", main = "Traceplot Numero Cluster (K)", col="blue")
abline(h = mean(res_final$k_out), col="red", lty=2)

# B. Barplot PIP di Consenso
colors <- rep("grey80", n_cov)
colors[pred_final] <- "firebrick"
barplot(pip_media, names.arg = colnames(x), las = 2, col = colors,
        main = "PIP Consenso LN (Soglia 0.5)", ylim = c(0,1))
abline(h = 0.5, col = "blue", lty = 2)

# C. Heatmap PSM
pheatmap::pheatmap(psm, main = "Matrice di Similarità Posteriore",
                   color = colorRampPalette(c("white", "orange", "red"))(100))

# D. Analisi Inclusione per Cluster (Heatmap PIP)
visualizza_pip_cluster <- function(res, nomi_covariate) {
  G <- length(res$gamma_out)
  max_k <- max(res$k_out)
  pip_matrix <- matrix(0, max_k, length(nomi_covariate))
  counts <- rep(0, max_k)
  for(g in 1:G) {
    k_g <- res$k_out[g]
    pip_matrix[1:k_g, ] <- pip_matrix[1:k_g, ] + res$gamma_out[[g]]
    counts[1:k_g] <- counts[1:k_g] + 1
  }
  pip_final <- pip_matrix / counts
  pip_final[is.nan(pip_final)] <- 0
  rownames(pip_final) <- paste("Cluster", 1:max_k)
  colnames(pip_final) <- nomi_covariate
  pheatmap::pheatmap(pip_final[counts > (G*0.1), , drop=FALSE], 
                     display_numbers = TRUE, main="PIP per Cluster (LN)")
}

visualizza_pip_cluster(res_final, colnames(x))

# ---------------------------------------------------------
# 6. VISUALIZZAZIONE FISICA DEI CLUSTER
# ---------------------------------------------------------

# A. BOXPLOT
cat("\nGenerazione Boxplot per tutte le covariate...\n")
colori_cluster <- rainbow(length(unique(partizione_ottimale)))
par(mfrow = c(ceiling(n_cov/3), 3), mar = c(4, 4, 3, 1))

for (l in 1:n_cov) {
  nome_var <- colnames(x)[l]
  colore_titolo <- if(l %in% pred_final) "firebrick" else "black"
  
  boxplot(x[, l] ~ partizione_ottimale, 
          main = paste("Var:", nome_var),
          xlab = "Cl", ylab = "Val",
          col = colori_cluster, 
          col.main = colore_titolo,
          las = 1)
  abline(h = 0, lty = 2, col = "gray")
}

# B. SCATTERPLOT con Regressione locale
cat("\nGenerazione Scatterplot Y vs Covariate...\n")
par(mfrow = c(ceiling(n_cov/3), 3), mar = c(4, 4, 3, 1))

for (l in 1:n_cov) {
  plot(x[, l], y, 
       col = colori_cluster[partizione_ottimale], 
       pch = 16,
       main = paste("Y vs", colnames(x)[l]), 
       xlab = colnames(x)[l], ylab = "Y")
  
  for(k in unique(partizione_ottimale)){
    idx_k <- which(partizione_ottimale == k)
    if(length(idx_k) > 5){ 
      abline(lm(y[idx_k] ~ x[idx_k, l]), col = colori_cluster[k], lwd = 1.5)
    }
  }
}

par(mfrow = c(1, 1))
