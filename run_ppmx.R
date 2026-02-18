# --- run_model.R (Modificato con stampe di diagnostica) ---

# -----------------------------------------------
# 1. CARICAMENTO DATI E FUNZIONI
# -----------------------------------------------
source("load_data.R") 
source("PPMx_kernel_gaussiano.R")
source("predictive_posterior_ppmx.R")
source("cluster_salso.R")

library(stats) 
library(salso) 
library(progress) 


# -----------------------------------------------
# 2. DEFINIZIONE IPERPARAMETRI E RUN
# -----------------------------------------------

n_covariate <- ncol(x) 
nomi_covariate <- colnames(x)
hyperpar_P0_y <- list(mu0 = mean(y), kappa0 = 0.01, nu0 = 2, sig_square0 = var(y))
hyperpar_P0_x <- vector("list", n_covariate) 

for (l in 1:n_covariate) {
  mean_x_l <- mean(x[, l])
  var_x_l  <- var(x[, l])
  hyperpar_P0_x[[l]] <- list(
    mu0 = mean_x_l, kappa0 = 0.01, nu0 = 2,
    sig_square0 = ifelse(var_x_l == 0, 1e-6, var_x_l)
  )
}
names(hyperpar_P0_x) <- nomi_covariate
hyperpar_Dir <- list(alpha = 1)
init_var_y <- var(y)
init <- list(
  k = 1L, clu = rep(1L, length(y)), nj = c(length(y)),
  tauj = c(1.0 / ifelse(init_var_y == 0, 1e-6, init_var_y)),
  muj = c(mean(y))
)

# --- Esegui il campionatore PPMx ---
res <- PPMx_kernel_gaussiano(
  y = y,
  x = x,
  n_iter = 50000,
  burn_in = 10000,
  thin = 1,
  hyperpar_P0_y = hyperpar_P0_y,
  hyperpar_P0_x = hyperpar_P0_x,
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------
# 3. DIAGNOSTICA E ANALISI
# -----------------------------------------------

par(mfrow = c(1, 1))
last_idx <- length(res$k_out)
clu_final <- res$clu_out[last_idx, ]
cat("Sommario delle dimensioni dei cluster (SOLO ULTIMA ITERAZIONE):\n")
print(table(clu_final))

### 1) Numero di cluster (k) nel tempo
plot(res$k_out, type = "l", lwd = 2,
     xlab = "Iterazione salvata",
     ylab = "Numero di cluster (k)",
     main = "Traceplot del numero di cluster")

### 2) Dati colorati per clustering "plug-in" (ultima iterazione)
par(mfrow = c(3, 3))
# Loop per plottare Y contro ciascuna delle 7 covariate
for (i in 1:7) {
  
  plot(x[,i], y, col = clu_final, pch = 16,
       main = paste0("Allocazione cluster (vs X", i, ")"),
       ylab = "Y (max90, standardizzata)", 
       xlab = paste0("Covariata (", nomi_covariate[i], ", standardizzata)")
  )
}
par(mfrow = c(1, 1))

### 3) Evoluzione delle dimensioni dei cluster
# (Questa parte funziona ancora come prima)
sizes_mat <- matrix(0, nrow = last_idx, ncol = max(unlist(res$nj_out, use.names = FALSE)))
for(i in 1:last_idx){
  nj <- res$nj_out[[i]]
  nj_sorted <- sort(nj, decreasing = TRUE)
  sizes_mat[i, 1:length(nj_sorted)] <- nj_sorted
}
matplot(sizes_mat, type = "l", lty = 1, lwd = 1.5,
        main = "Evoluzione dimensioni cluster (ordinate)",
        xlab = "Iterazione salvata", ylab = "Dimensione cluster")


###4) STIMA PARTIZIONE OTTIMALE (SALSO)
cat("Calcolo della partizione ottimale con SALSO\n")
partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

n_cluster_ottimali <- length(unique(partizione_ottimale))
cat("Numero ottimale di cluster (SALSO):", n_cluster_ottimali, "\n")
cat("Dimensioni dei cluster ottimali:\n")
print(table(partizione_ottimale))

# Grafici della partizione ottimale
# Loop per plottare Y contro ciascuna delle 7 covariate
par(mfrow = c(3,3))
for (i in 1:7) {
  
  plot(x[,i], y, col = partizione_ottimale , pch = 16,
       main = paste0("Allocazione cluster con SALSO (vs X", i, ")"),
       ylab = "Y (max90, standardizzata)", 
       xlab = paste0("Covariata (", nomi_covariate[i], ", standardizzata)")
  )
}
par(mfrow = c(1, 1))


###5) Similarity matrix
similarity= psm(res$clu_out)
heatmap(similarity)

### 6) TABELLA RIASSUNTIVA MEDIA/DEV.STD PER CLUSTER

if (exists("partizione_ottimale") && exists("x")) {
  
  # --- 1. Preparazione ---
  # Ottieni i cluster unici e il numero di covariate
  clusters_unici <- sort(unique(partizione_ottimale))
  n_cluster <- length(clusters_unici)
  n_cov <- ncol(x)
  
  
  # --- 2. Inizializzazione Tabella ---
  
  # Creiamo una matrice di caratteri per contenere "Media (Dev.Std)"
  tabella_riassuntiva <- matrix("", nrow = n_cov, ncol = n_cluster)
  
  # Imposta i nomi di righe e colonne
  rownames(tabella_riassuntiva) <- nomi_covariate
  colnames(tabella_riassuntiva) <- paste("Cluster", clusters_unici)
  
  # --- 3. Calcolo e Popolamento Tabella ---
  
  # Itera su ogni cluster (le colonne della nostra tabella)
  for (j in 1:n_cluster) {
    
    cl_id <- clusters_unici[j]
    
    # Trova gli indici (righe) delle osservazioni in questo cluster
    idx <- which(partizione_ottimale == cl_id)
    
    # Estrai il sottoinsieme di dati per questo cluster
    x_subset <- x[idx, , drop = FALSE]
    
    # Itera su ogni covariata (le righe della nostra tabella)
    for (i in 1:n_cov) {
      
      # Calcola media e deviazione standard
      media <- mean(x_subset[, i], na.rm = TRUE)
      dev_std <- sd(x_subset[, i], na.rm = TRUE)
      
      # Se la dev. std. è NA (es. cluster con 1 sola osservazione), impostala a 0
      if (is.na(dev_std)) {
        dev_std <- 0.0
      }
      
      # Formatta la stringa 
      testo_cella <- sprintf("%.3f (%.3f)", media, dev_std)
      
      # Inserisci nella tabella
      tabella_riassuntiva[i, j] <- testo_cella
    }
  }
  
  # --- 4. Stampa ---
  
  cat("\n--- Tabella Riassuntiva (Media (Deviazione Standard)) per Covariata e Cluster ---\n")
  # 'print.noquote' stampa la tabella senza virgolette attorno a ogni stringa
  print.noquote(tabella_riassuntiva)
} 
