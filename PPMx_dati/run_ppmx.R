# --- run_model.R (Analisi Dati Reali Laghi con PPMx) ---
# Questo script coordina l'intera analisi: caricamento, campionamento MCMC 
# ad alta intensità e interpretazione dei risultati tramite SALSO.

# -----------------------------------------------------------------
# 1. CARICAMENTO DATI E MODULI
# -----------------------------------------------------------------
source("load_data.R")           # Esegue preprocessing e scaling
source("PPMx_kernel_gaussiano.R")
source("predictive_posterior_ppmx.R")
source("cluster_salso.R")

library(stats) 
library(salso) 
library(progress) 

# -----------------------------------------------------------------
# 2. DEFINIZIONE IPERPARAMETRI E RUN DEL MODELLO
# -----------------------------------------------------------------

n_covariate <- ncol(x) 
nomi_covariate <- colnames(x)

# Impostazione Empirica degli Iperparametri
# Usiamo i momenti dei dati (mean/var) per centrare la misura di base P0.
# Questo aiuta l'MCMC a esplorare zone dello spazio dei parametri supportate dai dati.
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

hyperpar_Dir <- list(alpha = 1) # Parametro di concentrazione standard

# Inizializzazione basata sulla varianza totale
init_var_y <- var(y)
init <- list(
  k = 1L, clu = rep(1L, length(y)), nj = c(length(y)),
  tauj = c(1.0 / ifelse(init_var_y == 0, 1e-6, init_var_y)),
  muj = c(mean(y))
)

# --- ESECUZIONE MCMC ---
# Utilizziamo un numero elevato di iterazioni (50.000) per garantire
# la convergenza in uno spazio multi-dimensionale complesso.
cat("Avvio campionatore PPMx ad alta risoluzione...\n")
res <- PPMx_kernel_gaussiano(
  y = y, x = x,
  n_iter = 5000,
  burn_in = 1000,
  thin = 1, # Salviamo ogni iterazione per una matrice di similarità accurata
  hyperpar_P0_y = hyperpar_P0_y,
  hyperpar_P0_x = hyperpar_P0_x,
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------------------------
# 3. DIAGNOSTICA E ANALISI
# -----------------------------------------------------------------

last_idx <- length(res$k_out)
clu_final <- res$clu_out[last_idx, ]

### 1) Traceplot del numero di cluster (k)
# Verifichiamo la stazionarietà della catena. Se il grafico "oscilla"
# stabilmente, la convergenza è verosimile.
plot(res$k_out, type = "l", lwd = 2, col = "darkblue",
     xlab = "Iterazione salvata", ylab = "Numero di cluster (k)",
     main = "Traceplot: Evoluzione del numero di cluster")

### 2) Analisi Dimensioni Cluster
# Visualizziamo come le masse dei dati si spostano tra i cluster principali.
sizes_mat <- matrix(0, nrow = last_idx, ncol = max(unlist(res$nj_out, use.names = FALSE)))
for(i in 1:last_idx){
  nj_sorted <- sort(res$nj_out[[i]], decreasing = TRUE)
  sizes_mat[i, 1:length(nj_sorted)] <- nj_sorted
}
matplot(sizes_mat, type = "l", lty = 1, lwd = 1.5,
        main = "Evoluzione dimensioni cluster (ordinate)",
        xlab = "Iterazione", ylab = "n_j")

# -----------------------------------------------------------------
# 4. STIMA PARTIZIONE OTTIMALE (SALSO)
# -----------------------------------------------------------------
# SALSO risolve il label switching trovando il clustering che sintetizza 
# meglio l'intera distribuzione a posteriori.

cat("Calcolo della partizione ottimale con SALSO...\n")
partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

# Visualizzazione: Y vs Covariate colorate per Cluster Ottimale
par(mfrow = c(3,3))
for (i in 1:min(n_covariate, 9)) {
  plot(x[,i], y, col = partizione_ottimale , pch = 16,
       main = paste0("Cluster vs ", nomi_covariate[i]),
       ylab = "Y (standard)", xlab = nomi_covariate[i])
}
par(mfrow = c(1, 1))

# -----------------------------------------------------------------
# 5. MATRICE DI SIMILARITÀ (Posterior Similarity Matrix)
# -----------------------------------------------------------------
# Mostra la probabilità a posteriori che ogni coppia di osservazioni 
# appartenga allo stesso cluster.
cat("Generazione Matrice di Similarità...\n")
similarity = salso::psm(res$clu_out)
heatmap(similarity, main = "Posterior Similarity Matrix", 
        col = gray.colors(100, start = 1, end = 0))

# -----------------------------------------------------------------
# 6. TABELLA RIASSUNTIVA (Profilazione dei Cluster)
# -----------------------------------------------------------------
# Fondamentale per dare un significato fisico ai cluster trovati.
# Mostra Media e Deviazione Standard per ogni covariata in ogni gruppo.

if (exists("partizione_ottimale") && exists("x")) {
  clusters_unici <- sort(unique(partizione_ottimale))
  n_cluster <- length(clusters_unici)
  n_cov <- ncol(x)
  
  tabella_riassuntiva <- matrix("", nrow = n_cov, ncol = n_cluster)
  rownames(tabella_riassuntiva) <- nomi_covariate
  colnames(tabella_riassuntiva) <- paste("Cluster", clusters_unici)
  
  for (j in 1:n_cluster) {
    idx <- which(partizione_ottimale == clusters_unici[j])
    x_subset <- x[idx, , drop = FALSE]
    
    for (i in 1:n_cov) {
      media <- mean(x_subset[, i], na.rm = TRUE)
      dev_std <- sd(x_subset[, i], na.rm = TRUE)
      if (is.na(dev_std)) dev_std <- 0.0
      
      tabella_riassuntiva[i, j] <- sprintf("%.3f (%.3f)", media, dev_std)
    }
  }
  
  cat("\n--- Profilazione Cluster: Media (Deviazione Standard) ---\n")
  print.noquote(tabella_riassuntiva)
}

