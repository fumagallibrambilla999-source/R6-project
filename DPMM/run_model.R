# --- run_model.R ---
# Script principale per l'esecuzione del modello DPMM e l'analisi dei risultati.
# Gestisce l'intero workflow: caricamento, campionamento MCMC e diagnostica.

# -----------------------------------------------------------------
# 1. CARICAMENTO DATI E DIPENDENZE
# -----------------------------------------------------------------

# Caricamento del dataset generato precedentemente
source("generate_data.R")

# Caricamento dei moduli custom (funzioni definite negli altri file)
source("DPMM_gaussian_kernel.R")   # Campionatore Gibbs per DPMM
source("predictive_posterior.R")   # Funzione per densità Monte Carlo
source("cluster_salso.R")          # Integrazione con pacchetto SALSO

# Librerie necessarie
library(stats) 
library(salso) # Essenziale per la stima puntuale del clustering

# -----------------------------------------------------------------
# 2. CONFIGURAZIONE MODELLO E ESECUZIONE MCMC
# -----------------------------------------------------------------

# Definizione iperparametri P0 (Base Measure G0: Normal-Gamma)
# kappa0 piccolo = prior non informativa sulla media
# nu0 e sig_square0 definiscono la conoscenza a priori sulla varianza
hyperpar_P0 <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)

# Iperparametro alpha del Dirichlet Process (Concentrazione)
# Valori bassi favoriscono pochi cluster, valori alti ne favoriscono molti
hyperpar_Dir <- list(alpha = 1)

# Inizializzazione catena (Partiamo da un unico grande cluster)
init <- list(
  k = 1L,
  clu = rep(1L, length(y)),
  nj = c(length(y)),
  tauj = c(1.0),
  muj = c(mean(y))
)

# Esecuzione del campionatore Gibbs
# Parametri consigliati per produzione: n_iter > 20000, burn_in 5000
cat("Avvio campionatore MCMC...\n")
res <- DPMM_kernel_gaussiano(
  y = y,
  n_iter = 5000,   # Numero di iterazioni totali
  burn_in = 1000,  # Periodo di riscaldamento da scartare
  thin = 10,       # Diradamento per ridurre l'autocorrelazione
  hyperpar_P0 = hyperpar_P0,
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------------------------
# 3. DIAGNOSTICA DELLA CATENA MCMC
# -----------------------------------------------------------------

# Numero di campioni salvati dopo burn-in e thinning
last_idx <- length(res$k_out)

### A) Convergenza del numero di cluster (k)
# Fondamentale per verificare se la catena ha raggiunto la stazionarietà
plot(res$k_out, type = "l", lwd = 2, col = "darkgreen",
     xlab = "Iterazione (post-thinning)",
     ylab = "Numero di cluster (k)",
     main = "Traceplot di k (Convergenza)")

### B) Analisi delle dimensioni dei cluster nel tempo
# Questo plot aiuta a capire se la struttura dei cluster è stabile o fluttuante
sizes_mat <- matrix(0, nrow = last_idx, ncol = max(unlist(res$nj_out)))
for(i in 1:last_idx){
  nj_sorted <- sort(res$nj_out[[i]], decreasing = TRUE)
  sizes_mat[i, 1:length(nj_sorted)] <- nj_sorted
}

matplot(sizes_mat, type = "l", lty = 1, lwd = 1.2,
        main = "Stabilità delle dimensioni dei cluster",
        xlab = "Iterazione", ylab = "Numero osservazioni (n_j)")

# -----------------------------------------------------------------
# 4. STIME PUNTUALI E VISUALIZZAZIONE (Plug-in vs SALSO)
# -----------------------------------------------------------------

### C) Clustering "Plug-in" (Ultima iterazione)
# NOTA: Soggetto a label switching, usare solo per check rapidi
clu_final <- res$clu_out[last_idx, ]

### D) Clustering Ottimale (SALSO)
# Risolve il label switching minimizzando la funzione di perdita di Binder
cat("Ottimizzazione della partizione tramite SALSO...\n")
partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

# Confronto visivo
par(mfrow = c(2,1))
plot(y, col = clu_final, pch = 16, main = "Clustering: Ultima Iterazione (MCMC)")
plot(y, col = partizione_ottimale, pch = 16, main = "Clustering: Partizione Ottimale (SALSO)")
par(mfrow = c(1,1))

# -----------------------------------------------------------------
# 5. DENSITÀ PREDITTIVA POSTERIORE
# -----------------------------------------------------------------
# Rappresenta la probabilità di osservare un nuovo dato y* # integrando su tutta l'incertezza dei parametri.

cat("Calcolo densità predittiva posteriore...\n")
xs <- seq(min(y) - 2*sd(y), max(y) + 2*sd(y), length.out = 500)

dens_predittiva <- calcola_densita_predittiva(
  risultati_mcmc = res,
  hyperpar_P0    = hyperpar_P0,
  hyperpar_Dir   = hyperpar_Dir,
  n_obs          = length(y),
  grid_points    = xs
)

# Plot Finale di sintesi
hist(y, probability = TRUE, breaks = 45, col = "gray90", border = "white",
     main = "Fit Finale: Dati vs Densità Predittiva", xlab = "y")
lines(xs, dens_predittiva, lwd = 3, col = "royalblue")
legend("topright", legend = c("Densità Predittiva MCMC"), 
       col = "royalblue", lwd = 3, bty = "n")

cat("Analisi completata con successo.\n")

