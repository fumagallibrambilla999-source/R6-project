# --- run_model.R ---

# -----------------------------------------------
# 1. CARICAMENTO DATI E FUNZIONI
# -----------------------------------------------

# Carica i dati (assume che esista 'y')
load("data_generated.RData")

# Carica le nostre funzioni custom
source("DPMM_gaussian_kernel.R")       # Il campionatore MCMC
source("predictive_posterior.R") # Per il grafico della densità (Punto 8)
source("cluster_salso.R")    # Per la partizione ottimale (Punto 9)

# Carica le librerie necessarie
library(stats) # Per dnorm, dt, etc.
library(salso) # Per riassumere il clustering

# -----------------------------------------------
# 2. DEFINIZIONE IPERPARAMETRI E RUN
# -----------------------------------------------

# Iperparametri per P0 (Normal-Gamma)
hyperpar_P0 <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)
# Iperparametro per il DP
hyperpar_Dir <- list(alpha = 1)

# Stato iniziale (un singolo cluster)
init <- list(
  k = 1L,
  clu = rep(1L, length(y)),
  nj = c(length(y)),
  tauj = c(1.0),
  muj = c(mean(y))
)

# Esegui il campionatore MCMC
# NOTA: 5000 iterazioni sono poche, per un run reale usare valori > 20000
res <- DPMM_kernel_gaussiano(
  y = y,
  n_iter = 5000,
  burn_in = 1000,
  thin = 10,
  hyperpar_P0 = hyperpar_P0,
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------
# 3. DIAGNOSTICA E ANALISI
# -----------------------------------------------

### ===============================
###   STATISTICHE PRELIMINARI
### ===============================

# Indice dell'ultima iterazione salvata
last_idx <- length(res$k_out)
# Allocazione "plug-in" (solo dall'ultima iterazione,
# utile per check veloci ma soffre di label switching)
clu_final <- res$clu_out[last_idx, ]

cat("Sommario delle dimensioni dei cluster (SOLO ULTIMA ITERAZIONE):\n")
print(table(clu_final))


### 1) Numero di cluster (k) nel tempo
plot(res$k_out, type = "l", lwd = 2,
     xlab = "Iterazione salvata",
     ylab = "Numero di cluster (k)",
     main = "Traceplot del numero di cluster")


### 2) Dati colorati per clustering "plug-in" (ultima iterazione)
plot(y, col = clu_final, pch = 16,
     main = "Allocazione cluster (solo ultima iterazione)",
     ylab = "Valore", xlab = "Indice")


### 3) Istogramma per ogni cluster (ultima iterazione)
# Utile solo se i cluster sono ben separati
par(mfrow = c(1, max(clu_final)))
for(k in sort(unique(clu_final))){
  hist(y[clu_final == k], breaks = 20, col = k,
       main = paste("Cluster", k), xlab = "y")
}
par(mfrow = c(1,1))


### 4) Medie dei cluster (ultima iterazione)
cluster_means <- res$muj_out[[last_idx]]
barplot(cluster_means,
        main = "Medie cluster (solo ultima iterazione)",
        xlab = "Cluster", ylab = "μ_j",
        col = rainbow(length(cluster_means)))


### 5) Precisioni dei cluster (ultima iterazione)
cluster_tau <- res$tauj_out[[last_idx]]
barplot(cluster_tau,
        main = "Precisioni cluster τ_j (solo ultima iterazione)",
        xlab = "Cluster", ylab = "τ_j",
        col = rainbow(length(cluster_tau)))


### 6) Evoluzione delle dimensioni dei cluster
# Questo grafico è corretto perché usa nj_out e non soffre
# di label switching (la dimensione [50, 50] è = a [50, 50])
sizes_mat <- matrix(0, nrow = last_idx, ncol = max(unlist(res$nj_out, use.names = FALSE)))
for(i in 1:last_idx){
  nj <- res$nj_out[[i]]
  # Ordiniamo per dimensione per stabilizzare il plot
  nj_sorted <- sort(nj, decreasing = TRUE)
  sizes_mat[i, 1:length(nj_sorted)] <- nj_sorted
}

matplot(sizes_mat, type = "l", lty = 1, lwd = 1.5,
        main = "Evoluzione dimensioni cluster (ordinate)",
        xlab = "Iterazione salvata", ylab = "Dimensione cluster")


### 7) Heatmap delle allocazioni (dati × iterazioni)
# Utile per vedere il label switching e la stabilità
assign_mat <- res$clu_out
image(t(assign_mat), col = rainbow(max(clu_final)),
      main = "Heatmap allocazioni cluster",
      xlab = "Iterazione salvata", ylab = "Indice Dati")


### ==========================================================
###   8) DENSITÀ PREDITTIVA POSTERIORE (METODO CORRETTO)
### ==========================================================
# Sostituisce il vecchio metodo "plug-in"
# Questo implementa la formula dell'immagine:
# p(y|D) = (1/G) * sum[ p(y | stato_g) ]

cat("Calcolo della densità predittiva posteriore (metodo Monte Carlo)...\n")

# Definisci la griglia di punti su cui plottare
xs <- seq(min(y) - 2*sd(y), 
          max(y) + 2*sd(y), 
          length.out = 500)

# Calcola la densità chiamando la nostra funzione
dens_predittiva <- calcola_densita_predittiva(
  risultati_mcmc = res,
  hyperpar_P0    = hyperpar_P0,
  hyperpar_Dir   = hyperpar_Dir,
  n_obs          = length(y),
  grid_points    = xs
)

# Plotta il risultato
hist(y, probability = TRUE, breaks = 40, col = "gray80",
     main = "Densità Predittiva Posteriore (Corretta)",
     xlab = "y")
lines(xs, dens_predittiva, lwd = 3, col = "blue")
legend("topright", "Densità Predittiva (integrata)", col="blue", lwd=3, bty="n")


### ==========================================================
###   9) STIMA PARTIZIONE OTTIMALE (SALSO)
### ==========================================================
# Aggiunto per trovare il "miglior" clustering singolo
# risolvendo il problema del label switching.

cat("Calcolo della partizione ottimale con SALSO...\n")

# Chiama la nostra funzione wrapper per salso
partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

# Statistiche sulla partizione ottimale
n_cluster_ottimali <- length(unique(partizione_ottimale))
cat("Numero ottimale di cluster (SALSO):", n_cluster_ottimali, "\n")
cat("Dimensioni dei cluster ottimali:\n")
print(table(partizione_ottimale))

# Grafico della partizione ottimale
# (Questo è un grafico 1D, usiamo jitter per vedere i punti)
plot(y, 
     y = jitter(rep(0, length(y)), amount = 0.5),
     col = partizione_ottimale, 
     pch = 20, 
     cex = 1.5,
     main = "Clustering Ottimale (stimato con SALSO)",
     xlab = "Valore di y",
     ylab = "", yaxt = "n") # Nasconde asse y

cat("Analisi completata.\n")