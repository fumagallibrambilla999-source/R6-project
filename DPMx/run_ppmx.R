# --- run_model.R (Modificato per PPMx) ---

# -----------------------------------------------
# 1. CARICAMENTO DATI E FUNZIONI
# -----------------------------------------------

# Carica i dati (assume che esista 'y')
load("data_generated.RData") 

# Carica le nostre funzioni custom
source("PPMx_kernal_gaussiano.R")     # !!! NUOVO: Il campionatore PPMx
source("predictive_posterior_ppmx.R") # (DA CORREGGERE)
source("cluster_salso.R")        # (OK)

# Carica le librerie necessarie
library(stats) # Per dnorm, dt, etc.
library(salso) # Per riassumere il clustering

# -----------------------------------------------
# 2. DEFINIZIONE IPERPARAMETRI E RUN
# -----------------------------------------------

# --- NUOVO: Iperparametri per Y e X ---
# Iperparametri per la risposta Y
hyperpar_P0_y <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)

# Iperparametri per le covariate X (una lista per ogni covariata)
# Iperparametri per x1
hyperpar_P0_x1 <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)
# Iperparametri per x2
hyperpar_P0_x2 <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)
# Lista che le contiene
hyperpar_P0_x <- list(hyperpar_P0_x1, hyperpar_P0_x2)

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

# --- NUOVO: Esegui il campionatore PPMx ---
res <- PPMx_kernel_gaussiano(
  y = y,
  x = x,                     # Argomento 'x' aggiunto
  n_iter = 1000,
  burn_in = 100,
  thin = 10,
  hyperpar_P0_y = hyperpar_P0_y, # Argomento aggiornato
  hyperpar_P0_x = hyperpar_P0_x, # Argomento 'x' aggiunto
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------
# 3. DIAGNOSTICA E ANALISI
# -----------------------------------------------

### ===============================
###   STATISTICHE PRELIMINARI
### ===============================
# (Questa parte funziona ancora come prima)

# Indice dell'ultima iterazione salvata
last_idx <- length(res$k_out)
# Allocazione "plug-in"
clu_final <- res$clu_out[last_idx, ]

cat("Sommario delle dimensioni dei cluster (SOLO ULTIMA ITERAZIONE):\n")
print(table(clu_final))


### 1) Numero di cluster (k) nel tempo
plot(res$k_out, type = "l", lwd = 2,
     xlab = "Iterazione salvata",
     ylab = "Numero di cluster (k)",
     main = "Traceplot del numero di cluster")


### 2) Dati colorati per clustering "plug-in" (ultima iterazione)
# NOTA: Ora potremmo plottare y vs x1, colorato per cluster!
plot(x[,1], y, col = clu_final, pch = 16,
     main = "Allocazione cluster (vs X1)",
     ylab = "Valore (Y)", xlab = "Covariata (X1)")


### 3) Istogramma per ogni cluster (ultima iterazione)
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
        xlab = "Cluster", ylab = "μ_j (per Y)",
        col = rainbow(length(cluster_means)))


### 5) Precisioni dei cluster (ultima iterazione)
cluster_tau <- res$tauj_out[[last_idx]]
barplot(cluster_tau,
        main = "Precisioni cluster τ_j (solo ultima iterazione)",
        xlab = "Cluster", ylab = "τ_j (per Y)",
        col = rainbow(length(cluster_tau)))


### 6) Evoluzione delle dimensioni dei cluster
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


### 7) Heatmap delle allocazioni (dati × iterazioni)
# (Questa parte funziona ancora come prima)
assign_mat <- res$clu_out
image(t(assign_mat), col = rainbow(max(clu_final)),
      main = "Heatmap allocazioni cluster",
      xlab = "Iterazione salvata", ylab = "Indice Dati")


### ==========================================================
###   8) DENSITÀ PREDITTIVA POSTERIORE
### ==========================================================
###
### !!! ATTENZIONE - PROSSIMO PASSO !!!
###
### Il codice sottostante è ora ERRATO.
### La funzione 'calcola_densita_predittiva' è stata scritta per il modello
### DPMM, che calcola p(y_nuovo | y_dati).
###
### Per il modello PPMx, dobbiamo calcolare p(y_nuovo | x_nuovo, dati),
### che richiede una funzione diversa.
###
### Per ora, questa sezione è disattivata.
###
# 
# cat("Calcolo della densità predittiva posteriore...\n")
# 
# xs <- seq(min(y) - 2*sd(y), 
#           max(y) + 2*sd(y), 
#           length.out = 500)
# 
# dens_predittiva <- calcola_densita_predittiva(
#   risultati_mcmc = res,
#   hyperpar_P0    = hyperpar_P0_y, # <-- USA L'IPERPARAMETRO DI Y
#   hyperpar_Dir   = hyperpar_Dir,
#   n_obs          = length(y),
#   grid_points    = xs
# )
# 
# hist(y, probability = TRUE, breaks = 40, col = "gray80",
#      main = "Densità Predittiva Posteriore (METODO VECCHIO - NON ESEGUITO)",
#      xlab = "y")
# lines(xs, dens_predittiva, lwd = 3, col = "blue")
# 


### ==========================================================
###   9) STIMA PARTIZIONE OTTIMALE (SALSO)
### ==========================================================
# (Questa parte funziona ancora come prima)

cat("Calcolo della partizione ottimale con SALSO...\n")

partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

n_cluster_ottimali <- length(unique(partizione_ottimale))
cat("Numero ottimale di cluster (SALSO):", n_cluster_ottimali, "\n")
cat("Dimensioni dei cluster ottimali:\n")
print(table(partizione_ottimale))

# Grafico della partizione ottimale
plot(x[,1], y, 
     col = partizione_ottimale, 
     pch = 20, 
     cex = 1.5,
     main = "Clustering Ottimale (stimato con SALSO)",
     xlab = "Covariata (X1)",
     ylab = "Dati (Y)")

cat("Analisi completata.\n")
