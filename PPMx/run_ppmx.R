# --- run_model.R (Versione Finale Rifinita) ---

# -----------------------------------------------------------------
# 1. CARICAMENTO DATI E MODULI
# -----------------------------------------------------------------
if (!file.exists("data_generated.RData")) stop("Genera prima i dati!")
load("data_generated.RData") 

# Caricamento funzioni
source("PPMx_kernel_gaussiano.R")      
source("predictive_posterior_ppmx.R")  
source("cluster_salso.R")              

library(stats)
library(salso)
library(progress) # Assicuriamoci che sia caricata

# -----------------------------------------------------------------
# 2. CONFIGURAZIONE E RUN
# -----------------------------------------------------------------
hyperpar_P0_y <- list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)
hyperpar_P0_x <- list(
  x1 = list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1),
  x2 = list(mu0 = 0, kappa0 = 0.01, nu0 = 2, sig_square0 = 1)
)
hyperpar_Dir <- list(alpha = 1)

init <- list(
  k = 1L,
  clu = rep(1L, length(y)),
  nj = c(length(y)),
  tauj = c(1.0),
  muj = c(mean(y))
)

cat("Avvio campionatore PPMx...\n")
res <- PPMx_kernel_gaussiano(
  y = y, x = x,
  n_iter = 500, burn_in = 100, thin = 5,
  hyperpar_P0_y = hyperpar_P0_y,
  hyperpar_P0_x = hyperpar_P0_x,
  hyperpar_Dir = hyperpar_Dir,
  init = init
)

# -----------------------------------------------------------------
# 3. ANALISI PREDITTIVA
# -----------------------------------------------------------------
cat("Calcolo densità predittiva...\n")

# Definiamo un punto di test (es. X1 = 3, X2 = 0)
x_test <- c(3, 0) 
grid_y <- seq(min(y)-2, max(y)+2, length.out = 500)

dens_ppmx <- calcola_densita_predittiva_ppmx(
  x_new = x_test,
  grid_points_y = grid_y,
  risultati_mcmc = res,
  x_dati = x,
  hyperpar_P0_y = hyperpar_P0_y,
  hyperpar_P0_x = hyperpar_P0_x,
  hyperpar_Dir = hyperpar_Dir
)

# Plot densità predittiva
plot(grid_y, dens_ppmx, type = "l", lwd = 3, col = "royalblue",
     main = paste("Densità Predittiva per X =", paste(x_test, collapse=", ")),
     xlab = "y", ylab = "Densità")

# NOTA: Se vuoi vedere le medie vere, dovresti definirle manualmente qui 
# o salvarle nello script di generazione. 
# Esempio: abline(v = c(-3, -1, 0.5, 2, 3.5), col = "red", lty = 3)

# -----------------------------------------------------------------
# 4. PARTIZIONE OTTIMALE (SALSO)
# -----------------------------------------------------------------
cat("Ricerca partizione ottimale con SALSO...\n")
partizione_ottimale <- riassumi_cluster_salso(res, loss_type = "binder")

# Plot finale: Y vs X1 colorato per cluster ottimale
plot(x[,1], y, col = as.numeric(partizione_ottimale) + 1, pch = 20, 
     main = "Clustering Ottimale (SALSO) - Spazio Y vs X1",
     xlab = "X1 (Informativa)", ylab = "Y (Risposta)")

cat("Analisi completata con successo.\n")
