# --- generate_data_complex.R (Modificato per PPMx) ---
set.seed(48)

n_main <- 1000

# --- 1. Definiamo i parametri per 5 cluster ---

# Probabilità di appartenenza ai cluster
probs <- c(0.15, 0.25, 0.3, 0.2, 0.1)

# Parametri per la risposta Y
y_means <- c(-3, -1, 0.5, 2, 3.5)
y_sds   <- c(0.3, 0.6, 0.4, 0.5, 0.3)

# Parametri per le COVARIATE X (P=2)

# X1: Covariata INFORMATIVA (medie diverse)
x1_means <- c(-5, -2, 0, 3, 5) 
x1_sds   <- c(1, 1, 1, 1, 1)

# X2: Covariata NON INFORMATIVA (medie identiche = rumore)
x2_means <- c(0, 0, 0, 0, 0)
x2_sds   <- c(1, 1, 1, 1, 1)


# --- 2. Generiamo i dati principali ---

# Vettori per contenere i risultati
y_main <- numeric(n_main)
x1_main <- numeric(n_main)
x2_main <- numeric(n_main)
true_clusters <- numeric(n_main) # Teniamo traccia della "verità"

for(i in 1:n_main){
  # 1. Scegli un cluster
  u <- runif(1)
  cluster <- which(u <= cumsum(probs))[1]
  true_clusters[i] <- cluster
  
  # 2. Genera y, x1, e x2 da quel cluster
  y_main[i]  <- rnorm(1, mean = y_means[cluster], sd = y_sds[cluster])
  x1_main[i] <- rnorm(1, mean = x1_means[cluster], sd = x1_sds[cluster])
  x2_main[i] <- rnorm(1, mean = x2_means[cluster], sd = x2_sds[cluster])
}


# --- 3. Aggiungiamo qualche outlier ---
# (li aggiungiamo sia a y che a x)
n_out <- 10
y_out  <- rnorm(n_out, mean = 8, sd = 1)   # Outlier in Y
x1_out <- rnorm(n_out, mean = 10, sd = 1)  # Outlier in X1
x2_out <- rnorm(n_out, mean = 0, sd = 1)   # Non-outlier in X2

# --- 4. Combiniamo dati principali e outlier ---
y <- c(y_main, y_out)
x <- cbind(x1 = c(x1_main, x1_out), 
           x2 = c(x2_main, x2_out))

# Combiniamo anche i cluster (etichettiamo gli outlier come cluster 0)
true_clusters_final <- c(true_clusters, rep(0, n_out))

# --- 5. Salviamo i dati ---
# Ora salviamo SIA y CHE x
save(y, x, file = "data_generated.RData")

# --- 6. Quick plot per visualizzare la distribuzione ---
par(mfrow = c(1, 2))

# Istogramma di Y
hist(y, breaks = 50, col = "lightblue", 
     main = "Dati Sintetici (Y)", xlab = "y")

# Plot Y vs X1 (informativa), colorato per cluster VERO
plot(x[,1], y, col = true_clusters_final + 1, pch = 20, 
     main = "Y vs X1 (Colorato per Cluster Veri)", xlab = "x1 (informativa)")
legend("topleft", legend = c("Outlier", 1:5), col = 1:6, pch = 20, bty = "n")

par(mfrow = c(1, 1))
