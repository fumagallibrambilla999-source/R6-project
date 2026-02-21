# --- generate_data_complex.R (Modificato per PPMx) ---
# Script per la generazione di dati sintetici con struttura latente legata a covariate.
# Ideale per testare modelli PPMx dove il clustering è influenzato dalle variabili X.

set.seed(48)

n_main <- 1000

# --- 1. Configurazione Parametri ---

# Probabilità di appartenenza ai cluster (Pesi della mistura)
probs <- c(0.15, 0.25, 0.3, 0.2, 0.1)

# Parametri per la variabile risposta Y (Likelihood)
y_means <- c(-3, -1, 0.5, 2, 3.5)
y_sds   <- c(0.3, 0.6, 0.4, 0.5, 0.3)

# Parametri per le COVARIATE X (Dimensione P=2)
# In un PPMx, ci aspettiamo che il modello riconosca X1 come guida per il partizionamento.

# X1: Covariata INFORMATIVA 
# Le medie sono ben separate tra i cluster; X1 "indica" a quale cluster appartiene Y.
x1_means <- c(-5, -2, 0, 3, 5) 
x1_sds   <- c(1, 1, 1, 1, 1)

# X2: Covariata NON INFORMATIVA (Rumore)
# Le medie sono identiche per tutti i cluster; X2 non aiuta a distinguere i gruppi.
x2_means <- c(0, 0, 0, 0, 0)
x2_sds   <- c(1, 1, 1, 1, 1)


# --- 2. Generazione dei Dati Principali ---

# Inizializzazione vettori
y_main <- numeric(n_main)
x1_main <- numeric(n_main)
x2_main <- numeric(n_main)
true_clusters <- numeric(n_main) # Ground Truth per validazione finale

for(i in 1:n_main){
  # 1. Campionamento del cluster (Categorical distribution via Inverse Transform)
  u <- runif(1)
  cluster <- which(u <= cumsum(probs))[1]
  true_clusters[i] <- cluster
  
  # 2. Generazione congiunta: y, x1, e x2 dipendono dallo stesso indice di cluster
  y_main[i]  <- rnorm(1, mean = y_means[cluster], sd = y_sds[cluster])
  x1_main[i] <- rnorm(1, mean = x1_means[cluster], sd = x1_sds[cluster])
  x2_main[i] <- rnorm(1, mean = x2_means[cluster], sd = x2_sds[cluster])
}


# --- 3. Inserimento Outliers ---
# Creiamo un gruppo di osservazioni isolate sia nello spazio della risposta che delle covariate.
n_out <- 10
y_out  <- rnorm(n_out, mean = 8, sd = 1)   # Outlier alto in Y
x1_out <- rnorm(n_out, mean = 10, sd = 1)  # Outlier estremo in X1
x2_out <- rnorm(n_out, mean = 0, sd = 1)   # Rumore standard in X2

# --- 4. Consolidamento del Dataset ---
y <- c(y_main, y_out)
x <- cbind(x1 = c(x1_main, x1_out), 
           x2 = c(x2_main, x2_out))

# Definiamo il vettore dei cluster reali (gli outlier vengono etichettati come cluster 0)
true_clusters_final <- c(true_clusters, rep(0, n_out))

# --- 5. Salvataggio ---
# Esportiamo y (risposta) e x (matrice covariate)
save(y, x, file = "data_generated.RData")


# --- 6. Analisi Grafica Esplorativa ---
par(mfrow = c(1, 2))

# Plot A: Distribuzione marginale di Y
hist(y, breaks = 50, col = "lightblue", border = "white",
     main = "Distribuzione di Y", xlab = "y", prob = TRUE)
lines(density(y), col = "blue", lwd = 2)

# Plot B: Relazione Y-X1 (Visualizzazione della separabilità)
# Se il modello PPMx funziona, i gruppi di colori simili dovrebbero essere raggruppati insieme.
plot(x[,1], y, col = true_clusters_final + 1, pch = 20, 
     main = "Relazione Informativa: Y vs X1", 
     xlab = "x1 (covariata informativa)", ylab = "y")
legend("topleft", legend = c("Outlier", paste("Cluster", 1:5)), 
       col = 1:6, pch = 20, bty = "n", cex = 0.8)

par(mfrow = c(1, 1))