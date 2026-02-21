# --- generate_data_complex.R ---
# Script per la generazione di dati sintetici complessi (Mistura di Gaussiane + Outliers)
# Utile per testare modelli di clustering non parametrici (es. DPMM)

set.seed(48) # Fissiamo il seed per la riproducibilità dei risultati

# --- 1. Configurazione Parametri ---
n <- 1000  # Numero di osservazioni principali

# Definiamo i parametri per i 5 cluster (componenti della mistura)
means <- c(-3, -1, 0.5, 2, 3.5)       # Medie (centroidi)
sds   <- c(0.3, 0.6, 0.4, 0.5, 0.3)  # Deviazioni standard (eteroschedasticità)
probs <- c(0.15, 0.25, 0.3, 0.2, 0.1) # Pesi della mistura (devono sommare a 1)

# --- 2. Generazione dei Dati ---
y <- numeric(n)

for(i in 1:n) {
  # Estrazione del cluster di appartenenza tramite campionamento multinomiale
  # Usiamo il metodo della cumulativa per mappare un valore Uniforme(0,1) su un cluster
  u <- runif(1)
  cluster <- which(u <= cumsum(probs))[1]
  
  # Generazione del valore y data l'appartenenza al cluster (Likelihood Gaussiana)
  y[i] <- rnorm(1, mean = means[cluster], sd = sds[cluster])
}

# --- 3. Inserimento Outliers ---
# Aggiungiamo una piccola frazione di dati (1%) che non appartiene alla mistura principale
# Questo serve a testare la robustezza del modello di clustering
n_out <- 10
outliers <- rnorm(n_out, mean = 6, sd = 1)
y <- c(y, outliers)

# --- 4. Esportazione e Visualizzazione ---
# Salvataggio del dataset in formato binario RData per analisi successive
save(y, file = "data_generated.RData")

# Visualizzazione rapida della densità empirica dei dati generati
hist(y, 
     breaks = 50, 
     col = "lightblue", 
     border = "white",
     main = "Synthetic univariate data (complex)", 
     xlab = "Valore osservato (y)",
     prob = TRUE) # prob=TRUE per vedere la densità invece della frequenza

# Aggiunta di una curva di densità kernel per evidenziare i picchi
lines(density(y), col = "darkblue", lwd = 2)

