# --- load_data.R (Analisi dei dati reali dei Laghi) ---
# Questo script gestisce l'importazione, la pulizia e la trasformazione 
# dei dati idrologici per renderli idonei al modello PPMx.

# 1. Caricamento Dataset
# Il file "reg_data.csv" contiene misurazioni reali (es. indici di portata o livelli)
dati_laghi_raw <- read.csv("reg_data.csv")
cat("Dati grezzi caricati con successo.\n")

# Escludiamo STAID (ID stazione) dal calcolo statistico ma lo conserviamo come metadato
nomi_colonne <- setdiff(names(dati_laghi_raw), "STAID")

# -----------------------------------------------------------------
# 2. ANALISI ESPLORATIVA (Dati Originali)
# -----------------------------------------------------------------
# Molti dati ambientali presentano forti asimmetrie (skewness) e code pesanti.
# Questi grafici servono a giustificare la necessità di una trasformazione.

cat("Analisi visiva delle distribuzioni originali...\n")

# Istogrammi per verificare l'asimmetria
par(mfrow = c(1, 8), oma = c(0, 0, 2, 0)) 
for (col_name in nomi_colonne) {
  hist(dati_laghi_raw[[col_name]], main = col_name, 
       xlab = col_name, col = "salmon", breaks = 30, border = "white")
}
mtext("Istogrammi - Dati GREZZI", outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1))

# Q-Q Plots per verificare la normalità
par(mfrow = c(3, 3), oma = c(0, 0, 2, 0))
for (col_name in nomi_colonne) {
  qqnorm(dati_laghi_raw[[col_name]], main = col_name)
  qqline(dati_laghi_raw[[col_name]], col = "red")
}
mtext("Q-Q Plots - Dati GREZZI", outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1))

# -----------------------------------------------------------------
# 3. TRASFORMAZIONE LOGARITMICA ROBUSTA
# -----------------------------------------------------------------
# Poiché il kernel del modello è Gaussiano, la normalizzazione dei dati 
# migliora drasticamente la convergenza dell'MCMC.

cat("Applicazione trasformazione Log-Shift (log1p)...\n")

# Calcolo del minimo globale per gestire eventuali valori negativi
# Questo assicura che l'argomento del logaritmo sia sempre strettamente positivo (> 0)
matrice_dati <- as.matrix(dati_laghi_raw[nomi_colonne])
min_globale <- min(matrice_dati, na.rm = TRUE)

cat("Valore minimo globale individuato:", min_globale, "\n")

# Nuovo dataframe per i dati trasformati
dati_laghi_log <- dati_laghi_raw

# Applichiamo: log( (x - min) + 1 )
# log1p è numericamente più stabile per valori vicini allo zero
for (col_name in nomi_colonne) {
  dati_laghi_log[[col_name]] <- log1p(dati_laghi_raw[[col_name]] - min_globale)
}

cat("Trasformazione completata.\n")

# -----------------------------------------------------------------
# 4. STANDARDIZZAZIONE (Scaling)
# -----------------------------------------------------------------
# In un PPMx, la similarità delle covariate X è basata su distanze.
# Standardizzare (media=0, SD=1) è fondamentale affinché ogni variabile
# pesi equamente nel processo di clustering, indipendentemente dall'unità di misura.

cat("Standardizzazione (Z-score scaling)...\n")

dati_laghi_scaled <- dati_laghi_log 
dati_laghi_scaled[nomi_colonne] <- scale(dati_laghi_scaled[nomi_colonne])

cat("Dati log-trasformati e standardizzati pronti.\n")

# -----------------------------------------------------------------
# 5. PREPARAZIONE PER IL MODELLO
# -----------------------------------------------------------------

# Definizione Variabile Risposta (Y)
# Usiamo 'max90' come variabile target per il clustering
y <- dati_laghi_scaled$max90  

# Definizione Covariate (X)
# Escludiamo la risposta e l'identificativo dalle covariate
x_df <- subset(dati_laghi_scaled, select = -c(max90, STAID))
x <- as.matrix(x_df)

# Riepilogo finale
n <- length(y)
n_covariate <- ncol(x) 
cat("\n--- Setup Finale ---\n")
cat("Osservazioni (n):", n, "\n")
cat("Covariate (P):", n_covariate, "\n")
cat("Target Y: max90\n")
cat("--------------------\n")

# Salvataggio dati processati per lo script run_model.R
save(y, x, nomi_colonne, min_globale, file = "data_laghi_processed.RData")