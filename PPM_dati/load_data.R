# --- load_data.R (Modificato con Trasformazione Logaritmica ROBUSTA) ---

# 1. Carica i dati dal tuo CSV
dati_laghi_raw <- read.csv("reg_data.csv")
cat("Dati grezzi caricati.\n")

# Colonne da analizzare e trasformare (tutte tranne STAID)
nomi_colonne <- setdiff(names(dati_laghi_raw), "STAID")

# -----------------------------------------------
# --- 2. ANALISI (Dati GREZZI) ---
# -----------------------------------------------
cat("Visualizzazione dati GREZZI (chiudere la finestra grafica per continuare)...\n")
# ... (il tuo codice per i plot dei dati grezzi va qui, l'ho omesso per brevità) ...
par(mfrow = c(1, 8), oma = c(0, 0, 2, 0)) 
for (col_name in nomi_colonne) {
  hist(dati_laghi_raw[[col_name]], main = col_name, xlab = col_name, col = "salmon", breaks = 30)
}
mtext("Istogrammi - Dati GREZZI", outer = TRUE, cex = 1.5)
par(mfrow = c(1, 1))

par(mfrow = c(3, 3), oma = c(0, 0, 2, 0))
for (col_name in nomi_colonne) {
  qqnorm(dati_laghi_raw[[col_name]], main = col_name)
  qqline(dati_laghi_raw[[col_name]], col = "red")
}
mtext("Q-Q Plots - Dati GREZZI", outer = TRUE, cex = 1.5)
par(mfrow = c(1, 1))


# -----------------------------------------------
# --- 3. TRASFORMAZIONE LOGARITMICA (ROBUSTA) ---
# -----------------------------------------------
cat("Applicazione trasformazione logaritmica robusta (shift + log1p)...\n")

# --- INIZIO MODIFICA ROBUSTA ---

# 1. Trova il valore minimo globale tra tutte le colonne di dati
#    convertiamo prima in matrice per assicurarci che min() funzioni correttamente
matrice_dati <- as.matrix(dati_laghi_raw[nomi_colonne])
min_globale <- min(matrice_dati, na.rm = TRUE)

cat("Valore minimo globale trovato:", min_globale, "\n")

# 2. Definisci lo "spostamento". Se il min è >= 0, non spostiamo nulla.
#    Se il min è negativo (es. -5), lo spostamento sarà 5.
shift <- 0
if (min_globale < 0) {
  # Usiamo 'abs(min_globale)' o semplicemente '-min_globale'
  shift <- -min_globale 
  cat("I dati verranno spostati di:", shift, "prima del log.\n")
}

# Creiamo un nuovo dataframe per i dati trasformati
dati_laghi_log <- dati_laghi_raw

# 3. Applichiamo log1p(x + shift)
#    Se shift = 0, questo è log1p(x) (come prima)
#    Se shift = 5, questo è log1p(x + 5)
for (col_name in nomi_colonne) {
  # Nota: NON usiamo log1p(x - min_globale) perché se min_globale = 0,
  # avremmo log1p(x), che è corretto.
  # Se min_globale = -5, avremmo log1p(x - (-5)) = log1p(x + 5). Corretto.
  # Quindi, la formula più semplice è log1p(x - min_globale)
  
  dati_laghi_log[[col_name]] <- log1p(dati_laghi_raw[[col_name]] - min_globale)
}

# --- FINE MODIFICA ROBUSTA ---

cat("Trasformazione completata.\n")


# -----------------------------------------------
# --- 4. ANALISI (Dati TRASFORMATI) ---
# -----------------------------------------------
cat("Visualizzazione dati TRASFORMATI (chiudere la finestra grafica per continuare)...\n")

# Istogrammi (Trasformati)
par(mfrow = c(1, 8), oma = c(0, 0, 2, 0))
for (col_name in nomi_colonne) {
  hist(dati_laghi_log[[col_name]], 
       main = paste("log(..", col_name, "..+k)"), 
       xlab = col_name, 
       col = "lightblue", 
       breaks = 30)
}
mtext("Istogrammi - Dati TRASFORMATI (Shift + log1p)", outer = TRUE, cex = 1.5)
par(mfrow = c(1, 1))

# Q-Q Plots (Trasformati)
par(mfrow = c(3, 3), oma = c(0, 0, 2, 0))
for (col_name in nomi_colonne) {
  qqnorm(dati_laghi_log[[col_name]], main = paste("log(..", col_name, "..+k)"))
  qqline(dati_laghi_log[[col_name]], col = "blue")
}
mtext("Q-Q Plots - Dati TRASFORMATI (Shift + log1p)", outer = TRUE, cex = 1.5)
par(mfrow = c(1, 1))


# -----------------------------------------------
# --- 5. PREPARAZIONE FINALE (Scaling & Split) ---
# -----------------------------------------------
cat("Standardizzazione dati TRASFORMATI...\n")

colonne_da_scalare <- nomi_colonne
dati_laghi_scaled <- dati_laghi_log 
dati_laghi_scaled[colonne_da_scalare] <- scale(dati_laghi_scaled[colonne_da_scalare])

cat("Dati log-trasformati e standardizzati pronti.\n")

# 1. Seleziona la tua variabile RISPOSTA (Y)
y <- dati_laghi_scaled$max90  

# 2. Seleziona le tue COVARIATE (X)
x_df <- subset(dati_laghi_scaled, select = -c(max90, STAID))
x <- as.matrix(x_df)

# Controllo finale
n <- length(y)
n_covariate <- ncol(x) 
cat("----------------------------------\n")
cat("Setup finale pronto per il modello:\n")
cat(n, "osservazioni.\n")
cat("Risposta 'y': max90 (spostata, log-trasformata e standardizzata)\n")
cat("Covariate 'x':", n_covariate, "colonne (spostate, log-trasformate e standardizzate).\n")
cat("----------------------------------\n")