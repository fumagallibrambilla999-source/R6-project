# --- load_data.R ---

# 1. Carica i dati dal tuo CSV
dati_laghi_raw <- read.csv("reg_data.csv")
cat("Dati grezzi caricati.\n")

# --- MODIFICA QUI: Definiamo le colonne da ESCLUDERE ---
colonne_escluse <- c("STAID", "T_AVG_SITE", "RH_BASIN", "RRMEDIAN","T_AVG_BASIN")

# Colonne da analizzare e trasformare (tutte tranne quelle escluse)
nomi_colonne <- setdiff(names(dati_laghi_raw), colonne_escluse)

# -----------------------------------------------
# --- 3. TRASFORMAZIONE LOGARITMICA (ROBUSTA) ---
# -----------------------------------------------
cat("Applicazione trasformazione logaritmica robusta (shift + log1p)...\n")

# Calcoliamo il minimo solo sulle colonne che andremo effettivamente a trasformare
matrice_dati <- as.matrix(dati_laghi_raw[nomi_colonne])
min_globale <- min(matrice_dati, na.rm = TRUE)

cat("Valore minimo globale trovato nelle colonne da trasformare:", min_globale, "\n")

shift <- 0
if (min_globale < 0) {
  shift <- -min_globale 
  cat("I dati verranno spostati di:", shift, "prima del log.\n")
}

dati_laghi_log <- dati_laghi_raw

# Applichiamo log1p(x + shift) SOLO alle colonne in nomi_colonne
for (col_name in nomi_colonne) {
  dati_laghi_log[[col_name]] <- log1p(dati_laghi_raw[[col_name]] + shift)
}

cat("Trasformazione completata.\n")


# -----------------------------------------------
# --- 5. PREPARAZIONE FINALE (Scaling & Split) ---
# -----------------------------------------------
cat("Standardizzazione dati...\n")

# Decidi se vuoi scalare TUTTE le numeriche o solo quelle log-trasformate.
# Di solito si scalano tutte le covariate (tranne l'ID) per il modello.
colonne_da_scalare <- setdiff(names(dati_laghi_raw), "STAID")

dati_laghi_scaled <- dati_laghi_log 
dati_laghi_scaled[colonne_da_scalare] <- scale(dati_laghi_scaled[colonne_da_scalare])

cat("Dati pronti (alcuni log-trasformati, tutti standardizzati).\n")

# 1. Seleziona la tua variabile RISPOSTA (Y)
y <- dati_laghi_scaled$max90  

# 2. Seleziona le tue COVARIATE (X)
# Escludiamo la Y e l'ID
x_df <- subset(dati_laghi_scaled, select = -c(max90, STAID))
x <- as.matrix(x_df)

# Controllo finale
n <- length(y)
n_covariate <- ncol(x) 
cat("----------------------------------\n")
cat("Setup finale pronto per il modello:\n")
cat(n, "osservazioni.\n")
cat("Risposta 'y': max90\n")
cat("Covariate 'x':", n_covariate, "colonne.\n")
cat("Colonne NON trasformate logaritmicamente:", paste(setdiff(colonne_escluse, "STAID"), collapse=", "), "\n")
cat("----------------------------------\n")