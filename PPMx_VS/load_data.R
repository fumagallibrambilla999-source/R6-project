# --- load_data.R ---

# 1. CARICAMENTO DATI
# Carica il dataset CSV. Assicurati che il file sia nella cartella di lavoro.
dati_laghi_raw <- read.csv("reg_data.csv")
cat("Dati grezzi caricati.\n")

# 2. DEFINIZIONE DELLE COLONNE
# Identifichiamo le variabili che NON devono subire trasformazioni logaritmiche (ID o variabili già trattate)
colonne_escluse <- c("STAID", "T_AVG_SITE", "RH_BASIN", "RRMEDIAN", "T_AVG_BASIN")

# Selezioniamo per differenza le colonne da analizzare e trasformare
nomi_colonne <- setdiff(names(dati_laghi_raw), colonne_escluse)

# -----------------------------------------------
# --- 3. TRASFORMAZIONE LOGARITMICA (ROBUSTA) ---
# -----------------------------------------------
cat("Applicazione trasformazione logaritmica robusta (shift + log1p)...\n")

# Estraiamo i dati da trasformare in una matrice per calcolare il minimo globale
# Nota: na.rm = TRUE è fondamentale se il dataset contiene celle vuote
matrice_dati <- as.matrix(dati_laghi_raw[nomi_colonne])
min_globale <- min(matrice_dati, na.rm = TRUE)

cat("Valore minimo globale trovato nelle colonne da trasformare:", min_globale, "\n")

# Calcolo dello SHIFT: 
# Se il valore minimo è negativo, trasliamo i dati per renderli positivi prima del logaritmo.
shift <- 0
if (min_globale < 0) {
  shift <- -min_globale 
  cat("I dati verranno spostati di:", shift, "per evitare log di numeri negativi.\n")
}

# Creiamo una copia del dataset per non sovrascrivere i dati originali
dati_laghi_log <- dati_laghi_raw

# Applicazione della trasformazione: log1p(x) calcola log(1 + x) garantendo stabilità per valori vicini a zero
for (col_name in nomi_colonne) {
  dati_laghi_log[[col_name]] <- log1p(dati_laghi_raw[[col_name]] + shift)
}

cat("Trasformazione completata.\n")


# -----------------------------------------------
# --- 4. PREPARAZIONE FINALE (Scaling & Split) ---
# -----------------------------------------------
cat("Standardizzazione dati...\n")

# Definiamo quali colonne scalare (solitamente tutte le numeriche tranne l'ID)
# Lo scaling (Z-score) porta media a 0 e deviazione standard a 1, essenziale per modelli tipo LASSO/Ridge.
colonne_da_scalare <- setdiff(names(dati_laghi_raw), "STAID")

dati_laghi_scaled <- dati_laghi_log 
dati_laghi_scaled[colonne_da_scalare] <- scale(dati_laghi_scaled[colonne_da_scalare])

cat("Dati pronti (log-trasformati e standardizzati).\n")

# --- DEFINIZIONE VARIABILI PER IL MODELLO ---

# 1. Variabile Target (Y): La variabile che vogliamo predire (es. max90)
y <- dati_laghi_scaled$max90  

# 2. Covariate (X): Tutte le altre variabili eccetto il Target e l'identificatore (STAID)
x_df <- subset(dati_laghi_scaled, select = -c(max90, STAID))

# Trasformazione in Matrice: Molte librerie di regressione (es. glmnet) richiedono una matrice come input
x <- as.matrix(x_df)

# -----------------------------------------------
# --- 5. REPORT FINALE DI CONTROLLO ---
# -----------------------------------------------
n <- length(y)
n_covariate <- ncol(x) 

cat("----------------------------------\n")
cat("Setup finale pronto per il modello:\n")
cat("N. Osservazioni (righe):", n, "\n")
cat("Variabile Risposta 'y': max90\n")
cat("Numero Covariate 'x':", n_covariate, "\n")
cat("Variabili ESCLUSE dal log:", paste(setdiff(colonne_escluse, "STAID"), collapse=", "), "\n")
cat("Verifica NA: ", anyNA(x), "(se TRUE, controlla i dati prima di procedere)\n")
cat("----------------------------------\n")