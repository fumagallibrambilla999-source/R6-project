# ==============================================================================
# SCRIPT DI VALIDAZIONE: CLUSTER VS NOISE
# ==============================================================================
# Questo script verifica la qualità dei dati simulati attraverso tre pilastri:
# 1. Analisi della varianza (ANOVA) per la separazione dei cluster.
# 2. Regressione lineare classica (LM) come benchmark di selezione variabili.
# 3. Visualizzazione comparativa Segnale vs Rumore.
# ==============================================================================

if(TRUE){
  # --- Caricamento librerie grafiche ---
  if(!require(ggplot2)) install.packages("ggplot2")
  if(!require(gridExtra)) install.packages("gridExtra")
  library(ggplot2)
  library(gridExtra)
  
  # 1. Caricamento dati generati
  load("Dataset_Simulati_AppendixE_WithNoise.RData")
  
  # Utilizziamo il primo dataset della lista come campione rappresentativo
  d1 <- DATA[[1]]
  df_plot <- data.frame(
    Y = d1$Y,
    X1 = d1$X[, 1], # Variabile con segnale (Informativa)
    X5 = d1$X[, 5], # Variabile senza segnale (Noise)
    Cluster = factor(d1$true_clusters)
  )
  
  cat("\n==============================================================================\n")
  cat(" ANALISI STATISTICA DI VALIDAZIONE (Dataset 1)\n")
  cat("==============================================================================\n")
  
  # --- ANOVA: Validazione dei Cluster ---
  # Verifichiamo se la variabile risposta Y differisce significativamente tra i cluster veri.
  # Un p-value basso qui conferma che il clustering basato su Y ha senso.
  cat("\n[1. ANOVA] Differenza di Y tra i Cluster:\n")
  print(summary(aov(Y ~ Cluster, data = df_plot)))
  
  # --- REGRESSIONE: Benchmark di Variable Selection ---
  # Eseguiamo una regressione OLS. Se X5-X7 risultano non significative (p-value alto),
  # il dataset è costruito correttamente per "sfidare" il modello PPMx-VS.
  cat("\n[2. REGRESSIONE OLS] Impatto delle variabili (X1-X4 attese vere, X5-X7 attese noise):\n")
  lm_res <- lm(d1$Y ~ d1$X)
  names(lm_res$coefficients) <- c("Intercept", paste0("X", 1:7))
  print(summary(lm_res))
  
  cat("\n==============================================================================\n")
  cat(" GENERAZIONE GRAFICI DI CONFRONTO\n")
  cat("==============================================================================\n")
  
  # --- GRAFICO A: Boxplot Y per Cluster ---
  # Conferma visivamente la separazione delle medie dei cluster.
  p1 <- ggplot(df_plot, aes(x = Cluster, y = Y, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(aes(color = Cluster), width = 0.2, size = 1.2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "A. Distribuzione Y per Cluster",
         subtitle = "Separazione verticale: indica la capacità di Y di discriminare i gruppi",
         x = "Cluster ID", y = "Risposta Y") +
    scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
    theme(legend.position = "none")
  
  # --- GRAFICO B: Segnale (Y vs X1) ---
  # Mostra la relazione lineare all'interno dei cluster (pendenze diverse).
  p2 <- ggplot(df_plot, aes(x = X1, y = Y, color = Cluster)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, size = 1) +
    theme_minimal() +
    labs(title = "B. Segnale: Y vs X1 (Informativa)",
         subtitle = "Interazioni: la relazione X->Y cambia a seconda del cluster",
         x = "X1 (Informativa)", y = "Risposta Y") +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))
  
  # --- GRAFICO C: Rumore (Y vs X5) ---
  # Conferma l'assenza di pattern. Se le linee sono sovrapposte o piatte, la noise è corretta.
  p3 <- ggplot(df_plot, aes(x = X5, y = Y, color = Cluster)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, size = 1, linetype = "dotted") +
    theme_minimal() +
    labs(title = "C. Rumore: Y vs X5 (Noise)",
         subtitle = "Assenza di pattern: X5 non contribuisce alla spiegazione di Y",
         x = "X5 (Noise)", y = "Risposta Y") +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))
  
  # Visualizzazione combinata
  grid.arrange(p1, p2, p3, ncol = 1)
  
  cat("\nValidazione completata. I dati sono pronti per l'MCMC.\n")
}