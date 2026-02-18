# ==============================================================================
# SCRIPT DI VALIDAZIONE: CLUSTER VS NOISE
# ==============================================================================

#7 COOVARIATE
if(TRUE){
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# 1. Caricamento dati
load("Dataset_Simulati_AppendixE_WithNoise.RData")

# Selezioniamo il primo dataset
d1 <- DATA[[1]]
df_plot <- data.frame(
  Y = d1$Y,
  X1 = d1$X[, 1], # Informativa
  X5 = d1$X[, 5], # Noise
  Cluster = factor(d1$true_clusters)
)

cat("\n==============================================================================\n")
cat(" ANALISI STATISTICA (Dataset 1)\n")
cat("==============================================================================\n")

# ANOVA
cat("\n[ANOVA] Differenza di Y tra i Cluster:\n")
print(summary(aov(Y ~ Cluster, data = df_plot)))

# REGRESSIONE
cat("\n[REGRESSIONE] Impatto delle variabili (X1-X4 vere, X5-X7 noise):\n")
lm_res <- lm(d1$Y ~ d1$X)
names(lm_res$coefficients) <- c("Intercept", paste0("X", 1:7))
print(summary(lm_res))

cat("\n==============================================================================\n")
cat(" GENERAZIONE GRAFICI DI CONFRONTO\n")
cat("==============================================================================\n")

# GRAFICO A: Boxplot Y per Cluster (Richiesto)
p1 <- ggplot(df_plot, aes(x = Cluster, y = Y, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = Cluster), width = 0.2, size = 1.2, alpha = 0.6) +
  theme_minimal() +
  labs(title = "A. Distribuzione Y per Cluster",
       subtitle = "I gruppi sono chiaramente separati verticalmente",
       x = "Cluster", y = "Risposta Y") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  theme(legend.position = "none")

# GRAFICO B: Y vs X1 (Variabile INFORMATIVA)
p2 <- ggplot(df_plot, aes(x = X1, y = Y, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8) +
  theme_minimal() +
  labs(title = "B. Segnale: Y vs X1 (Informativa)",
       subtitle = "Esiste una relazione chiara e distinta per gruppo",
       x = "X1 (Informativa)", y = "Risposta Y") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

# GRAFICO C: Y vs X5 (Variabile NOISE)
p3 <- ggplot(df_plot, aes(x = X5, y = Y, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8, linetype = "dotted") +
  theme_minimal() +
  labs(title = "C. Rumore: Y vs X5 (Noise)",
       subtitle = "Nessun pattern: i colori sono mescolati e le rette piatte",
       x = "X5 (Noise / Rumore)", y = "Risposta Y") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

# Organizzazione finale in una griglia
grid.arrange(p1, p2, p3, ncol = 1)

cat("\nVisualizzazione completata.\n")
}

# ==============================================================================
# SCRIPT DI VALIDAZIONE: CLUSTER VS NOISE (Versione 5 Covariate)
# ==============================================================================

if(FALSE){
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# 1. Caricamento dati
load("Dataset_Simulati_AppendixE_WithNoise.RData")

# Selezioniamo il primo dataset
d1 <- DATA[[1]]
df_plot <- data.frame(
  Y = d1$Y,
  X1 = d1$X[, 1], # Informativa
  X5 = d1$X[, 5], # Noise
  Cluster = factor(d1$true_clusters)
)

cat("\n==============================================================================\n")
cat(" ANALISI STATISTICA (Dataset 1)\n")
cat("==============================================================================\n")

# ANOVA
cat("\n[ANOVA] Differenza di Y tra i Cluster:\n")
print(summary(aov(Y ~ Cluster, data = df_plot)))

# REGRESSIONE
# cat("\n[REGRESSIONE] Impatto delle variabili (X1-X4 vere, X5-X7 noise):\n")
cat("\n[REGRESSIONE] Impatto delle variabili (X1-X4 vere, X5 noise):\n")
lm_res <- lm(d1$Y ~ d1$X)

# names(lm_res$coefficients) <- c("Intercept", paste0("X", 1:7))
names(lm_res$coefficients) <- c("Intercept", paste0("X", 1:5)) # Aggustato per 5 covariate

print(summary(lm_res))

cat("\n==============================================================================\n")
cat(" GENERAZIONE GRAFICI DI CONFRONTO\n")
cat("==============================================================================\n")

# GRAFICO A: Boxplot Y per Cluster
p1 <- ggplot(df_plot, aes(x = Cluster, y = Y, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = Cluster), width = 0.2, size = 1.2, alpha = 0.6) +
  theme_minimal() +
  labs(title = "A. Distribuzione Y per Cluster",
       subtitle = "I gruppi sono chiaramente separati verticalmente",
       x = "Cluster", y = "Risposta Y") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  theme(legend.position = "none")

# GRAFICO B: Y vs X1 (Variabile INFORMATIVA)
p2 <- ggplot(df_plot, aes(x = X1, y = Y, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8) +
  theme_minimal() +
  labs(title = "B. Segnale: Y vs X1 (Informativa)",
       subtitle = "Esiste una relazione chiara e distinta per gruppo",
       x = "X1 (Informativa)", y = "Risposta Y") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

# GRAFICO C: Y vs X5 (Variabile NOISE)
p3 <- ggplot(df_plot, aes(x = X5, y = Y, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, size = 0.8, linetype = "dotted") +
  theme_minimal() +
  labs(title = "C. Rumore: Y vs X5 (Noise)",
       subtitle = "Nessun pattern: i colori sono mescolati e le rette piatte",
       x = "X5 (Noise / Rumore)", y = "Risposta Y") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

# Organizzazione finale in una griglia
grid.arrange(p1, p2, p3, ncol = 1)

cat("\nVisualizzazione completata.\n")
}