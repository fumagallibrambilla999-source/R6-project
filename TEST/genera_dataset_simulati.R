# ==============================================================================
# GENERAZIONE DATASET SIMULATI - Appendice E + Noise Covariates
# ==============================================================================



# Struttura:
# - 3 gruppi con dimensioni (75, 75, 50) = 200 osservazioni totali
# - 4 covariate informative: 2 continue (x1, x2) + 2 binarie (x3, x4)
# - 3 covariate noise: continue indipendenti da Y
# - Totale: 7 covariate (4 vere + 3 noise)
# ==============================================================================


if(TRUE){
set.seed(12345)  # Per riproducibilità

# ==============================================================================
# PARAMETRI GLOBALI
# ==============================================================================
N_DATASETS <- 10      # Numero di dataset da generare
N_TOTAL <- 200        # Osservazioni totali per dataset
N_GROUPS <- 3         # Numero di gruppi veri
GROUP_SIZES <- c(75, 75, 50)  # Dimensioni dei gruppi

# ==============================================================================
# PARAMETRI DEI GRUPPI (da Appendice E)
# ==============================================================================

# Gruppo 1
mu_1 <- c(-3, 3)      # Medie per (x1, x2)
q_1 <- 0.1            # Prob. per binarie (x3, x4)
beta_1 <- c(1, 5, 2, 1, 0)  # Coefficienti: intercetta + 4 covariate

# Gruppo 2
mu_2 <- c(0, 0)
q_2 <- 0.5
beta_2 <- c(4, 2, -2, 1, -1)

# Gruppo 3
mu_3 <- c(3, 3)
q_3 <- 0.9
beta_3 <- c(-1, -5, -2, -1, 1)

# Matrice di covarianza per covariate continue (identità scalata)
Sigma_cov <- 0.5 * diag(2)

# Varianza residua per Y
sigma_y <- 0.5

# Varianze per le covariate noise
sigma_noise <- c(1, 2, 3)  # Varianze crescenti

# ==============================================================================
# FUNZIONE PER GENERARE UN SINGOLO DATASET
# ==============================================================================

genera_dataset_simulato <- function(dataset_id, verbose = TRUE) {
  
  if(verbose) {
    cat("\n========================================")
    cat("\n Generazione Dataset", dataset_id)
    cat("\n========================================\n")
  }
  
  # ----------------------------------------------------------------------------
  # 1. Inizializzazione
  # ----------------------------------------------------------------------------
  n_total <- sum(GROUP_SIZES)
  
  # Matrici di output
  X_all <- matrix(NA, nrow = n_total, ncol = 7)  # 4 vere + 3 noise
  Y_all <- numeric(n_total)
  true_labels <- integer(n_total)
  
  colnames(X_all) <- paste0("X", 1:7)
  
  # ----------------------------------------------------------------------------
  # 2. Generazione per ogni gruppo
  # ----------------------------------------------------------------------------
  row_idx <- 1
  
  for(group in 1:N_GROUPS) {
    
    n_group <- GROUP_SIZES[group]
    
    # Selezione parametri del gruppo
    if(group == 1) {
      mu_g <- mu_1
      q_g <- q_1
      beta_g <- beta_1
    } else if(group == 2) {
      mu_g <- mu_2
      q_g <- q_2
      beta_g <- beta_2
    } else {
      mu_g <- mu_3
      q_g <- q_3
      beta_g <- beta_3
    }
    
    # ------------------------------------------------------------------------
    # 2a. Covariate CONTINUE INFORMATIVE (x1, x2)
    # ------------------------------------------------------------------------
    X_cont <- MASS::mvrnorm(n = n_group, mu = mu_g, Sigma = Sigma_cov)
    
    # ------------------------------------------------------------------------
    # 2b. Covariate BINARIE INFORMATIVE (x3, x4)
    # ------------------------------------------------------------------------
    X_bin <- matrix(
      rbinom(n = 2 * n_group, size = 1, prob = q_g),
      nrow = n_group,
      ncol = 2
    )
    
    # ------------------------------------------------------------------------
    # 2c. Covariate NOISE (x5, x6, x7) - INDIPENDENTI DA Y
    # ------------------------------------------------------------------------
    X_noise <- cbind(
      rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[1])),  # x5 ~ N(0, 1)
      rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[2])),  # x6 ~ N(0, 2)
      rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[3]))   # x7 ~ N(0, 3)
    )
    
    # ------------------------------------------------------------------------
    # 2d. Combina tutte le covariate
    # ------------------------------------------------------------------------
    X_group <- cbind(X_cont, X_bin, X_noise)
    
    # ------------------------------------------------------------------------
    # 2e. Genera risposta Y (dipende SOLO da x1, x2, x3, x4)
    # ------------------------------------------------------------------------
    # Design matrix: intercetta + 4 covariate informative
    X_design <- cbind(1, X_cont, X_bin)  # 5 colonne
    
    # Media condizionale
    mu_y <- X_design %*% beta_g
    
    # Risposta con errore gaussiano
    Y_group <- rnorm(n_group, mean = mu_y, sd = sqrt(sigma_y))
    
    # ------------------------------------------------------------------------
    # 2f. Salva nel dataset completo
    # ------------------------------------------------------------------------
    idx_range <- row_idx:(row_idx + n_group - 1)
    X_all[idx_range, ] <- X_group
    Y_all[idx_range] <- Y_group
    true_labels[idx_range] <- group
    
    row_idx <- row_idx + n_group
    
    if(verbose) {
      cat("  Gruppo", group, "| N =", n_group, "| rows:", 
          min(idx_range), "-", max(idx_range), "\n")
    }
  }
  
  # ----------------------------------------------------------------------------
  # 3. Permutazione casuale (importante per MCMC)
  # ----------------------------------------------------------------------------
  perm_idx <- sample(1:n_total)
  X_all <- X_all[perm_idx, ]
  Y_all <- Y_all[perm_idx]
  true_labels <- true_labels[perm_idx]
  
  # ----------------------------------------------------------------------------
  # 4. Output
  # ----------------------------------------------------------------------------
  dataset <- list(
    X = X_all,
    Y = Y_all,
    true_clusters = true_labels,
    selected.cols = 1:4,  # Le prime 4 covariate sono quelle vere
    noise.cols = 5:7,     # Le ultime 3 sono noise
    group_params = list(
      mu = list(mu_1, mu_2, mu_3),
      q = c(q_1, q_2, q_3),
      beta = list(beta_1, beta_2, beta_3)
    ),
    info = list(
      n_total = n_total,
      n_groups = N_GROUPS,
      group_sizes = GROUP_SIZES,
      dataset_id = dataset_id
    )
  )
  
  if(verbose) {
    cat("\n  Dataset generato:")
    cat("\n  - Dimensioni X:", nrow(X_all), "x", ncol(X_all))
    cat("\n  - Covariate vere: 1, 2, 3, 4")
    cat("\n  - Covariate noise: 5, 6, 7")
    cat("\n  - Gruppi veri:", unique(true_labels))
    cat("\n")
  }
  
  return(dataset)
}

# ==============================================================================
# GENERAZIONE DI TUTTI I DATASET
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat(" GENERAZIONE", N_DATASETS, "DATASET SIMULATI\n")
cat("==============================================================================\n")
cat(" Struttura per dataset:\n")
cat("  - N osservazioni:", N_TOTAL, "\n")
cat("  - N gruppi veri:", N_GROUPS, "\n")
cat("  - Dimensioni gruppi:", paste(GROUP_SIZES, collapse = ", "), "\n")
cat("  - Covariate informative: 4 (2 continue + 2 binarie)\n")
cat("  - Covariate noise: 3 (indipendenti da Y)\n")
cat("  - Totale covariate: 7\n")
cat("==============================================================================\n")

DATA <- vector("list", N_DATASETS)

for(d in 1:N_DATASETS) {
  DATA[[d]] <- genera_dataset_simulato(dataset_id = d, verbose = TRUE)
}

# ==============================================================================
# DIAGNOSTICA RAPIDA
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat(" DIAGNOSTICA DATASET\n")
cat("==============================================================================\n")

# Esempio: primo dataset
d1 <- DATA[[1]]

cat("\n Dataset 1 - Statistiche descrittive:\n")
cat(" -------------------------------------------\n")
cat(" Y: mean =", round(mean(d1$Y), 2), 
    ", sd =", round(sd(d1$Y), 2), "\n")

for(j in 1:7) {
  is_informative <- j %in% d1$selected.cols
  cat(sprintf(" X%d (%s): mean = %6.2f, sd = %6.2f\n", 
              j, 
              ifelse(is_informative, "INFO", "NOISE"),
              mean(d1$X[, j]),
              sd(d1$X[, j])))
}

cat("\n Correlazioni Y con covariate:\n")
cat(" -------------------------------------------\n")
for(j in 1:7) {
  is_informative <- j %in% d1$selected.cols
  cor_yj <- cor(d1$Y, d1$X[, j])
  cat(sprintf(" cor(Y, X%d) = %6.3f %s\n", 
              j, 
              cor_yj,
              ifelse(is_informative, "✓", " ")))
}

# ==============================================================================
# SALVATAGGIO
# ==============================================================================

output_filename <- "Dataset_Simulati_AppendixE_WithNoise.RData"
save(DATA, file = output_filename)

cat("\n")
cat("==============================================================================\n")
cat(" SALVATAGGIO COMPLETATO\n")
cat("==============================================================================\n")
cat(" File:", output_filename, "\n")
cat(" Numero di dataset:", length(DATA), "\n")
cat("\n Per caricare:\n")
cat("   load('", output_filename, "')\n", sep = "")
cat("\n Struttura oggetto DATA:\n")
cat("   DATA[[i]]$X            - Matrice covariate (200 x 7)\n")
cat("   DATA[[i]]$Y            - Vettore risposta (200)\n")
cat("   DATA[[i]]$true_clusters - Labels veri dei gruppi\n")
cat("   DATA[[i]]$selected.cols - Indici covariate informative (1:4)\n")
cat("   DATA[[i]]$noise.cols    - Indici covariate noise (5:7)\n")
cat("==============================================================================\n\n")

# ==============================================================================
# VISUALIZZAZIONE ESEMPIO (opzionale)
# ==============================================================================

if(require(ggplot2, quietly = TRUE)) {
  
  cat("Creazione grafici diagnostici...\n")
  
  # Dataset 1
  d1 <- DATA[[1]]
  df_plot <- data.frame(
    Y = d1$Y,
    X1 = d1$X[, 1],
    X2 = d1$X[, 2],
    X5_noise = d1$X[, 5],
    X6_noise = d1$X[, 6],
    Cluster = factor(d1$true_clusters)
  )
  
  # Plot 1: Y vs X1 (informativa)
  p1 <- ggplot(df_plot, aes(x = X1, y = Y, color = Cluster)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_minimal() +
    labs(title = "Dataset 1: Y vs X1 (informativa)",
         subtitle = "I 3 gruppi veri sono visibili") +
    theme(legend.position = "bottom")
  
  # Plot 2: Y vs X5 (noise)
  p2 <- ggplot(df_plot, aes(x = X5_noise, y = Y, color = Cluster)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_minimal() +
    labs(title = "Dataset 1: Y vs X5 (noise)",
         subtitle = "Nessuna relazione con Y") +
    theme(legend.position = "bottom")
  
  # Salva plots
  ggsave("Dataset1_Y_vs_X1_informative.png", plot = p1, 
         width = 8, height = 6)
  ggsave("Dataset1_Y_vs_X5_noise.png", plot = p2, 
         width = 8, height = 6)
  
  cat("Grafici salvati:\n")
  cat("  - Dataset1_Y_vs_X1_informative.png\n")
  cat("  - Dataset1_Y_vs_X5_noise.png\n\n")
  
} else {
  cat("(Installa ggplot2 per visualizzazioni)\n\n")
}

# ==============================================================================
# VERIFICA INDIPENDENZA NOISE
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat(" VERIFICA INDIPENDENZA COVARIATE NOISE DA Y\n")
cat("==============================================================================\n")
cat(" Test di correlazione per tutti i dataset\n\n")

correlations_summary <- matrix(NA, nrow = N_DATASETS, ncol = 7)
colnames(correlations_summary) <- paste0("X", 1:7)

for(d in 1:N_DATASETS) {
  for(j in 1:7) {
    correlations_summary[d, j] <- cor(DATA[[d]]$Y, DATA[[d]]$X[, j])
  }
}

cat(" Media correlazioni (su 10 dataset):\n")
cat(" -------------------------------------------\n")
mean_cors <- colMeans(correlations_summary)
for(j in 1:7) {
  is_informative <- j <= 4
  cat(sprintf(" X%d (%s): %6.3f (sd = %.3f)\n", 
              j,
              ifelse(is_informative, "INFO", "NOISE"),
              mean_cors[j],
              sd(correlations_summary[, j])))
}

cat("\n Le covariate noise (X5, X6, X7) dovrebbero avere\n")
cat(" correlazioni vicine a zero.\n")
cat("==============================================================================\n\n")

cat("GENERAZIONE COMPLETATA CON SUCCESSO!\n\n")
}


# Struttura:
# - 3 gruppi con dimensioni (75, 75, 50) = 200 osservazioni totali
# - 4 covariate informative: 2 continue (x1, x2) + 2 binarie (x3, x4)
# - 1 covariata noise: continue indipendenti da Y
# - Totale: 5 covariate (4 vere + 1 noise)
# ==============================================================================

if(FALSE){
  set.seed(12345)  # Per riproducibilità
  
  # ==============================================================================
  # PARAMETRI GLOBALI
  # ==============================================================================
  N_DATASETS <- 10      # Numero di dataset da generare
  N_TOTAL <- 100        # Osservazioni totali per dataset
  N_GROUPS <- 3         # Numero di gruppi veri
  GROUP_SIZES <- c(35, 40, 25)  # Dimensioni dei gruppi
  
  # ==============================================================================
  # PARAMETRI DEI GRUPPI (da Appendice E)
  # ==============================================================================
  
  # Gruppo 1
  mu_1 <- c(-3, 3)      # Medie per (x1, x2)
  q_1 <- 0.1            # Prob. per binarie (x3, x4)
  beta_1 <- c(1, 5, 2, 1, 0)  # Coefficienti: intercetta + 4 covariate
  
  # Gruppo 2
  mu_2 <- c(0, 0)
  q_2 <- 0.5
  beta_2 <- c(4, 2, -2, 1, -1)
  
  # Gruppo 3
  mu_3 <- c(3, 3)
  q_3 <- 0.9
  beta_3 <- c(-1, -5, -2, -1, 1)
  
  # Matrice di covarianza per covariate continue (identità scalata)
  Sigma_cov <- 0.5 * diag(2)
  
  # Varianza residua per Y
  sigma_y <- 0.5
  
  # Varianze per le covariate noise
  sigma_noise <- c(1, 2, 3)  # Varianze crescenti
  
  # ==============================================================================
  # FUNZIONE PER GENERARE UN SINGOLO DATASET
  # ==============================================================================
  
  genera_dataset_simulato <- function(dataset_id, verbose = TRUE) {
    
    if(verbose) {
      cat("\n========================================")
      cat("\n Generazione Dataset", dataset_id)
      cat("\n========================================\n")
    }
    
    # ----------------------------------------------------------------------------
    # 1. Inizializzazione
    # ----------------------------------------------------------------------------
    n_total <- sum(GROUP_SIZES)
    
    # Matrici di output
    X_all <- matrix(NA, nrow = n_total, ncol = 5)  # 4 vere + 1 noise
    #X_all <- matrix(NA, nrow = n_total, ncol = 7)  # 4 vere + 3 noise
    Y_all <- numeric(n_total)
    true_labels <- integer(n_total)
    
    colnames(X_all) <- paste0("X", 1:5)
    #colnames(X_all) <- paste0("X", 1:7)
    
    # ----------------------------------------------------------------------------
    # 2. Generazione per ogni gruppo
    # ----------------------------------------------------------------------------
    row_idx <- 1
    
    for(group in 1:N_GROUPS) {
      
      n_group <- GROUP_SIZES[group]
      
      # Selezione parametri del gruppo
      if(group == 1) {
        mu_g <- mu_1
        q_g <- q_1
        beta_g <- beta_1
      } else if(group == 2) {
        mu_g <- mu_2
        q_g <- q_2
        beta_g <- beta_2
      } else {
        mu_g <- mu_3
        q_g <- q_3
        beta_g <- beta_3
      }
      
      # ------------------------------------------------------------------------
      # 2a. Covariate CONTINUE INFORMATIVE (x1, x2)
      # ------------------------------------------------------------------------
      X_cont <- MASS::mvrnorm(n = n_group, mu = mu_g, Sigma = Sigma_cov)
      
      # ------------------------------------------------------------------------
      # 2b. Covariate BINARIE INFORMATIVE (x3, x4)
      # ------------------------------------------------------------------------
      X_bin <- matrix(
        rbinom(n = 2 * n_group, size = 1, prob = q_g),
        nrow = n_group,
        ncol = 2
      )
      
      # ------------------------------------------------------------------------
      # 2c. Covariate NOISE (x5, x6, x7) - INDIPENDENTI DA Y
      # ------------------------------------------------------------------------
      X_noise <- cbind(
        rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[1]))  # x5 ~ N(0, 1)
        #rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[2])),  # x6 ~ N(0, 2)
        #rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[3]))   # x7 ~ N(0, 3)
      )
      
      # ------------------------------------------------------------------------
      # 2d. Combina tutte le covariate
      # ------------------------------------------------------------------------
      X_group <- cbind(X_cont, X_bin, X_noise)
      
      # ------------------------------------------------------------------------
      # 2e. Genera risposta Y (dipende SOLO da x1, x2, x3, x4)
      # ------------------------------------------------------------------------
      # Design matrix: intercetta + 4 covariate informative
      X_design <- cbind(1, X_cont, X_bin)  # 5 colonne
      
      # Media condizionale
      mu_y <- X_design %*% beta_g
      
      # Risposta con errore gaussiano
      Y_group <- rnorm(n_group, mean = mu_y, sd = sqrt(sigma_y))
      
      # ------------------------------------------------------------------------
      # 2f. Salva nel dataset completo
      # ------------------------------------------------------------------------
      idx_range <- row_idx:(row_idx + n_group - 1)
      X_all[idx_range, ] <- X_group
      Y_all[idx_range] <- Y_group
      true_labels[idx_range] <- group
      
      row_idx <- row_idx + n_group
      
      if(verbose) {
        cat("  Gruppo", group, "| N =", n_group, "| rows:", 
            min(idx_range), "-", max(idx_range), "\n")
      }
    }
    
    # ----------------------------------------------------------------------------
    # 3. Permutazione casuale (importante per MCMC)
    # ----------------------------------------------------------------------------
    perm_idx <- sample(1:n_total)
    X_all <- X_all[perm_idx, ]
    Y_all <- Y_all[perm_idx]
    true_labels <- true_labels[perm_idx]
    
    # ----------------------------------------------------------------------------
    # 4. Output
    # ----------------------------------------------------------------------------
    dataset <- list(
      X = X_all,
      Y = Y_all,
      true_clusters = true_labels,
      selected.cols = 1:4,  # Le prime 4 covariate sono quelle vere
      noise.cols = 5,       # Solo la 5 è noise (modificato da 5:7)
      group_params = list(
        mu = list(mu_1, mu_2, mu_3),
        q = c(q_1, q_2, q_3),
        beta = list(beta_1, beta_2, beta_3)
      ),
      info = list(
        n_total = n_total,
        n_groups = N_GROUPS,
        group_sizes = GROUP_SIZES,
        dataset_id = dataset_id
      )
    )
    
    if(verbose) {
      cat("\n  Dataset generato:")
      cat("\n  - Dimensioni X:", nrow(X_all), "x", ncol(X_all))
      cat("\n  - Covariate vere: 1, 2, 3, 4")
      cat("\n  - Covariate noise: 5")
      cat("\n  - Gruppi veri:", unique(true_labels))
      cat("\n")
    }
    
    return(dataset)
  }
  
  # ==============================================================================
  # GENERAZIONE DI TUTTI I DATASET
  # ==============================================================================
  
  cat("\n")
  cat("==============================================================================\n")
  cat(" GENERAZIONE", N_DATASETS, "DATASET SIMULATI\n")
  cat("==============================================================================\n")
  cat(" Struttura per dataset:\n")
  cat("  - N osservazioni:", N_TOTAL, "\n")
  cat("  - N gruppi veri:", N_GROUPS, "\n")
  cat("  - Dimensioni gruppi:", paste(GROUP_SIZES, collapse = ", "), "\n")
  cat("  - Covariate informative: 4 (2 continue + 2 binarie)\n")
  cat("  - Covariate noise: 1 (indipendenti da Y)\n") # Modificato commento logico
  cat("  - Totale covariate: 5\n") # Modificato commento logico
  cat("==============================================================================\n")
  
  DATA <- vector("list", N_DATASETS)
  
  for(d in 1:N_DATASETS) {
    DATA[[d]] <- genera_dataset_simulato(dataset_id = d, verbose = TRUE)
  }
  
  # ==============================================================================
  # DIAGNOSTICA RAPIDA
  # ==============================================================================
  
  cat("\n")
  cat("==============================================================================\n")
  cat(" DIAGNOSTICA DATASET\n")
  cat("==============================================================================\n")
  
  # Esempio: primo dataset
  d1 <- DATA[[1]]
  
  cat("\n Dataset 1 - Statistiche descrittive:\n")
  cat(" -------------------------------------------\n")
  cat(" Y: mean =", round(mean(d1$Y), 2), 
      ", sd =", round(sd(d1$Y), 2), "\n")
  
  for(j in 1:5) { # Modificato da 1:7
    is_informative <- j %in% d1$selected.cols
    cat(sprintf(" X%d (%s): mean = %6.2f, sd = %6.2f\n", 
                j, 
                ifelse(is_informative, "INFO", "NOISE"),
                mean(d1$X[, j]),
                sd(d1$X[, j])))
  }
  
  cat("\n Correlazioni Y con covariate:\n")
  cat(" -------------------------------------------\n")
  for(j in 1:5) { # Modificato da 1:7
    is_informative <- j %in% d1$selected.cols
    cor_yj <- cor(d1$Y, d1$X[, j])
    cat(sprintf(" cor(Y, X%d) = %6.3f %s\n", 
                j, 
                cor_yj,
                ifelse(is_informative, "✓", " ")))
  }
  
  # ==============================================================================
  # SALVATAGGIO
  # ==============================================================================
  
  output_filename <- "Dataset_Simulati_AppendixE_WithNoise.RData"
  save(DATA, file = output_filename)
  
  # ==============================================================================
  # VISUALIZZAZIONE ESEMPIO (opzionale)
  # ==============================================================================
  
  if(require(ggplot2, quietly = TRUE)) {
    
    cat("Creazione grafici diagnostici...\n")
    
    # Dataset 1
    d1 <- DATA[[1]]
    df_plot <- data.frame(
      Y = d1$Y,
      X1 = d1$X[, 1],
      X2 = d1$X[, 2],
      X5_noise = d1$X[, 5],
      Cluster = factor(d1$true_clusters)
    )
    
    # Plot 1: Y vs X1 (informativa)
    p1 <- ggplot(df_plot, aes(x = X1, y = Y, color = Cluster)) +
      geom_point(alpha = 0.6, size = 2) +
      theme_minimal() +
      labs(title = "Dataset 1: Y vs X1 (informativa)",
           subtitle = "I 3 gruppi veri sono visibili") +
      theme(legend.position = "bottom")
    
    # Plot 2: Y vs X5 (noise)
    p2 <- ggplot(df_plot, aes(x = X5_noise, y = Y, color = Cluster)) +
      geom_point(alpha = 0.6, size = 2) +
      theme_minimal() +
      labs(title = "Dataset 1: Y vs X5 (noise)",
           subtitle = "Nessuna relazione con Y") +
      theme(legend.position = "bottom")
    
    # Salva plots
    ggsave("Dataset1_Y_vs_X1_informative.png", plot = p1, width = 8, height = 6)
    ggsave("Dataset1_Y_vs_X5_noise.png", plot = p2, width = 8, height = 6)
    
  }
  
  # ==============================================================================
  # VERIFICA INDIPENDENZA NOISE
  # ==============================================================================
  
  cat("\n")
  cat("==============================================================================\n")
  cat(" VERIFICA INDIPENDENZA COVARIATE NOISE DA Y\n")
  cat("==============================================================================\n")
  
  correlations_summary <- matrix(NA, nrow = N_DATASETS, ncol = 5) # Modificato da 7 a 5
  colnames(correlations_summary) <- paste0("X", 1:5)
  
  for(d in 1:N_DATASETS) {
    for(j in 1:5) {
      correlations_summary[d, j] <- cor(DATA[[d]]$Y, DATA[[d]]$X[, j])
    }
  }
  
  cat(" Media correlazioni (su 10 dataset):\n")
  cat(" -------------------------------------------\n")
  mean_cors <- colMeans(correlations_summary)
  for(j in 1:5) {
    is_informative <- j <= 4
    cat(sprintf(" X%d (%s): %6.3f (sd = %.3f)\n", 
                j,
                ifelse(is_informative, "INFO", "NOISE"),
                mean_cors[j],
                sd(correlations_summary[, j])))
  }
  
  cat("\n Le covariate noise (X5) dovrebbero avere\n")
  cat(" correlazioni vicine a zero.\n")
  cat("==============================================================================\n\n")
  
  cat("GENERAZIONE COMPLETATA CON SUCCESSO!\n\n")
}