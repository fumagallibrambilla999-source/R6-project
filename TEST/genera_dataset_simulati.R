# ==============================================================================
# GENERAZIONE DATASET SIMULATI - Modello PPMx con Variable Selection
# ==============================================================================
# Questo script genera dataset sintetici complessi per testare la capacità del 
# modello di identificare correttamente le covariate informative rispetto al noise.
#
# Struttura del segnale:
# - 3 popolazioni (cluster) con parametri di regressione distinti.
# - Covariate X1, X2: Continue, definiscono la separazione spaziale.
# - Covariate X3, X4: Binarie, definiscono la similarità categoriale.
# - Covariate X5, X6, X7: Noise (disturbatori), indipendenti dalla risposta Y.
# ==============================================================================

if(TRUE){
  set.seed(12345)  # Garantisce la riproducibilità degli esperimenti
  
  # --- PARAMETRI DI CONFIGURAZIONE ---
  N_DATASETS <- 10      
  N_TOTAL <- 200        
  N_GROUPS <- 3         
  GROUP_SIZES <- c(75, 75, 50)  
  
  # --- PARAMETRI DEI GRUPPI (Specifiche Appendice E) ---
  # Ogni gruppo ha una propria "logica" che lega X a Y tramite coefficienti beta.
  
  # Gruppo 1: Medie basse, probabilità binarie basse
  mu_1 <- c(-3, 3)      
  q_1 <- 0.1            
  beta_1 <- c(1, 5, 2, 1, 0)  
  
  # Gruppo 2: Medie centrali, probabilità binarie bilanciate
  mu_2 <- c(0, 0)
  q_2 <- 0.5
  beta_2 <- c(4, 2, -2, 1, -1)
  
  # Gruppo 3: Medie alte, probabilità binarie elevate
  mu_3 <- c(3, 3)
  q_3 <- 0.9
  beta_3 <- c(-1, -5, -2, -1, 1)
  
  Sigma_cov <- 0.5 * diag(2) # Covarianza intra-cluster per X continue
  sigma_y <- 0.5             # Incertezza sulla risposta (Rumore di misura)
  sigma_noise <- c(1, 2, 3)  # Diversi livelli di varianza per le noise variables
  
  # --- FUNZIONE CORE DI GENERAZIONE ---
  genera_dataset_simulato <- function(dataset_id, verbose = TRUE) {
    
    n_total <- sum(GROUP_SIZES)
    X_all <- matrix(NA, nrow = n_total, ncol = 7)  
    Y_all <- numeric(n_total)
    true_labels <- integer(n_total)
    colnames(X_all) <- paste0("X", 1:7)
    
    row_idx <- 1
    for(group in 1:N_GROUPS) {
      n_group <- GROUP_SIZES[group]
      
      # Assegnazione parametri dinamica in base al gruppo
      if(group == 1) { mu_g <- mu_1; q_g <- q_1; beta_g <- beta_1 }
      else if(group == 2) { mu_g <- mu_2; q_g <- q_2; beta_g <- beta_2 }
      else { mu_g <- mu_3; q_g <- q_3; beta_g <- beta_3 }
      
      # 1. Generazione Covariate INFORMATIVE (Continue e Binarie)
      X_cont <- MASS::mvrnorm(n = n_group, mu = mu_g, Sigma = Sigma_cov)
      X_bin <- matrix(rbinom(n = 2 * n_group, size = 1, prob = q_g), nrow = n_group, ncol = 2)
      
      # 2. Generazione Covariate NOISE (X5, X6, X7)
      # Queste variabili NON entrano nel calcolo di Y, testano la Variable Selection.
      X_noise <- cbind(
        rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[1])),
        rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[2])),
        rnorm(n_group, mean = 0, sd = sqrt(sigma_noise[3]))
      )
      
      X_group <- cbind(X_cont, X_bin, X_noise)
      
      # 3. Generazione Risposta Y (Likelihood)
      # Segue il modello lineare Y = Intercetta + X*Beta + epsilon
      X_design <- cbind(1, X_cont, X_bin) # Matrice di disegno (include intercetta)
      mu_y <- X_design %*% beta_g
      Y_group <- rnorm(n_group, mean = mu_y, sd = sqrt(sigma_y))
      
      # Archiviazione dati
      idx_range <- row_idx:(row_idx + n_group - 1)
      X_all[idx_range, ] <- X_group
      Y_all[idx_range] <- Y_group
      true_labels[idx_range] <- group
      row_idx <- row_idx + n_group
    }
    
    # 4. PERMUTAZIONE CASUALE
    # Essenziale per evitare che l'MCMC sia influenzato dall'ordine dei dati
    perm_idx <- sample(1:n_total)
    X_all <- X_all[perm_idx, ]
    Y_all <- Y_all[perm_idx]
    true_labels <- true_labels[perm_idx]
    
    return(list(X = X_all, Y = Y_all, true_clusters = true_labels, 
                selected.cols = 1:4, noise.cols = 5:7))
  }
  
  # --- ESECUZIONE GENERAZIONE ---
  DATA <- lapply(1:N_DATASETS, genera_dataset_simulato)
  
  # --- DIAGNOSTICA ---
  # Verifica che le covariate informative siano effettivamente correlate con Y
  # mentre le noise abbiano correlazione prossima allo zero.
  cat("\nVerifica Correlazioni Media (Segnale vs Rumore):\n")
  cor_results <- sapply(DATA, function(d) cor(d$Y, d$X))
  print(rowMeans(cor_results))
  
  # Salva l'oggetto DATA per l'utilizzo nel modello MCMC
  save(DATA, file = "Dataset_Simulati_AppendixE_WithNoise.RData")
}
  