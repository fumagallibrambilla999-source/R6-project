# -----------------------------------------------------------------
# FUNZIONE HELPER PER IL CALCOLO DELLA LOG-DENSITÀ PREDITTIVA
# -----------------------------------------------------------------
# Questa funzione è il "motore" del clustering. Calcola quanto un'osservazione
# è probabile all'interno di un cluster integrando i parametri (Marginal Likelihood).
# Il risultato è una distribuzione t-Student non standardizzata.

calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  
  # Estrazione iperparametri dalla Misura di Base P0 (Normal-Gamma)
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  
  njj <- length(z_j)
  
  if (njj == 0) {
    # --- CASO A: PREDITTIVA A PRIORI (Nuovo Cluster) ---
    # Se il cluster non esiste, la predittiva dipende solo dalla prior P0.
    df_loc  <- 2 * alpha0
    mu_loc  <- mu0
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
    
  } else {
    # --- CASO B: PREDITTIVA POSTERIORE (Cluster Esistente) ---
    # Aggiornamento Bayesiano basato sulle osservazioni già presenti nel cluster.
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    
    # Somma degli scarti quadratici (Sum of Squares)
    ss_j <- sum(z_j^2) - (sum_z^2) / njj 
    
    # Aggiornamento del parametro di scala (beta) della Gamma
    beta_loc  <- beta0 + 
      0.5 * ss_j + 
      0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    
    df_loc  <- 2 * alpha_loc
    # Varianza della t-Student predittiva
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) / (alpha_loc * kappa_loc)
  }
  
  # Controllo robustezza: evita errori numerici se la scala è nulla
  if (sig_sq_loc <= 0) {
    warning("Scala predittiva non positiva, restituisco -Inf")
    return(-Inf)
  }
  
  sig_loc <- sqrt(sig_sq_loc)
  
  # Calcolo della densità log-t: log(p(z_i | z_j))
  # Usiamo log = TRUE per mantenere la precisione numerica (evita underflow)
  log_dens <- dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc)
  
  return(log_dens)
}

# -----------------------------------------------------------------
# MODELLO PPMx (Product Partition Model with Covariates)
# -----------------------------------------------------------------
# Il modello PPMx estende il DPMM permettendo alle covariate X di influenzare
# direttamente la probabilità che due osservazioni finiscano nello stesso cluster.
# -----------------------------------------------------------------

#' Implementazione di un PPMx con kernel Gaussiano (Algoritmo 8 Neal)
#'
#' @param y Vettore della variabile risposta.
#' @param x Matrice delle covariate (n x P).
#' @param n_iter Numero totale di iterazioni MCMC.
#' @param burn_in Iterazioni scartate (fase di riscaldamento).
#' @param thin Fattore di diradamento (salvataggio periodico).
#' @param hyperpar_P0_y Iperparametri per la risposta Y.
#' @param hyperpar_P0_x Lista di liste di iperparametri (una per ogni colonna di X).
#' @param hyperpar_Dir Iperparametro alpha (concentrazione DP).
#' @param init Stato iniziale della catena.

PPMx_kernel_gaussiano <- function(y, x, 
                                  n_iter, burn_in, thin,
                                  hyperpar_P0_y, hyperpar_P0_x, 
                                  hyperpar_Dir, init) {
  
  # --- 1. SETUP E VALIDAZIONE ---
  n <- length(y)
  if(n == 0) stop("y deve contenere dati.")
  
  if(is.null(x)) stop("Matrice x necessaria per il modello PPMx.")
  if(!is.matrix(x)) x <- as.matrix(x)
  
  n_covariate <- ncol(x) 
  
  saved_iters <- seq(from = burn_in + 1, to = n_iter, by = thin)
  G <- length(saved_iters)
  
  # --- 2. PRE-ALLOCAZIONE OUTPUT ---
  k_out   <- rep(NA_integer_, G)
  clu_out <- matrix(NA_integer_, nrow = G, ncol = n)
  nj_out  <- vector("list", G)
  tauj_out <- vector("list", G)
  muj_out  <- vector("list", G)
  
  # --- 3. INIZIALIZZAZIONE ---
  k_current    <- init$k
  clu_current  <- as.integer(init$clu)
  nj_current   <- as.integer(init$nj)
  tauj_current <- as.numeric(init$tauj)
  muj_current  <- as.numeric(init$muj)
  
  # Struttura dati 'rho' per accesso rapido agli indici di ogni cluster
  rho <- vector("list", k_current)
  for(j in seq_len(k_current)) rho[[j]] <- which(clu_current == j)
  
  # Verifica integrità dello stato iniziale
  sanitize_state <- function() {
    if(length(nj_current) != k_current) stop("Mismatch lunghezza nj_current.")
    if(any(sapply(rho, length) != nj_current)) stop("Incoerenza tra rho e nj_current.")
  }
  sanitize_state()
  
  save_idx <- 0
  
  # --- 4. LOOP GIBBS SAMPLER ---
  for(g in 1:n_iter) {
    
    # --- 4a. AGGIORNAMENTO ALLOCAZIONI (c_i) ---
    for(i in 1:n) {
      
      # 1. Rimozione dell'osservazione i dal cluster attuale
      clu_i <- clu_current[i]
      idx_in_rho <- which(rho[[clu_i]] == i)
      rho[[clu_i]] <- rho[[clu_i]][-idx_in_rho]
      nj_current[clu_i] <- nj_current[clu_i] - 1
      
      # 2. Se il cluster rimane vuoto, lo eliminiamo (Housekeeping)
      if(nj_current[clu_i] == 0) {
        nj_current <- nj_current[-clu_i]
        rho <- rho[-clu_i]
        k_current <- k_current - 1
        # Ri-etichettiamo i cluster successivi per evitare "buchi" negli indici
        clu_current[clu_current > clu_i] <- clu_current[clu_current > clu_i] - 1
      }
      
      # 3. Calcolo LOG-PESI (Log-Probability) per il ricampionamento
      # Nel PPMx, il peso è: Coesione * Likelihood_Y * Similarity_X
      log_w <- numeric(k_current + 1)
      
      if(k_current > 0) {
        for(j in 1:k_current) {
          # (1) Prior di Coesione (Polya Urn / CRP)
          log_w[j] <- log(nj_current[j])
          
          # (2) Predittiva per Y (quanto i è simile a Y nel cluster j)
          log_w[j] <- log_w[j] + calcola_log_pred_dens_t(y[i], y[rho[[j]]], hyperpar_P0_y)
          
          # (3) Similarità per X (Prodotto delle predittive per ogni covariata)
          # 
          for (l in 1:n_covariate) {
            log_w[j] <- log_w[j] + calcola_log_pred_dens_t(x[i, l], x[rho[[j]], l], hyperpar_P0_x[[l]])
          }
        }
      }
      
      # Peso per la creazione di un NUOVO cluster (basato su alpha e P0)
      log_w[k_current + 1] <- log(hyperpar_Dir$alpha) + 
        calcola_log_pred_dens_t(y[i], numeric(0), hyperpar_P0_y)
      for (l in 1:n_covariate) {
        log_w[k_current + 1] <- log_w[k_current + 1] + 
          calcola_log_pred_dens_t(x[i, l], numeric(0), hyperpar_P0_x[[l]])
      }
      
      # --- 4a.4 Trasformazione da Log-Space a Probabilità ---
      # Usiamo il Log-Sum-Exp trick sottraendo il max per stabilità numerica
      max_log_w <- max(log_w)
      w <- exp(log_w - max_log_w)
      
      new_cluster <- sample.int(length(w), size = 1L, prob = w)
      clu_current[i] <- new_cluster
      
      # --- 4a.5 Aggiornamento Stato Corrente ---
      if(new_cluster <= k_current) {
        nj_current[new_cluster] <- nj_current[new_cluster] + 1
        rho[[new_cluster]] <- c(rho[[new_cluster]], i)
      } else {
        # Creazione nuovo cluster
        nj_current <- c(nj_current, 1L) 
        k_current <- k_current + 1 
        rho[[k_current]] <- i 
        # Inizializzazione placeholder per i parametri di Y
        tauj_current <- c(tauj_current, 1.0) 
        muj_current <- c(muj_current, 0.0)   
      }
    } # Fine loop osservazioni i
    
    # --- 4b. CAMPIONAMENTO PARAMETRI DI Y (mu, tau) ---
    # Dato il clustering, aggiorniamo i parametri della risposta per ogni gruppo.
    # 
    if(k_current > 0) {
      for(j in 1:k_current) {
        idx_j <- rho[[j]]
        sum_y  <- sum(y[idx_j])
        sum_y2 <- sum(y[idx_j]^2)
        njj    <- nj_current[j]
        
        kn <- hyperpar_P0_y$kappa0 + njj
        an <- hyperpar_P0_y$nu0/2 + njj/2
        mun <- (hyperpar_P0_y$kappa0 * hyperpar_P0_y$mu0 + sum_y) / kn
        bn <- 0.5 * (hyperpar_P0_y$nu0 * hyperpar_P0_y$sig_square0) +
          0.5 * (sum_y2 - (sum_y^2)/njj) +
          0.5 * (hyperpar_P0_y$kappa0 * njj / kn) * (sum_y/njj - hyperpar_P0_y$mu0)^2
        
        tauj_current[j] <- rgamma(1, shape = an, rate = bn)
        muj_current[j]  <- rnorm(1, mean = mun, sd = sqrt(1 / (kn * tauj_current[j])))
      }
    }
    
    # --- 4c. SALVATAGGIO DEI RISULTATI ---
    if(g %in% saved_iters) {
      save_idx <- save_idx + 1
      k_out[save_idx] <- k_current
      clu_out[save_idx, ] <- clu_current
      nj_out[[save_idx]] <- nj_current
      tauj_out[[save_idx]] <- tauj_current[1:k_current]
      muj_out[[save_idx]] <- muj_current[1:k_current]
    }
    
    if(g %% 500 == 0) sanitize_state() # Check periodico di coerenza
  } 
  
  return(list(
    k_out = k_out, clu_out = clu_out, nj_out = nj_out,
    tauj_out = tauj_out, muj_out = muj_out, saved_iters = saved_iters
  ))
}