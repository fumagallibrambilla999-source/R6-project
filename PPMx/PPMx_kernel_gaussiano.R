# -----------------------------------------------------------------
# FUNZIONE HELPER PER CALCOLO DENSITÀ PREDITTIVA
# -----------------------------------------------------------------

#' Calcola la log-densità predittiva t-Student
#'
#' Calcola la log-densità per una nuova osservazione 'z_i' data
#' un set di osservazioni 'z_j' da un cluster, usando una prior
#' Normal-Gamma (parametrizzata come Normal-InvGamma).
#'
#' @param z_i L'osservazione singola da valutare.
#' @param z_j Un vettore di osservazioni *già* nel cluster (può essere numeric(0)).
#' @param hyperpar_P0 Una lista con (mu0, kappa0, nu0, sig_square0).
#'
#' @return Il valore della log-densità (scalare).

calcola_log_pred_dens_t <- function(z_i, z_j, hyperpar_P0) {
  
  # Estrai iperparametri
  mu0     <- hyperpar_P0$mu0
  kappa0  <- hyperpar_P0$kappa0
  nu0     <- hyperpar_P0$nu0
  sig_sq0 <- hyperpar_P0$sig_square0
  
  # Calcola i parametri prior per la t-Student (per il caso "nuovo cluster")
  alpha0 <- nu0 / 2
  beta0  <- 0.5 * (nu0 * sig_sq0)
  
  # Statistiche sufficienti
  njj <- length(z_j)
  
  if (njj == 0) {
    # --- CASO A: PREDITTIVA A PRIORI (Nuovo Cluster) ---
    
    # Parametri t-Student (basati solo sulla prior)
    df_loc  <- 2 * alpha0
    mu_loc  <- mu0
    # Scala^2
    sig_sq_loc <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
    
  } else {
    # --- CASO B: PREDITTIVA POSTERIORE (Cluster Esistente) ---
    
    # Statistiche sufficienti
    sum_z  <- sum(z_j)
    mean_z <- sum_z / njj
    
    # Iperparametri "locali" (aggiornati) per la Normal-Gamma
    kappa_loc <- kappa0 + njj
    mu_loc    <- (kappa0 * mu0 + sum_z) / kappa_loc
    alpha_loc <- alpha0 + njj / 2
    
    # Calcolo di beta_loc
    # Usiamo sum( (z_j - mean_z)^2 )
    ss_j <- sum(z_j^2) - (sum_z^2) / njj 
    beta_loc  <- beta0 + 
      0.5 * ss_j + 
      0.5 * (kappa0 * njj / kappa_loc) * (mean_z - mu0)^2
    
    # Parametri t-Student (basati sulla posterior)
    df_loc  <- 2 * alpha_loc
    # Scala^2
    sig_sq_loc <- (beta_loc * (kappa_loc + 1)) / (alpha_loc * kappa_loc)
  }
  
  # Evita divisioni per zero o scale negative se i parametri sono degeneri
  if (is.na(sig_sq_loc) || sig_sq_loc <= 0) {
    warning("Scala predittiva non positiva, restituisco -Inf")
    return(-Inf)
  }
  
  sig_loc <- sqrt(sig_sq_loc)
  
  # Calcola la log-densità t-Student non-standardizzata
  # log( (1/sig_loc) * dt( (z_i - mu_loc) / sig_loc, df = df_loc) )
  log_dens <- dt((z_i - mu_loc) / sig_loc, df = df_loc, log = TRUE) - log(sig_loc)
  
  return(log_dens)
}


# -----------------------------------------------------------------
# -----------------------------------------------------------------
# MODELLO PPMx (Product Partition Model with Covariates)
# -----------------------------------------------------------------
# -----------------------------------------------------------------

# 1. LIVELLO DATI (LIKELIHOOD per Y)
# Y_i | theta_i ~ N(mu_i, tau_i^-1)
# theta_i = (mu_i, tau_i)

# 2. LIVELLO PROCESSO (Allocazioni)
# I parametri theta_i sono estratti da una distribuzione (discreta) P.
# P(c_i = j | c_-i, x) \propto c(S_j) * g(x_j^*) [formula del paper]
#
# Questo è implementato (via Polya Urn) come:
# P(c_i = j) \propto P(c_i=j | c_-i) * P(y_i | y_j) * P(x_i | x_j)
#
# P(c_i = j) \propto (n_j) * m(y_i | y_j) * g(x_i | x_j)
#
# - n_j: Coesione (Polya Urn, DP)
# - m(y_i | y_j): Predictive likelihood per Y (t-Student)
# [cite_start]- g(x_i | x_j): Similarità (Predictive likelihood per X) [cite: 681]

# 3. LIVELLO PRIOR (PROCESSO DI DIRICHLET)
# P ~ DP(alpha, P0)
# - alpha: è il parametro di concentrazione (scalare).
# - P0: è la misura di base (Normal-Gamma)

# [cite_start]4. MODELLO AUSILIARIO PER X (Similarità) [cite: 687]
# Per ogni covariata l=1...P:
# X_il | eta_il ~ N(mu_il, tau_il^-1)
# eta_il ~ P0_x,l
# P0_x,l = Normal-Gamma(...)
#
# La similarità g(x_j*) è il prodotto delle predittive (t-Student)
# per ciascuna covariata:
# [cite_start]g(x_j*) = \prod_{l=1}^P g_l(x_jl*) [cite: 1036]

# -----------------------------------------------------------------


#' Implementazione di un PPMx con kernel Gaussiano (Algoritmo 8 Neal)
#'
#' @param y Vettore (lungo n) delle risposte
#' @param x Matrice (n x P) delle covariate. P = numero di covariate.
#' @param n_iter Numero totale iterazioni MCMC
#' @param burn_in Numero iterazioni da scartare
#' @param thin Intervallo di salvataggio
#' @param hyperpar_P0_y Lista iperparametri per Y: (mu0, kappa0, nu0, sig_square0)
#' @param hyperpar_P0_x LISTA di liste. Lunga P.
#'                      hyperpar_P0_x[[l]] = lista iperpar. per la covariata l.
#' @param hyperpar_Dir Lista contenente il parametro 'alpha' del DP.
#' @param init Lista con lo stato iniziale (k, clu, nj, tauj, muj).
#'
#' @return Lista di output MCMC (k_out, clu_out, nj_out, tauj_out, muj_out)

PPMx_kernel_gaussiano <- function(y, x, 
                                  n_iter, burn_in, thin,
                                  hyperpar_P0_y, hyperpar_P0_x, 
                                  hyperpar_Dir, init) {
  
  # --- 1. SETUP INIZIALE E VALIDAZIONE INPUT ---
  
  # n = numero totale di osservazioni
  n <- length(y)
  if(n == 0) stop("y must be non-empty")
  
  # Gestione covariate
  if(is.null(x)) stop("x (covariate) matrix must be provided.")
  if(!is.matrix(x)) x <- as.matrix(x)
  if(nrow(x) != n) stop("Number of rows in x must match length of y.")
  
  n_covariate <- ncol(x) # Numero di covariate (P)
  if(length(hyperpar_P0_x) != n_covariate) {
    stop("hyperpar_P0_x must be a list of length equal to ncol(x).")
  }
  
  if(burn_in >= n_iter) stop("burn_in must be < n_iter")
  if(thin <= 0) stop("thin must be positive")
  
  saved_iters <- seq(from = burn_in + 1, to = n_iter, by = thin)
  G <- length(saved_iters)
  if(G == 0) stop("No iterations will be saved: check burn_in/thin/n_iter")
  
  # --- 2. PRE-ALLOCAZIONE DEGLI OGGETTI DI OUTPUT ---
  # (Output si riferisce ai parametri per Y, non X)
  k_out   <- rep(NA_integer_, G)
  clu_out <- matrix(NA_integer_, nrow = G, ncol = n)
  nj_out  <- vector("list", G)
  tauj_out<- vector("list", G)
  muj_out <- vector("list", G)
  
  # --- 3. INIZIALIZZAZIONE DELLO STATO CORRENTE ---
  k_current   <- init$k
  clu_current <- as.integer(init$clu)
  nj_current  <- as.integer(init$nj)
  tauj_current<- as.numeric(init$tauj)
  muj_current <- as.numeric(init$muj)
  
  # 'rho' è la lista (lunga k_current) degli INDICI
  rho <- vector("list", k_current)
  for(j in seq_len(k_current)) rho[[j]] <- which(clu_current == j)
  
  # Funzione di utilità per controllo coerenza
  sanitize_state <- function() {
    if(length(nj_current) != k_current) stop("nj_current length mismatch k_current")
    if(length(rho) != k_current) stop("rho length mismatch k_current")
    if(any(sapply(rho, length) != nj_current)) stop("rho and nj_current inconsistent")
    if(any(unlist(rho) < 1 | unlist(rho) > n)) stop("rho contains invalid indices")
  }
  sanitize_state()
  
  # Indice di salvataggio
  save_idx <- 0
  
  # --- 4. INIZIO DEL LOOP MCMC (GIBBS SAMPLER) ---
  
  cat("   -> Loop MCMC avviato...\n")
  pb <- progress::progress_bar$new( # <--- Aggiunto progress::
    format = "  Iterazione MCMC [:bar] :percent (ETA: :eta)",
    total = n_iter,
    clear = FALSE,
    width = 70
  )
  # --- FINE AGGIUNTA ---
  
  for(g in 1:n_iter) {
    
    # --- 4a. CAMPIONAMENTO DELLE ALLOCAZIONI (c_i) ---
    for(i in 1:n) {
      
      # --- 4a.1 Rimuovi l'osservazione 'i' dal suo cluster ---
      clu_i <- clu_current[i]
      idx_in_rho <- which(rho[[clu_i]] == i)
      if(length(idx_in_rho) != 1) stop("index_to_remove problem")
      rho[[clu_i]] <- rho[[clu_i]][-idx_in_rho]
      nj_current[clu_i] <- nj_current[clu_i] - 1
      
      # Estrai i dati dell'osservazione 'i'
      y_i <- y[i]
      x_i <- x[i, ] # Vettore lungo n_covariate
      
      # --- 4a.2 Gestione cluster vuoti ---
      if(nj_current[clu_i] == 0) {
        nj_current <- nj_current[-clu_i]
        rho <- rho[-clu_i]
        k_current <- k_current - 1
        clu_current[clu_current > clu_i] <- clu_current[clu_current > clu_i] - 1
      }
      
      # --- 4a.3 Calcolo Pesi per l'Allocazione (in LOG-SPACE) ---
      # w_j \propto n_j * m(y_i | y_j) * \prod_l g_l(x_il | x_jl)
      # log(w_j) = log(n_j) + log(m_y) + sum( log(g_l) )
      
      log_w <- numeric(k_current + 1)
      
      # Calcola i pesi per i k_current cluster *esistenti*
      if(k_current > 0) {
        for(j in 1:k_current) {
          
          # (1) Prior (Coesione DP)
          log_prior_j <- log(nj_current[j])
          
          # Indici delle osservazioni nel cluster j
          rho_current <- rho[[j]] 
          
          # (2) Predictive Likelihood per Y
          y_j <- y[rho_current]
          log_pred_y_j <- calcola_log_pred_dens_t(y_i, y_j, hyperpar_P0_y)
          
          # (3) Predictive Likelihood per X (Similarità)
          log_pred_x_j <- 0.0
          for (l in 1:n_covariate) {
            x_lj <- x[rho_current, l]
            x_il <- x_i[l]
            log_pred_x_j <- log_pred_x_j + 
              calcola_log_pred_dens_t(x_il, x_lj, hyperpar_P0_x[[l]])
          }
          
          # (4) Peso totale
          log_w[j] <- log_prior_j + log_pred_y_j + log_pred_x_j
        }
      }
      
      # --- Calcolo peso per NUOVO cluster (j = k_current + 1) ---
      # log(w_new) = log(alpha) + log(m_y(y_i)) + sum( log(g_l(x_il)) )
      
      # (1) Prior (Coesione DP)
      log_prior_new <- log(hyperpar_Dir$alpha)
      
      # (2) Predictive Likelihood per Y (a priori)
      log_pred_y_new <- calcola_log_pred_dens_t(y_i, numeric(0), hyperpar_P0_y)
      
      # (3) Predictive Likelihood per X (a priori)
      log_pred_x_new <- 0.0
      for (l in 1:n_covariate) {
        x_il <- x_i[l]
        log_pred_x_new <- log_pred_x_new + 
          calcola_log_pred_dens_t(x_il, numeric(0), hyperpar_P0_x[[l]])
      }
      
      # (4) Peso totale
      log_w[k_current + 1] <- log_prior_new + log_pred_y_new + log_pred_x_new
      
      # --- 4a.4 Campiona e assegna il nuovo cluster ---
      
      # Log-Sum-Exp trick per stabilità numerica
      # (evita exp(log_w) se i log_w sono molto piccoli)
      max_log_w <- max(log_w)
      w <- exp(log_w - max_log_w)
      
      # Controllo di stabilità
      if(all(w == 0)) w <- rep(1, length(w)) # Assegna a caso se tutti 0
      
      new_cluster <- sample.int(length(w), size = 1L, prob = w)
      
      clu_current[i] <- new_cluster
      
      # --- 4a.5 Aggiorna lo stato ---
      if(new_cluster <= k_current) {
        # Caso 1: Cluster esistente
        nj_current[new_cluster] <- nj_current[new_cluster] + 1
        rho[[new_cluster]] <- c(rho[[new_cluster]], i)
      } else {
        # Caso 2: Nuovo cluster
        nj_current <- c(nj_current, 1L) 
        k_current <- k_current + 1 
        rho[[k_current]] <- i 
        
        # Aggiungi placeholder per i parametri di Y
        # (verranno campionati correttamente in 4b)
        tauj_current <- c(tauj_current, 1.0) 
        muj_current <- c(muj_current, 0.0)   
      }
    } # --- Fine loop su i (osservazioni) ---
    
    
    # --- 4b. CAMPIONAMENTO PARAMETRI CLUSTER (per Y) ---
    # Questa parte è IDENTICA a prima.
    # Stiamo campionando i parametri P(theta_j | c, y)
    # [cite_start]I parametri per X sono ausiliari e non ci servono. [cite: 742]
    
    if(k_current > 0) {
      for(j in 1:k_current) {
        
        rho_current <- rho[[j]]
        sum_y  <- sum(y[rho_current])
        sum_y2 <- sum(y[rho_current]^2)
        njj    <- nj_current[j]
        
        # Iperparametri posterior (per Y)
        kappa_loc <- hyperpar_P0_y$kappa0 + njj
        mu_loc    <- (hyperpar_P0_y$kappa0 * hyperpar_P0_y$mu0 + sum_y) / kappa_loc
        alpha_loc <- hyperpar_P0_y$nu0/2 + njj/2
        
        mean_y_j <- sum_y / njj
        ss_y_j <- sum_y2 - (sum_y^2) / njj
        beta_loc  <- 0.5 * (hyperpar_P0_y$nu0 * hyperpar_P0_y$sig_square0) +
          0.5 * ss_y_j +
          0.5 * (hyperpar_P0_y$kappa0 * njj / kappa_loc) * (mean_y_j - hyperpar_P0_y$mu0)^2
        
        # Campiona tau_j ~ Gamma
        tauj_current[j] <- rgamma(1, shape = alpha_loc, rate = beta_loc)
        
        # Campiona mu_j ~ Normale
        muj_current[j] <- rnorm(1, 
                                mean = mu_loc, 
                                sd = sqrt(1 / (kappa_loc * tauj_current[j])))
      }
    }
    
    # --- 4c. SALVATAGGIO DEI RISULTATI (Thinning) ---
    
    if(g %in% saved_iters) {
      save_idx <- save_idx + 1
      
      k_out[save_idx] <- k_current
      clu_out[save_idx, ] <- clu_current
      nj_out[[save_idx]] <- nj_current
      # Salva solo i k_current parametri attivi
      tauj_out[[save_idx]] <- tauj_current[1:k_current]
      muj_out[[save_idx]] <- muj_current[1:k_current]
    }
    
    # Controllo di sanità occasionale
    if(g %% 500 == 0) {
      sanitize_state()
    }
    
    
    pb$tick()
    
  } # --- Fine loop MCMC (g) ---
  
  # --- 5. RITORNO DEI RISULTATI ---
  
  return(list(
    k_out = k_out,
    clu_out = clu_out,
    nj_out = nj_out,
    tauj_out = tauj_out,
    muj_out = muj_out,
    saved_iters = saved_iters
  ))
}