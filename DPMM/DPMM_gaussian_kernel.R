# -----------------------------------------------------------------
# Modello Gerarchico Bayesiano - Dirichlet Process Mixture Model (DPMM)
# -----------------------------------------------------------------
# Implementazione del Collapsed Gibbs Sampler (Neal, 2000 - Algorithm 8)
# con kernel Gaussiano e prior coniugata Normal-Gamma.
# -----------------------------------------------------------------

#' DPMM con Kernel Gaussiano e Prior Normal-Gamma
#'
#' @param y Vettore delle osservazioni.
#' @param n_iter Numero totale di iterazioni MCMC.
#' @param burn_in Numero di iterazioni iniziali da scartare.
#' @param thin Fattore di diradamento (salva ogni 'thin' iterazioni).
#' @param hyperpar_P0 Lista: mu0, kappa0, nu0, sig_square0 (parametri Base Measure P0).
#' @param hyperpar_Dir Lista: alpha (parametro di concentrazione DP).
#' @param init Lista: k, clu, nj, tauj, muj (stato iniziale).

DPMM_kernel_gaussiano <- function(y, n_iter, burn_in, thin,
                                  hyperpar_P0, hyperpar_Dir, init) {
  
  # --- 1. SETUP E VALIDAZIONE ---
  n <- length(y)
  if(n == 0) stop("Il vettore y non può essere vuoto.")
  
  # Calcolo delle iterazioni da salvare effettivamente
  saved_iters <- seq(from = burn_in + 1, to = n_iter, by = thin)
  G <- length(saved_iters)
  if(G == 0) stop("Parametri burn_in/thin/n_iter non validi: nessun campione salvato.")
  
  # --- 2. PRE-ALLOCAZIONE OUTPUT ---
  # Utilizziamo liste per parametri a dimensione variabile (nj, tauj, muj)
  k_out    <- rep(NA_integer_, G)
  clu_out  <- matrix(NA_integer_, nrow = G, ncol = n)
  nj_out   <- vector("list", G)
  tauj_out <- vector("list", G)
  muj_out  <- vector("list", G)
  
  # --- 3. INIZIALIZZAZIONE STATO CORRENTE ---
  k_current    <- init$k
  clu_current  <- as.integer(init$clu)
  nj_current   <- as.integer(init$nj)
  tauj_current <- as.numeric(init$tauj)
  muj_current  <- as.numeric(init$muj)
  
  # rho: struttura di supporto che indicizza quali 'i' appartengono a ogni cluster
  rho <- vector("list", k_current)
  for(j in seq_len(k_current)) rho[[j]] <- which(clu_current == j)
  
  # Funzione interna di validazione dello stato (debugging)
  sanitize_state <- function() {
    if(length(nj_current) != k_current) stop("Incoerenza lunghezza nj_current.")
    if(any(sapply(rho, length) != nj_current)) stop("Incoerenza tra rho e nj_current.")
  }
  
  save_idx <- 0 # Puntatore per il salvataggio degli output
  
  # --- 4. LOOP PRINCIPALE MCMC ---
  for(g in 1:n_iter) {
    
    # 4a. AGGIORNAMENTO DELLE ALLOCAZIONI (Sampling c_i)
    # Applichiamo il Polya Urn Scheme (Chinese Restaurant Process)
    for(i in 1:n) {
      
      # 1. Rimozione temporanea dell'osservazione i dallo stato
      clu_i <- clu_current[i]
      idx_in_rho <- which(rho[[clu_i]] == i)
      rho[[clu_i]] <- rho[[clu_i]][-idx_in_rho]
      nj_current[clu_i] <- nj_current[clu_i] - 1
      
      # 2. Se il cluster rimane vuoto, lo eliminiamo (Housekeeping)
      if(nj_current[clu_i] == 0) {
        nj_current <- nj_current[-clu_i]
        rho <- rho[-clu_i]
        tauj_current <- tauj_current[-clu_i]
        muj_current <- muj_current[-clu_i]
        k_current <- k_current - 1
        # Shift degli indici di allocazione per mantenere la continuità
        clu_current[clu_current > clu_i] <- clu_current[clu_current > clu_i] - 1
      }
      
      # 3. Calcolo dei pesi per il ricampionamento di c_i
      # w = [P(c_i = 1), ..., P(c_i = k), P(c_i = nuovo cluster)]
      w <- numeric(k_current + 1)
      
      if(k_current > 0) {
        for(j in 1:k_current) {
          # Statistiche sufficienti del cluster j (senza i)
          y_j <- y[rho[[j]]]
          njj <- nj_current[j]
          sum_y <- sum(y_j)
          sum_y2 <- sum(y_j^2)
          
          # Aggiornamento parametri posterior (Normal-Gamma)
          kn <- hyperpar_P0$kappa0 + njj
          an <- hyperpar_P0$nu0/2 + njj/2
          mun <- (hyperpar_P0$kappa0 * hyperpar_P0$mu0 + sum_y) / kn
          
          # Calcolo del parametro di scala della Gamma (bn)
          # Somma dei quadrati a priori + SS dati + termine di discrepanza medie
          bn <- 0.5 * (hyperpar_P0$nu0 * hyperpar_P0$sig_square0 + 
                         (sum_y2 - (sum_y^2)/njj) + 
                         (hyperpar_P0$kappa0 * njj / kn) * ((sum_y/njj) - hyperpar_P0$mu0)^2)
          
          # Densità predittiva t-Student: p(y_i | y_{-i}, c_i = j)
          # sig_n è il fattore di scala della t
          sig_n <- sqrt( bn * (kn + 1) / (an * kn) )
          df_n  <- 2 * an
          
          # Peso = Likelihood Marginale * Prior CRP (n_j)
          w[j] <- (dt((y[i] - mun) / sig_n, df = df_n) * (1 / sig_n)) * njj
        }
      }
      
      # 4. Peso per la creazione di un NUOVO cluster (basato su P0 e alpha)
      kappa0 <- hyperpar_P0$kappa0
      mu0    <- hyperpar_P0$mu0
      alpha0 <- hyperpar_P0$nu0 / 2
      beta0  <- 0.5 * (hyperpar_P0$nu0 * hyperpar_P0$sig_square0)
      sig0   <- sqrt( beta0 * (kappa0 + 1) / (alpha0 * kappa0) )
      df0    <- 2 * alpha0
      
      w[k_current + 1] <- (dt((y[i] - mu0) / sig0, df = df0) * (1 / sig0)) * hyperpar_Dir$alpha
      
      # 5. Campionamento della nuova allocazione
      if(all(w <= 0)) w <- rep(1e-10, length(w)) # Protezione da underflow numerico
      new_cluster <- sample.int(length(w), size = 1L, prob = w)
      
      # 6. Aggiornamento dello stato con la nuova allocazione
      clu_current[i] <- new_cluster
      if(new_cluster <= k_current) {
        nj_current[new_cluster] <- nj_current[new_cluster] + 1
        rho[[new_cluster]] <- c(rho[[new_cluster]], i)
      } else {
        # Creazione effettiva del nuovo cluster
        k_current <- k_current + 1
        nj_current <- c(nj_current, 1L)
        rho[[k_current]] <- i
        # Inizializzazione parametri per il nuovo cluster
        tauj_current <- c(tauj_current, 1.0)
        muj_current <- c(muj_current, y[i])
      }
    } # Fine campionamento c_i
    
    # 4b. CAMPIONAMENTO PARAMETRI CLUSTER (Sampling mu_j, tau_j)
    # Dato il clustering c, campioniamo i parametri dalla posterior Normal-Gamma
    for(j in 1:k_current) {
      y_j <- y[rho[[j]]]
      njj <- nj_current[j]
      sum_y <- sum(y_j)
      sum_y2 <- sum(y_j^2)
      
      kn  <- hyperpar_P0$kappa0 + njj
      mun <- (hyperpar_P0$kappa0 * hyperpar_P0$mu0 + sum_y) / kn
      an  <- hyperpar_P0$nu0/2 + njj/2
      bn  <- 0.5 * (hyperpar_P0$nu0 * hyperpar_P0$sig_square0 + 
                      (sum_y2 - (sum_y^2)/njj) + 
                      (hyperpar_P0$kappa0 * njj / kn) * ((sum_y/njj) - hyperpar_P0$mu0)^2)
      
      # Campionamento precisione (Gamma) e poi media (Normale condizionata)
      tau_j <- rgamma(1, shape = an, rate = bn)
      mu_j  <- rnorm(1, mean = mun, sd = sqrt(1 / (kn * tau_j)))
      
      tauj_current[j] <- tau_j
      muj_current[j]  <- mu_j
    }
    
    # 4c. SALVATAGGIO (Thinning)
    if(g %in% saved_iters) {
      save_idx <- save_idx + 1
      k_out[save_idx] <- k_current
      clu_out[save_idx, ] <- clu_current
      nj_out[[save_idx]]  <- nj_current
      tauj_out[[save_idx]] <- tauj_current
      muj_out[[save_idx]]  <- muj_current
    }
    
    # Feedback ogni 1000 iterazioni
    if(g %% 1000 == 0) cat("Iterazione:", g, "/", n_iter, "\n")
  }
  
  # --- 5. RETURN ---
  return(list(
    k_out = k_out, clu_out = clu_out, nj_out = nj_out,
    tauj_out = tauj_out, muj_out = muj_out, saved_iters = saved_iters
  ))
}