# Input dati dalla funzione:
# 1) y --> dati da clusterizzare
# 2) n_iter --> numero di iterazini che verranno fatte dalla funzione
# 3) burn_in --> numero di iterazioni da scartare per far convergere la catena di MArkov
# 4) thin --> numero di iterazioni dopo cui si salva lo stato della MC nell'output
# 5) hyperpar_P0 --> lista di iperparametri della distribuzione P0 (= NormalGamma(mu0, kappa0, nu0, sig_square0))
# 6) hyperpar_Dir --> lista di iperparametri della distribuzione di Dirichlet (sostanzialmente contiene alpha)
# 7) init --> lista delle condizioni iniziali della MC:
#                                                       k = numero iniziale di clusters
#                                                       clu = vettore che indica a quale cluster appartiene ogni dato
#                                                       nj = vettore che indica la numerosità di ogni cluster
#                                                       tauj = vettore dei valori iniziali di tau per ogni cluster
#                                                       muj = vettore dei valori iniziali di mu per ogni cluster



#Notazione: Peter Hoff

# -----------------------------------------------------------------
# Modello Gerarchico Bayesiano - Dirichlet Process Mixture Model (DPMM)
# -----------------------------------------------------------------

# 1. LIVELLO DATI (LIKELIHOOD)
# Le osservazioni Y_i (per i = 1, ..., n) sono condizionatamente indipendenti
# dato un set di parametri theta_i.
#
# Y_i | theta_i  ~  f(y_i | theta_i)
#
#f è la densità di una distribuzione Normale:
# theta_i = (mu_i, sigma_i^2)
# Y_i | mu_i, sigma_i^2  ~  N(mu_i, sigma_i^2)

#tau_i = 1/sigma_i^2

# 2. LIVELLO PROCESSO
# I parametri theta_i sono estratti in modo indipendente e identico (iid)
# da una distribuzione (discreta) P.
#
# theta_1, ..., theta_n | P  ~ iid P


# 3. LIVELLO PRIOR (PROCESSO DI DIRICHLET)
# La distribuzione P è a sua volta una variabile aleatoria che segue
# un Processo di Dirichlet (DP).
#
# P ~ DP(alpha, P_0)
#
# - alpha: è il parametro di concentrazione (scalare).
#(hyperpar_Dir$alpha)

# - P_0: è la misura di base (base measure), che rappresenta la nostra
#        "credenza a priori" sulla forma di P.


# 4. MISURA DI BASE (PRIORE SUI PARAMETRI)
# La misura di base P_0 è la distribuzione a priori per i parametri
# (mu, sigma^2) di ciascun cluster.
#
# P_0 è la distribuzione di una Normal-InverseGamma(mu0, kappa0, nu0, sig_square0)
#
# (Questa è la prior coniugata standard per la media (mu) e la
# varianza (sigma^2) di una distribuzione Normale).




# DPMM_kernel_gaussiano: corrected and robust implementation
DPMM_kernel_gaussiano <- function(y, n_iter, burn_in, thin,
                                  hyperpar_P0, hyperpar_Dir, init) {
  # y: vettore di osservazioni
  # hyperpar_P0: lista(mu0, kappa0, nu0, sig_square0) -> iperparametri di P0
  # hyperpar_Dir: lista(alpha) -> parametro di concentrazione del DP
  # init: lista(k, clu, nj, tauj, muj) -> stato iniziale della catena
  
  # --- 1. SETUP INIZIALE E VALIDAZIONE INPUT ---
  
  # n = numero totale di osservazioni
  n <- length(y)
  if(n == 0) stop("y must be non-empty")
  if(burn_in >= n_iter) stop("burn_in must be < n_iter")
  if(thin <= 0) stop("thin must be positive")
  
  # Calcola gli indici delle iterazioni che verranno salvate
  # Inizia dopo il burn_in, e salva un'iterazione ogni 'thin'
  saved_iters <- seq(from = burn_in + 1, to = n_iter, by = thin)
  
  # G = numero totale di campioni che verranno salvati
  G <- length(saved_iters)
  if(G == 0) stop("No iterations will be saved: check burn_in/thin/n_iter")
  
  # --- 2. PRE-ALLOCAZIONE DEGLI OGGETTI DI OUTPUT ---
  
  # k_out: vettore per salvare il numero di cluster (k) ad ogni iterazione salvata
  k_out   <- rep(NA_integer_, G)
  
  # clu_out: matrice G x n. Ogni riga [g,] conterrà il vettore delle
  #          allocazioni ai cluster (clu_current) all'iterazione salvata 'g' (appunti: c_i^(g))
  clu_out <- matrix(NA_integer_, nrow = G, ncol = n)
  
  # nj_out: LISTA di lunghezza G. Ogni elemento [[g]] conterrà il vettore
  #         delle numerosità dei cluster (nj_current). 
  nj_out  <- vector("list", G)
  
  # tauj_out, muj_out: LISTE di lunghezza G, come nj_out, per salvare
  #                  i parametri (precisione e media) di ciascun cluster. (sigma^2 = 1/tau)
  tauj_out<- vector("list", G)
  muj_out <- vector("list", G)
  
  # --- 3. INIZIALIZZAZIONE DELLO STATO CORRENTE ---
  # Questi sono i valori che verranno aggiornati ad *ogni* iterazione g
  
  # k_current = numero di cluster *attualmente* attivi
  k_current   <- init$k
  
  # clu_current = vettore (lungo n) delle allocazioni *attuali*
  clu_current <- as.integer(init$clu)
  
  # nj_current = vettore (lungo k_current) delle numerosità *attuali*
  nj_current  <- as.integer(init$nj)
  
  # tauj_current, muj_current = vettori (lunghi k_current) dei parametri *attuali*
  tauj_current<- as.numeric(init$tauj)
  muj_current <- as.numeric(init$muj)
  
  # 'rho' è una struttura dati di supporto. È una lista (lunga k_current)
  # dove rho[[j]] contiene gli INDICI (da 1 a n) delle osservazioni
  # che appartengono al cluster j.
  rho <- vector("list", k_current)
  # Popoliamo 'rho' in base alle allocazioni iniziali (clu_current)
  for(j in seq_len(k_current)) rho[[j]] <- which(clu_current == j)
  
  # Funzione di utilità per controllare la coerenza dello stato
  sanitize_state <- function() {
    if(length(nj_current) != k_current) stop("nj_current length mismatch k_current")
    if(length(rho) != k_current) stop("rho length mismatch k_current")
    if(any(sapply(rho, length) != nj_current)) stop("rho and nj_current inconsistent")
    if(any(unlist(rho) < 1 | unlist(rho) > n)) stop("rho contains invalid indices")
  }
  # Controllo iniziale per assicurarsi che 'init' sia valido
  sanitize_state()
  
  # save_idx: contatore per gli output. Indica in quale riga/elemento
  # degli oggetti *_out stiamo salvando. Va da 1 a G.
  save_idx <- 0
  
  # --- 4. INIZIO DEL LOOP MCMC (GIBBS SAMPLER) ---
  for(g in 1:n_iter) {
    
    # --- 4a. CAMPIONAMENTO DELLE ALLOCAZIONI (c_i) ---
    # Questo è il cuore del "Collapsed Gibbs Sampler" (noto come Algoritmo 8 di Neal).
    # Iteriamo su ogni osservazione i=1,...,n e campioniamo la sua
    # allocazione c_i dalla distribuzione condizionata P(c_i | c_{-i}, y).
    
    for(i in 1:n) {
      
      # --- 4a.1 Rimuovi l'osservazione 'i' dal suo cluster ---
      # Dobbiamo campionare c_i condizionato a c_{-i} (tutti gli altri).
      # Per farlo, prima "togliamo" l'osservazione i dallo stato corrente.
      
      # clu_i = cluster a cui 'i' è *attualmente* assegnato
      clu_i <- clu_current[i]
      
      # Trova la posizione di 'i' nel vettore di indici rho[[clu_i]]
      idx_in_rho <- which(rho[[clu_i]] == i)
      if(length(idx_in_rho) != 1) stop("index_to_remove problem")
      
      # Rimuovi l'indice 'i' da rho
      rho[[clu_i]] <- rho[[clu_i]][-idx_in_rho]
      
      # Riduci il contatore per quel cluster
      nj_current[clu_i] <- nj_current[clu_i] - 1
      
      # --- 4a.2 Gestione cluster vuoti ---
      # Se rimuovendo 'i' il cluster clu_i è diventato vuoto,
      # dobbiamo rimuoverlo completamente e ri-etichettare.
      if(nj_current[clu_i] == 0) {
        
        # Rimuovi il cluster (e i suoi parametri) dai vettori di stato
        # NOTA: qui si assume che tauj_current e muj_current vengano
        # ricampionati da zero alla fine del loop, quindi non
        # è strettamente necessario rimuoverli. Si rimuovono
        # nj_current e rho che servono per il calcolo dei pesi.
        nj_current <- nj_current[-clu_i]
        rho <- rho[-clu_i]
        
        # Riduci il numero totale di cluster
        k_current <- k_current - 1
        
        # Ri-etichetta: tutti i cluster con indice > clu_i
        # vengono scalati di -1 (es. se rimuovo il cluster 3, il 4 diventa 3, il 5 diventa 4, ...)
        clu_current[clu_current > clu_i] <- clu_current[clu_current > clu_i] - 1
      }
      
      # --- 4a.3 Calcolo Pesi per l'Allocazione (Polya Urn Scheme) ---
      # Dobbiamo calcolare la probabilità (non normalizzata) 'w'
      # che l'osservazione 'i' appartenga a:
      # 1. Ciascuno dei k_current cluster esistenti (w[1...k_current])
      # 2. Un *nuovo* cluster (w[k_current + 1])
      
      # w = vettore dei pesi (probabilità non normalizzate)
      w <- numeric(k_current + 1)
      
      # Calcola i pesi per i k_current cluster *esistenti*
      if(k_current > 0) {
        for(j in 1:k_current) {
          
          # Statistiche sufficienti del cluster j (SENZA l'osservazione 'i')
          rho_current <- rho[[j]] # Indici delle osservazioni in cluster j
          sum_y  <- sum(y[rho_current])
          sum_y2 <- sum(y[rho_current]^2)
          njj    <- nj_current[j] # Numerosità del cluster j
          
          # --- Calcolo della DISTRIBUZIONE PREDITTIVA POSTERIORE ---
          # P(c_i = j | c_{-i}, y) \propto P(c_i = j | c_{-i}) * P(y_i | y_j)
          # P(c_i = j | c_{-i}) = njj (dal Chinese Restaurant Process)
          # P(y_i | y_j) = densità predittiva t-Student
          
          # Iperparametri "locali" (aggiornati) per la Normal-Gamma
          # Questo è l'aggiornamento Bayesiano standard
          kappa_loc <- hyperpar_P0$kappa0 + njj
          mu_loc    <- (hyperpar_P0$kappa0 * hyperpar_P0$mu0 + sum_y) / kappa_loc
          alpha_loc <- hyperpar_P0$nu0/2 + njj/2
          
          # beta_loc è il parametro di 'rate' (tasso) della Gamma
          beta_loc  <- 0.5 * (
            hyperpar_P0$nu0 * hyperpar_P0$sig_square0 + # Dalla prior
              (sum_y2 - (sum_y^2) / njj) + # Dalla somma degli scarti quadratici dei dati
              (hyperpar_P0$kappa0 * njj / kappa_loc) * ((sum_y / njj) - hyperpar_P0$mu0)^2 # Dalla differenza tra medie
          )
          
          # Parametri della t-Student predittiva
          # (derivata integrando via mu e tau dalla Normal-Gamma)
          sig_loc <- sqrt( beta_loc * (kappa_loc + 1) / (alpha_loc * kappa_loc) )
          df_loc  <- 2 * alpha_loc
          
          # Calcolo del peso w[j]:
          # w[j] = (densità t-Student di y[i]) * (prob. prior del cluster, njj)
          # dt(...) * (1 / sig_loc) calcola la densità di una t-Student
          # con media mu_loc, scala sig_loc e g.d.l. df_loc
          w[j] <- (dt((y[i] - mu_loc) / sig_loc, df = df_loc) * (1 / sig_loc)) * njj
        }
      }
      
      # --- Calcolo della DISTRIBUZIONE PREDITTIVA PRIORE ---
      # P(c_i = k+1 | c_{-i}, y) \propto P(c_i = k+1 | c_{-i}) * P(y_i | P0)
      # P(c_i = k+1 | c_{-i}) = alpha (dal Chinese Restaurant Process)
      # P(y_i | P0) = densità predittiva *a priori* (una t-Student basata solo su P0)
      
      # Calcola il peso per un *nuovo* cluster (j = k_current + 1)
      # È la stessa formula di prima, ma con njj = 0
      kappa0 <- hyperpar_P0$kappa0
      mu0    <- hyperpar_P0$mu0
      alpha0 <- hyperpar_P0$nu0 / 2
      beta0  <- 0.5 * (hyperpar_P0$nu0 * hyperpar_P0$sig_square0)
      sig0   <- sqrt( beta0 * (kappa0 + 1) / (alpha0 * kappa0) )
      df0    <- 2 * alpha0
      
      # Calcolo del peso w[k_current + 1]:
      # w[k+1] = (densità t-Student di y[i]) * (prob. prior del cluster, alpha)
      w[k_current + 1] <- (dt((y[i] - mu0) / sig0, df = df0) * (1 / sig0)) * hyperpar_Dir$alpha
      
      # Controllo di stabilità numerica (se tutti i pesi sono 0 per underflow)
      if(all(w == 0)) w <- rep(1e-12, length(w))
      
      # --- 4a.4 Campiona e assegna il nuovo cluster ---
      # Campiona un indice (da 1 a k_current+1) usando i pesi 'w'
      new_cluster <- sample.int(length(w), size = 1L, prob = w)
      
      # Assegna l'osservazione 'i' al cluster campionato
      clu_current[i] <- new_cluster
      
      # --- 4a.5 Aggiorna lo stato ---
      if(new_cluster <= k_current) {
        # Caso 1: 'i' è stato assegnato a un cluster *esistente*
        nj_current[new_cluster] <- nj_current[new_cluster] + 1
        rho[[new_cluster]] <- c(rho[[new_cluster]], i)
      } else {
        # Caso 2: 'i' è stato assegnato a un *nuovo* cluster
        # Aggiungi il nuovo cluster allo stato
        nj_current <- c(nj_current, 1L) # Aggiungi 1 al vettore delle numerosità
        k_current <- k_current + 1 # Incrementa il conteggio dei cluster
        rho[[k_current]] <- i # Aggiungi un nuovo elemento alla lista rho
        
        # Aggiungiamo placeholder per i parametri (verranno campionati sotto)
        tauj_current <- c(tauj_current, 1.0) # Valore a caso, non importa
        muj_current <- c(muj_current, 0.0)   # Valore a caso, non importa
      }
    } # --- Fine loop su i (osservazioni) ---
    
    # A questo punto, tutte le allocazioni c_i sono state aggiornate
    
    # --- 4b. CAMPIONAMENTO DEI PARAMETRI DEL CLUSTER (theta_j) ---
    # Ora campioniamo i parametri (mu_j, tau_j) per ogni cluster j=1...k_current
    # dalla loro distribuzione condizionata P(theta_j | c, y)
    # Questa è la distribuzione Normal-Gamma a posteriori
    
    if(k_current > 0) {
      for(j in 1:k_current) {
        
        # Statistiche sufficienti per il cluster j (ora complete)
        rho_current <- rho[[j]]
        sum_y  <- sum(y[rho_current])
        sum_y2 <- sum(y[rho_current]^2)
        njj    <- nj_current[j]
        
        # Calcola gli iperparametri della posterior Normal-Gamma
        # (sono identici a quelli calcolati in 4a.3)
        kappa_loc <- hyperpar_P0$kappa0 + njj
        mu_loc    <- (hyperpar_P0$kappa0 * hyperpar_P0$mu0 + sum_y) / kappa_loc
        alpha_loc <- hyperpar_P0$nu0/2 + njj/2
        beta_loc  <- 0.5 * (
          hyperpar_P0$nu0 * hyperpar_P0$sig_square0 +
            (sum_y2 - (sum_y^2) / njj) +
            (hyperpar_P0$kappa0 * njj / kappa_loc) * ((sum_y / njj) - hyperpar_P0$mu0)^2
        )
        
        # Campionamento dalla P(tau | y, c)
        # tau_j ~ Gamma(shape=alpha_loc, rate=beta_loc)
        tauj_current[j] <- rgamma(1, shape = alpha_loc, rate = beta_loc)
        
        # Campionamento dalla P(mu | tau, y, c)
        # mu_j ~ Normal(mean=mu_loc, variance = 1 / (kappa_loc * tau_j))
        muj_current[j] <- rnorm(1, 
                                mean = mu_loc, 
                                sd = sqrt(1 / (kappa_loc * tauj_current[j])))
      }
    }
    
    # --- 4c. SALVATAGGIO DEI RISULTATI (Thinning) ---
    
    # Controlla se l'iterazione corrente 'g' deve essere salvata
    if(g %in% saved_iters) {
      
      # Incrementa l'indice di salvataggio
      save_idx <- save_idx + 1
      
      # Salva lo stato corrente negli oggetti di output
      k_out[save_idx] <- k_current
      clu_out[save_idx, ] <- clu_current
      nj_out[[save_idx]] <- nj_current
      # Salva solo i parametri dei k_current cluster attivi
      tauj_out[[save_idx]] <- tauj_current[1:k_current]
      muj_out[[save_idx]] <- muj_current[1:k_current]
    }
    
    # Controllo di sanità occasionale
    if(g %% 500 == 0) {
      sanitize_state()
    }
  } # --- Fine loop MCMC (g) ---
  
  # --- 5. RITORNO DEI RISULTATI ---
  
  return(list(
    k_out = k_out,       # Vettore (lungo G) del numero di cluster
    clu_out = clu_out,   # Matrice (G x n) delle allocazioni
    nj_out = nj_out,     # Lista (lunga G) dei vettori di numerosità
    tauj_out = tauj_out, # Lista (lunga G) dei vettori di precisioni
    muj_out = muj_out,   # Lista (lunga G) dei vettori di medie
    saved_iters = saved_iters # Vettore degli indici delle iterazioni salvate
  ))
}

