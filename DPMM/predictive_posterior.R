# -----------------------------------------------------------------
# FUNZIONE PER IL CALCOLO DELLA DENSITÀ PREDITTIVA POSTERIORE
# -----------------------------------------------------------------
# Implementa la formula di integrazione Monte Carlo:
# p(y_nuovo | y_dati) = (1/G) * sum_{g=1 a G} [ p(y_nuovo | stato_g) ]
#
# Dove p(y_nuovo | stato_g) è:
# sum_{j=1 a k_g} [ n_j / (n + alpha) * N(y | mu_j, tau_j) ] + [ alpha / (n + alpha) * p(y | P0) ]
# -----------------------------------------------------------------

#' Calcola la Densità Predittiva Posteriore (Approccio Corretto)
#'
#' @param risultati_mcmc La lista di output della funzione `DPMM_kernel_gaussiano`.
#' @param hyperpar_P0 Lista degli iperparametri originali (mu0, kappa0, nu0, sig_square0).
#' @param hyperpar_Dir Lista degli iperparametri del DP (alpha).
#' @param n_obs Numero totale di osservazioni nei dati originali (n).
#' @param grid_points Vettore di punti 'y' (la griglia 'xs') su cui valutare la densità.
#'
#' @return Un vettore numerico (lungo quanto grid_points) contenente i valori di densità.

calcola_densita_predittiva <- function(risultati_mcmc, 
                                       hyperpar_P0, 
                                       hyperpar_Dir, 
                                       n_obs, 
                                       grid_points) {
  
  # --- 1. Estrai parametri MCMC e Iperparametri ---
  G <- length(risultati_mcmc$k_out)
  if (G == 0) stop("Nessuna iterazione salvata in risultati_mcmc.")
  
  alpha_val <- hyperpar_Dir$alpha
  
  # Denominatore comune (costante)
  denom <- n_obs + alpha_val
  
  # --- 2. Calcola il "Termine 2" (Nuovo Cluster) ---
  # p(y | P0) = t-Student(y | 2*alpha0, mu0, sig0^2)
  # Questo termine è costante per tutte le G iterazioni, 
  # quindi lo calcoliamo una sola volta.
  
  # Calcola i parametri della t-Student a priori
  mu0    <- hyperpar_P0$mu0
  kappa0 <- hyperpar_P0$kappa0
  alpha0 <- hyperpar_P0$nu0 / 2
  beta0  <- 0.5 * hyperpar_P0$nu0 * hyperpar_P0$sig_square0
  df0    <- 2 * alpha0
  
  # Calcola la scala (sd) della t-Student
  sig0_sq <- (beta0 * (kappa0 + 1)) / (alpha0 * kappa0)
  sig0    <- sqrt(sig0_sq)
  
  # Calcola la densità t-Student sulla griglia
  dens_t_prior <- dt((grid_points - mu0) / sig0, df = df0) * (1 / sig0)
  
  # Calcola il termine completo (pesato) per il nuovo cluster
  # Questo è il contributo medio della "parte 2" della formula
  # Nota: (1/G) * sum_{g=1 a G} [alpha / (n+alpha) * dens_t] = alpha / (n+alpha) * dens_t
  dens_term2_final <- (alpha_val / denom) * dens_t_prior
  
  
  # --- 3. Calcola il "Termine 1" (Cluster Esistenti) ---
  # (1/G) * sum_{g=1 a G} [ sum_{j=1 a k_g} (n_j / denom) * N(y | mu_j, tau_j) ]
  
  # Accumulatore per la somma esterna (su G)
  dens_term1_accumulator <- numeric(length(grid_points))
  
  for (g in 1:G) {
    # Estrai lo stato per l'iterazione g
    k_g   <- risultati_mcmc$k_out[g]
    
    # Assicurati che nj_out sia presente. Se l'hai rimosso,
    # dovresti ricalcolarlo prima di chiamare questa funzione.
    if (is.null(risultati_mcmc$nj_out)) {
      stop("risultati_mcmc$nj_out non trovato. È necessario per il calcolo dei pesi.")
    }
    
    nj_g  <- risultati_mcmc$nj_out[[g]]
    muj_g <- risultati_mcmc$muj_out[[g]]
    tauj_g<- risultati_mcmc$tauj_out[[g]]
    
    # Accumulatore per la somma interna (su j = 1...k_g)
    dens_g_inner_sum <- numeric(length(grid_points))
    
    if (k_g > 0) {
      for (j in 1:k_g) {
        # Parametri del cluster j all'iterazione g
        n_j     <- nj_g[j]
        mu_j    <- muj_g[j]
        sigma_j <- 1 / sqrt(tauj_g[j])
        
        # Peso per questo componente j
        weight_j <- n_j / denom
        
        # Calcola la densità Normale per questo componente
        dens_j_normal <- dnorm(grid_points, mean = mu_j, sd = sigma_j)
        
        # Aggiungi alla somma interna
        dens_g_inner_sum <- dens_g_inner_sum + (weight_j * dens_j_normal)
      }
    }
    
    # Aggiungi la somma interna (per l'iterazione g) all'accumulatore principale
    dens_term1_accumulator <- dens_term1_accumulator + dens_g_inner_sum
  }
  
  # Completa la media per il Termine 1
  dens_term1_final <- dens_term1_accumulator / G
  
  # --- 4. Somma i due termini ---
  densita_totale <- dens_term1_final + dens_term2_final
  
  return(densita_totale)
}


# -----------------------------------------------------------------
# ESEMPIO DI UTILIZZO
# -----------------------------------------------------------------

# Supponiamo di avere già eseguito il tuo MCMC e di avere:
#
# 1. y_dati: I tuoi dati originali
# 2. risultati_mcmc: L'output della tua funzione DPMM_kernel_gaussiano
# 3. iper_P0: La lista di iperparametri P0 usata per l'MCMC
# 4. iper_Dir: La lista di iperparametri Dir usata per l'MCMC

# Per testare, creiamo dei risultati fittizi:
# (NON ESEGUIRE QUESTA PARTE SE HAI GIÀ I TUOI RISULTATI VERI)
if (FALSE) {
  y_dati <- c(rnorm(50, -3), rnorm(50, 4))
  iper_P0 <- list(mu0 = 0, kappa0 = 1, nu0 = 4, sig_square0 = 1)
  iper_Dir<- list(alpha = 1)
  
  # Creiamo un finto 'risultati_mcmc' (con G=2 iterazioni)
  risultati_mcmc <- list(
    k_out = c(2, 2),
    nj_out = list(c(50, 50), c(49, 51)),
    muj_out = list(c(-3.1, 4.1), c(-2.9, 3.9)),
    tauj_out = list(c(1.1, 0.9), c(1.0, 1.1))
  )
}


# --- ESECUZIONE VERA E PROPRIA ---

# 1. Definisci la griglia di punti su cui plottare la densità
xs <- seq(min(y_dati) - 2*sd(y_dati), 
          max(y_dati) + 2*sd(y_dati), 
          length.out = 500)

# 2. Calcola la densità predittiva chiamando la nuova funzione
dens_predittiva <- calcola_densita_predittiva(
  risultati_mcmc = risultati_mcmc,
  hyperpar_P0 = iper_P0,
  hyperpar_Dir = iper_Dir,
  n_obs = length(y_dati),
  grid_points = xs
)

# 3. Crea il grafico
hist(y_dati, probability = TRUE, breaks = 40, col = "gray80",
     main = "Densità Predittiva Posteriore (Corretta)",
     xlab = "y")
lines(xs, dens_predittiva, lwd = 3, col = "blue")
legend("topright", "Densità Predittiva", col="blue", lwd=3, bty="n")
