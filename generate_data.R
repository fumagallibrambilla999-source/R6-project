# --- generate_data_complex.R ---
set.seed(48)

n <- 1000

# Definiamo più cluster con medie e deviazioni diverse
means <- c(-3, -1, 0.5, 2, 3.5)       # 5 cluster
sds   <- c(0.3, 0.6, 0.4, 0.5, 0.3)  # deviazioni diverse
probs <- c(0.15, 0.25, 0.3, 0.2, 0.1) # probabilità di ciascun cluster

y <- numeric(n)
for(i in 1:n){
  u <- runif(1)
  cluster <- which(u <= cumsum(probs))[1]
  y[i] <- rnorm(1, mean = means[cluster], sd = sds[cluster])
}

# Aggiungiamo qualche outlier
n_out <- 10
outliers <- rnorm(n_out, mean = 6, sd = 1)
y <- c(y, outliers)

# Salviamo i dati
save(y, file = "data_generated.RData")

# Quick plot per visualizzare la distribuzione
hist(y, breaks = 50, col = "lightblue", main = "Synthetic univariate data (complex)", xlab = "y")
