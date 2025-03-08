data {
    int<lower=0> N;                // Número de observaciones
    array[N] real y;       // Datos
    real<lower=1e-3> lambda;       // Parámetro de la normal
    real<lower=1e-3> alpha;        // Parámetro de la gamma
    real<lower=1e-3> beta;         // Parámetro de la gamma
}
parameters {
    real mu;           // Parámetro latente
    real<lower=1e-3> tau;          // Parámetro latente
}
model {
    mu ~ normal(0, lambda);          // Prior
    tau ~ gamma(alpha, beta);      // Prior
    y ~ normal(mu, tau);           // Likelihood
}