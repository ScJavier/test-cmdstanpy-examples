data {
    int<lower=0> N;                  // Número de observaciones
    array[N] int<lower=0> y;         // Datos de conteo (Corrección aquí)
    real<lower=1e-3> alpha;              // Parámetro de la gamma
    real<lower=1e-3> beta;               // Parámetro de la gamma
}
parameters {
    real<lower=1e-3> lambda;             // Parámetro latente
}
model {
    lambda ~ gamma(alpha, beta);       // Prior
    y ~ poisson(lambda);               // Likelihood
}
