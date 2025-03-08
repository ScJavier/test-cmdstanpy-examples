data {
    int<lower=0> N;                   // Número de observaciones
    array[N] vector[2] y;             // Datos
    real<lower=1e-3> lambda1;         // Paramétro de la normal 1
    real<lower=1e-3> lambda2;         // Paramétro de la normal 2
    real<lower=1e-3> alpha1;          // Parámetro de la gamma 1
    real<lower=1e-3> beta1;           // Parámetro de la gamma 1
    real<lower=1e-3> alpha2;          // Parámetro de la gamma 2
    real<lower=1e-3> beta2;           // Parámetro de la gamma 2
}
parameters {
    real<lower=1e-3> mu1;                 // Parámetro latente
    real<lower=1e-3> mu2;                 // Parámetro latente
    real<lower=1e-3> tau1;                // Parámetro latente
    real<lower=1e-3> tau2;                // Parámetro latente
    real<lower=-0.999,upper=0.999> rho;   // Parámetro latente
}

transformed parameters {
    vector[2] mu = [mu1, mu2]';          // Vector de medias
    matrix[2,2] Sigma;                   // Matriz de covarianza
    Sigma[1,1] = tau1^2;
    Sigma[1,2] = rho * tau1 * tau2;
    Sigma[2,1] = rho * tau1 * tau2;
    Sigma[2,2] = tau2^2;
}

model {
    mu1 ~ normal(0, lambda1);          // Prior
    tau1 ~ gamma(alpha1, beta1);       // Prior
    mu2 ~ normal(0, lambda2);          // Prior
    tau2 ~ gamma(alpha2, beta2);       // Prior
    rho ~ uniform(-1, 1);

    y ~ multi_normal(mu, Sigma);       // Likelihood
}

generated quantities {
    vector[2] y_new;                   // Muestras simuladas
    y_new = multi_normal_rng(mu, Sigma); 
}