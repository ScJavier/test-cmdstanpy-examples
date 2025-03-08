data {
    int<lower=0> N;                   // Número de observaciones
    array[N] vector[2] y;             // Datos
    matrix[2,2] Lambda;               // Parametro de la normal para mu
    real<lower=1e-3> alpha;           // Parametro de la matriz de covarianzas
    real<lower=1e-3> beta;            // Parametro de la matriz de covarianzas
}

parameters {
    vector[2] mu;                     // Vector de medias
    real<lower=1e-3> nu;              // Grados de libertad
    cholesky_factor_corr[2] L_Omega;  // Matriz de correlaciones de Cholesky
    vector<lower=1e-3>[2] tau;           // Desviaciones estándar
}

transformed parameters {
    matrix[2,2] Sigma;                          // Matriz de covarianza
    matrix[2,2] L_Sigma;
    L_Sigma = diag_pre_multiply(tau, L_Omega);  // Construcción de la covarianza
    Sigma = L_Sigma * L_Sigma';

    

}

model {
    vector[2] mu_0 = rep_vector(0, 2);  // Vector de medias en 0

    // Prior sobre mu (multivariada)
    mu ~ multi_normal(mu_0, Lambda);

     // Priors sobre la matriz de covarianza (factorizada con LKJ)
    nu ~ gamma(alpha, beta);             // Priors en grados de libertad
    tau ~ gamma(alpha, beta);            // Priors en desviaciones estándar
    L_Omega ~ lkj_corr_cholesky(nu);     // Prior en la correlación
    y ~ multi_normal(mu, Sigma);         // Likelihood
}

generated quantities {
    vector[2] y_new;  // Muestras simuladas
    real tau1;  // Raíz cuadrada de sigma_11
    real tau2;  // Raíz cuadrada de sigma_22
    real rho;  // Cálculo de la correlacion

    y_new = multi_normal_rng(mu, Sigma);

    // Calcular sqrt de los elementos diagonales de Sigma
    tau1 = sqrt(Sigma[1,1]);
    tau2 = sqrt(Sigma[2,2]);

    // Calcular rho basado en la matriz de covarianza
    rho = Sigma[1,2] / (tau1 * tau2);
}