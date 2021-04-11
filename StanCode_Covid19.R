###################Naive Model########################
Naive_poisson = "data{
    int<lower=0> N;
    int<lower=0> z[N];//Observed counts
    vector<lower=0>[N] E;//exposure
    int<lower=1> K;
    matrix[N,K] X;//design matrix for counts
}
transformed data{
    vector[N] logE=log(E);
}
parameters{
//coefficients for Poisson regression
  real gamma0;
  vector[K] gammas;
 // vector[N] theta; //random effects
  real<lower=0> sigma;
}
model{
//likelihood
  z ~ poisson_log(gamma0+X*gammas + logE);
 // theta ~ normal(0,sigma);
  gamma0 ~ normal(-5,2);
  gammas ~ normal(0,5);
 // sigma ~ normal(0,2);
}
generated quantities{
 
}
"
###############
compiled_naive = rstan::stan_model(model_code = Naive_poisson)

##########M2: Under_Non_Spatial#########
Under_Non_Spatial = "data{
    int<lower=0> N;
    int<lower=0> z[N];//Observed counts
    vector<lower=0>[N] E;//exposure
    int<lower=1> K;
    matrix[N,K] X;//design matrix for counts
    //
    int<lower=1> J;
    matrix[N,J] W;//design matrix for reporting
}
transformed data{
    vector[N] logE=log(E);
}
parameters {
//coefficients for Poisson
  real beta0;
  vector[K] gammas;
//coefficients for Logistic
  real gamma0;
  vector[J] betas;
//real logit_rho;
//vector[N] phi;// spatial random effects
//vector[N] theta; //unstructured random effects
//real<lower=0> sigma;
}
transformed parameters{
  real<lower=0,upper=1> p0;//reporting rate when W take average
  vector<lower=0,upper=1>[N] p;
  vector<lower=0>[N] lambda;
  vector<lower=0>[N] mu;
//real<lower=0,upper=1> rho = inv_logit(logit_rho);
//vector[N] convolved_re;
//convolved_re = sqrt(1-rho)*theta +sqrt(rho/scaling_factor)*phi;
  lambda = exp(logE + gamma0 + X*gammas);
  p0 = inv_logit(beta0);
  p = inv_logit(beta0 + W*betas);
  mu = lambda .* p;
}
model {
//likelihood
  z ~ poisson(mu);
  target += beta0-2*log(1+exp(beta0));//Jacobian adjustment
//priors
  p0 ~ beta(7,55);//will induce prior on beta0
  //p0 ~ beta(21,181);
  gammas ~ normal(-5,2);
  gamma0 ~ normal(0,5);
  betas ~ normal(0,5);
//logit_rho ~ normal(0,1);
//sigma ~ normal(0,1);
//theta ~ normal(0,1);
//phi ~ icar_normal(N, node1, node2);
}
generated quantities{
 
}
"
###############

compiled_Under_Non_Spatial = rstan::stan_model(model_code = Under_Non_Spatial)

##########M3: Spa_Non_Under#########
Spa_Non_Under = "
functions{
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[]node2) {
  return -0.5 * dot_self(phi[node1] - phi[node2])+ normal_lpdf(sum(phi) | 0,0.001 * N);}
  }
data{
    int<lower=0> N;
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];
    int<lower=1, upper=N> node2[N_edges];
    real<lower=0> scaling_factor;
    int<lower=0> z[N];//Observed counts
    vector<lower=0>[N] E;//exposure
    int<lower=1> K;
    matrix[N,K] X;//design matrix for counts
}
transformed data{
    vector[N] logE=log(E);
}
parameters{
//coefficients for Poisson regression
  real gamma0;
  vector[K] gammas;
  vector[N] theta; //random effects
  real logit_rho;
  vector[N] phi;// spatial random effects
  real<lower=0> sigma;
}
transformed parameters{
  //vector<lower=0,upper=1>[N] P;
  vector<lower=0>[N] lambda;
  //vector<lower=0>[N] mu;
  
  real<lower=0,upper=1> rho = inv_logit(logit_rho);
  vector[N] convolved_re;
  convolved_re = sqrt(1-rho)*theta +sqrt(rho/scaling_factor)*phi;
  lambda = exp(logE + gamma0+X*gammas + convolved_re*sigma);
  //P = inv_logit(W*betas);
  //mu = lambda .* P;
}
model{
//likelihood
  z ~ poisson(lambda);
  
// priors  
// betas ~ normal(0,5);
  gamma0 ~ normal(0,5);
  gammas ~ normal(0,5);
  
  logit_rho ~ normal(0,1);
  sigma ~ normal(0,1);
  theta ~ normal(0,1);
  phi ~ icar_normal(N, node1, node2);
}
generated quantities{
 
}
"

compiled_Spa_Non_Under = rstan::stan_model(model_code = Spa_Non_Under)


##################Full_BYM2 model#############
covid19_bym2_under <- "functions{
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[]node2) {
  return -0.5 * dot_self(phi[node1] - phi[node2])+ normal_lpdf(sum(phi) | 0,0.001 * N);}
  }
  data {
    int<lower=0> N;
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];
    int<lower=1, upper=N> node2[N_edges];
    int<lower=0> z[N];//Observed counts
    vector<lower=0>[N] E;//exposure
    int<lower=1> K;
    matrix[N,K] x;//design matrix for counts
    int<lower=1> J;
    matrix[N,J] w;//design matrix for reporting rate
    real<lower=0> scaling_factor;
}
transformed data{
    vector[N] logE=log(E);
}
parameters {
//coefficients for Poisson
  real beta0;
  vector[K] gammas;
//coefficients for Logistic
  real gamma0;
  vector[J] betas;
  real logit_rho;
  vector[N] phi;// spatial random effects
  vector[N] theta; //unstructured random effects
  real<lower=0> sigma;
}
transformed parameters{
  real<lower=0,upper=1> p0;//reporting rate when W take average
  vector<lower=0,upper=1>[N] p;
  vector<lower=0>[N] lambda;
  vector<lower=0>[N] mu;
  real<lower=0,upper=1> rho = inv_logit(logit_rho);
  vector[N] convolved_re;
  convolved_re = sqrt(1-rho)*theta +sqrt(rho/scaling_factor)*phi;
  lambda = exp(logE + gamma0 + x*gammas + convolved_re*sigma);
  p0 = inv_logit(beta0);
  p = inv_logit(beta0 + w*betas);
  mu = lambda .* p;
}
model {
//likelihood
  z ~ poisson(mu);
  target += beta0-2*log(1+exp(beta0));//Jacobian adjustment
//priors
  p0 ~ beta(7,55);//will induce prior on gamma0
  //p0 ~ beta(21,181);
  gammas ~ normal(0,5);
  gamma0 ~ normal(0,5);
  betas ~ normal(0,5);
  logit_rho ~ normal(0,1);
  sigma ~ normal(0,1);
  theta ~ normal(0,1);
  phi ~ icar_normal(N, node1, node2);
}
generated quantities{
  vector[N] eta = logE + gamma0 + x*gammas + convolved_re*sigma;
  vector[N] lambda_rep = exp(eta);
  vector[N] p_rep = inv_logit(beta0 + w*betas);
  int y_rep[N];
  int z_rep[N];
  if(max(eta)>20){
  for ( n in 1:N) {
  y_rep[n] = -1;
  z_rep[n] = -1;
  }
  }else{
  for ( n in 1:N) {
  y_rep[n] = poisson_log_rng(eta[n]);
  z_rep[n] = poisson_rng(lambda_rep[n]*p_rep[n]);
  }
  }
}
"
compiled_full_BYM2 = rstan::stan_model(model_code=covid19_bym2_under)

#####################