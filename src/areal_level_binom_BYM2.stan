/* 
Some resources 
BYM(2) in Stan
http://www.stat.columbia.edu/~gelman/research/published/bym_article_SSTEproof.pdf
https://www.sciencedirect.com/science/article/pii/S1877584518301175
SAE, variance smoothing model - the some parts of the code is not compatible with newer stan versions hence updated
https://github.com/peteragao/VSALM/blob/main/inst/stan/spatial_joint_smooth_logit.stan 

*/
functions {
    real icar_normal_lpdf(vector theta, int N, array[] int node1, array[] int node2) {
        return -0.5 * dot_self(theta[node1] - theta[node2]) + 
        normal_lpdf(sum(theta) | 0, 0.001 * N); 
        //the second term added for soft sum-to-zero contraits
    }
}
data{
    int<lower=1> N; // number of admn2 areas
    int<lower=1> NS; // number of admin2 areas with valid estimates 
    //array[NS] int<lower=1,upper=N> adm2_index; //index
    array[NS] int<lower=1, upper=N> adm2_index;
    array[NS] real<lower=0, upper =1> p_hat; // direct estimate of prevalence, could be omitted
    array[NS] int<lower=0> n_star; //effective sample size
    array[NS] int<lower=0> y_star;

    int<lower=1> N_edges ; 
    array[N_edges] int<lower=1, upper=N> node1 ;
    array[N_edges] int<lower=1, upper=N> node2 ;
    real<lower=0> scaling_factor; 
}
transformed data {
    real delta = 1e-9;
    //vector[NS] ystar = floor(to_vector(p_hat) .* to_vector(nstar));
}
parameters{
    // For mean model (BYM2)
    vector[N] u1;// structured effect
    vector[N] u2;// random effect
    real<lower=0> sigma_u;
    real<lower=0,upper=1> rho;
    real b0;
}
transformed parameters {
    vector[N] u = (sigma_u)*(sqrt(1-rho) * (u1) + sqrt(rho/scaling_factor) * (u2) );
    vector[N] p = inv_logit(b0+u);
}
model{
    // likelihood
    //target += lgamma(to_vector(nstar)) - lgamma(ystar) - lgamma(to_vector(nstar)-ystar) + log(to_vector(nstar)) - log(ystar) - log(to_vector(nstar)-ystar) + (ystar).*log(p[adm2_index])+(to_vector(nstar)-ystar).*log(1-p[adm2_index]);
    target += binomial_lpmf(y_star|n_star,p[adm2_index]) ;
    // mean model prior
    target += normal_lpdf(u1|0,1);   
    target += icar_normal_lpdf(u2|N, node1,node2); 
    target += normal_lpdf(sigma_u|0,1); //change this to penalised complexity prior 
    target += beta_lpdf(rho| 0.5, 0.5); 
    target +=normal_lpdf(b0|0,5);
}