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
    array[NS] int<lower=1,upper=N> adm2_index; //index
    array[NS] real<lower=0, upper =1> p_hat; // direct estimate of prevalence
    array[NS] real<lower=0> v_hat; // direct estimate of variance
    array[NS] real<lower=0> d; // degree of freedom in each area, at admin2 level ; 
    array[NS] int<lower=0> k; //number of respondants in each area, at admin2 level 

    int<lower=1> N_edges ; 
    array[N_edges] int<lower=1, upper=N> node1 ;
    array[N_edges] int<lower=1, upper=N> node2 ;
    real<lower=0> scaling_factor; 
}
transformed data {
    real delta = 1e-9;
    vector[NS] d_vhat = (to_vector(d).* to_vector(v_hat));
}
parameters{
    // For mean model (BYM2)
    vector[N] u1;// structured effect
    vector[N] u2;// random effect
    real<lower=0> sigma_u;
    real<lower=0,upper=1> rho;
    // For variance model
    vector[NS] tau;
    real<lower=0> sigma_tau ;
    real gamma0;
    real gamma1;
    real gamma2;
}
transformed parameters {
    vector[N] u = (sigma_u)*(sqrt(rho/scaling_factor) * (u1) + sqrt(1-rho) * (u2));
    vector[N] p = inv_logit(u);
    vector[NS] v = exp(gamma0 + gamma1*log(p[adm2_index].*(1-p[adm2_index])) + gamma2*log(to_vector(k)) + square(sigma_tau)*tau);
    vector[NS] scaled_vhat= d_vhat./v;
}
model{
    // likelihood
    target += normal_lpdf(p_hat|p[adm2_index],sqrt(v));
    target += chi_square_lpdf(scaled_vhat|d);
    // mean model prior
    target += icar_normal_lpdf(u1|N, node1,node2); 
    target += normal_lpdf(u2|0,1);   
    target += normal_lpdf(sigma_u|0,1); //change this to penalised complexity prior 
    target += beta_lpdf(rho| 0.5,0.5); 
    // var model prior
    target += normal_lpdf(tau|0,1); //instead of target += normal_lpdf(log(v)|f, sigma_tau); 
    target += normal_lpdf(sigma_tau|0,1); //change this to penalised prior
    target += normal_lpdf(gamma0|0,1);
    target += normal_lpdf(gamma1|1,0.5);
    target += normal_lpdf(gamma2|-1, 0.5);
}
