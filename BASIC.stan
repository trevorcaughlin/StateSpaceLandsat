data {
  int<lower=1> Nyear; //number of years of annual time series
  int n; // number of pixels
  int yr[Nyear]; //years in annual time series
  int ID[n]; //pixel identity in plot for random effect
  int Nmis; //number of missing observations
  int Nobs; //number of observed observations
  matrix[n,Nyear] ndvi_obs; //matrix of NDVI (including missing obs)
  real<lower=0> a; //prior for K (carrying capacity)
  real<lower=0> b;//prior for K (carrying capacity)
  real<lower=0> base_sigma_m; //prior for standard deviation of NDVI due to measurement error
  real<lower=0> base_sigma_sd; //prior for standard deviation of NDVI due to measurement error
  int missers_row[Nmis]; //index for missing values
  int presents_row[Nobs]; //index for observed values
  int missers_col[Nmis]; //index for missing values
  int presents_col[Nobs]; //index for observed values
}


transformed data{
  int max_plots;
  
  max_plots = max(ID); //maximum number of plots for random effect
  
}

parameters {
  real<lower=-1,upper=1> ndvi_mis[Nmis]; //missing values that will be interpolated
  real<lower=0> mu_r1; //mean growth rate
  real<lower=0> r1_raw[max_plots]; //random effect for growth rate (r)
  real B_int_raw[max_plots]; //random effect for initial starting NDVI
  real<lower=0> tau_r1; //non-centered standard deviation for growth random effect (r)
  real<lower=0.6,upper=1> K; //carrying capacity; NDVI saturation
  real<lower=0,upper=1> mu_int; //mean value of initial starting NDVI 
  real<lower=0> tau_int; //standard deviation for non-centered version of initial starting conditions
  real<lower=0,upper=1> sigma_mm; //measurement error 
  matrix<lower=-1,upper=1>[n,Nyear] ndvi_true; //NDVI as latent state variable
}

transformed parameters{
  
  matrix<lower=-1,upper=1>[n,Nyear] ndvi; //ndvi matrix as response variable
  real<lower=0> r1[max_plots]; //random effect for growth rate (r)
  real B_int[max_plots]; //random effect for initial conditions

for(m in 1:Nmis){
ndvi[missers_row[m],missers_col[m]] = ndvi_mis[m]; //interpolate missing values
}

for(p in 1:Nobs){
ndvi[presents_row[p],presents_col[p]] = ndvi_obs[presents_row[p],presents_col[p]]; //let observed values remain untouched
}

  for(jj in 1:max_plots){ //non-centered parameterization. this helped convergence a lot!
    r1[jj] = mu_r1 + tau_r1* r1_raw[jj]; //random effect for growth
     B_int[jj] = mu_int + tau_int* B_int_raw[jj]; //random effect for initial starting conditions
  }

  

}

model {
    sigma_mm~normal(base_sigma_m,base_sigma_sd); //standard deviation of measurement error 
    K ~ beta(a,b); //"carrying capacity" or NDVI saturation
    mu_r1 ~ normal(1,1); //mean growth rate
    tau_r1 ~ exponential(5); //non-centered standard deviation for growth rate (r)
    tau_int~exponential(50);  //non-centered standard deviation for initial condition
    mu_int~beta(10,10); //mean initial condition
   
   for(kk in 1:max_plots){
   B_int_raw[kk]~normal(0,1); //drawing from standard normal for non-centered parameterization
   r1_raw[kk]~normal(0,1);
   }


  for(j in 1:n) {
    
     ndvi_true[j,1] ~ normal(B_int[ID[j]],0.01); //initial condition of state variable
      
        ndvi[j,1] ~ normal(ndvi_true[j,1],sigma_mm); //initial condition of observed data
      
    for (i in 2:Nyear) {
   
     ndvi_true[j,i] ~ normal(ndvi_true[j,i-1]*(1+r1[ID[j]]*(1-ndvi_true[j,i-1]/K)),0.01); //process model
      
      ndvi[j,i]~normal(ndvi_true[j,i],sigma_mm); //measurement model
    }
  }
 
}