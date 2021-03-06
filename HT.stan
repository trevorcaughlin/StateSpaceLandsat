data {
  int<lower = 0> Nlidar;//Number of Canopy Height Model pixels from Lidar 
  vector[Nlidar] mmLidar; //Day of year for Landsat pixels in Landsat-lidar model
  real<lower = 0> HT[Nlidar];//Height from Canopy Height Model
  vector[Nlidar] ndvi_lidar;//NDVI values for Landsat-lidar model

  int<lower=1> Nyear; //number of years in annual time series
  int n; // number of pixels
  int yr[Nyear];//years in annual time series
  int ID[n]; //pixel identity in plot for random effect
  int Nmis;//number of missing observations
  int Nobs; //number of observed observations
  matrix[n,Nyear] ndvi_obs;//matrix of NDVI (including missing obs)
  matrix[n,Nyear] mm; //matrix of date-of-acquisition per pixel, to pair with ndvi_obs
  
  real<lower=0> a_HT; //priors for K (carrying capacity)
  real<lower=0> b_HT; //priors for K (carrying capacity)
  
  //priors for NDVI observation model
  real<lower=0> year_mm_sd; //prior for effect of date-of-acquisition on measurement error
  real<upper=0> year_mm_m; //prior for effect of date-of-acquisition on measurement error
  
  real<lower=0> year_sigma_m; //prior for standard deviation of NDVI due to measurement error
  real<lower=0> year_sigma_sd; //prior for standard deviation of NDVI due to measurement error
  
  int missers_row[Nmis]; //index for missing values
  int presents_row[Nobs]; //index for observed values
  int missers_col[Nmis]; //index for missing values
  int presents_col[Nobs]; //index for observed values
  
  real mean_years[Nyear]; //prior for effect of year, as categorical variable, on measurement error
  real sd_years[Nyear]; //prior for effect of year, as categorical variable, on measurement error
  
  int<lower=0,upper=1> dum_years[Nyear]; //dummy variable for whether year is year of lidar acquisition or not 
}

transformed data{
  int max_plots; //maximum number of plots for random effect
  
  max_plots = max(ID);
  
}

parameters {
  real HT_intercept;//intercept in gamma glm for Landsat-lidar fusion
  real<lower=0> HT_NDVIslope; //slope for effect of NDVI in Landsat-lidar fusion
  real year_effect[Nyear]; //effect of year as categorical variable
  real<lower=-1,upper=1> ndvi_mis[Nmis]; //missing values that will be interpolated
  real<lower=0> mu_r1; //mean growth rate
  real<lower=0> r1_raw[max_plots]; //random effect for growth rate (r)
  real B_int_raw[max_plots]; //random effect for initial starting NDVI
  real<lower=0> tau_r1;//non-centered standard deviation for growth random effect (r)
  real<lower=0> K; //carrying capacity; NDVI saturation
  real<lower=0> sigma; //standard deviation for process model
  real<lower=0> mu_int;//mean value of initial starting NDVI 
  real<lower=0> tau_int; //standard deviation for non-centered version of initial starting conditions
  real<lower=0,upper=1> sigma_mm; //measurement error 
  matrix<lower=0>[n,Nyear] ndvi_true; //NDVI as latent state variable
  real<upper=0> slope_mm; //effect of date-of-acquisition in measurement error model
  real<lower=0> b_LiDARmm; //effect of date-of-acquisition in Landsat-lidar fusion
  real<lower=0> gamma_shape; //shape parameter for gamma glm (Landsat-lidar fusion)
}

transformed parameters{
  
  matrix<lower=-1,upper=1>[n,Nyear] ndvi;//ndvi matrix as response variable
  matrix[n,Nyear] ndvi1; //matrix to back-transform Landsat-lidar fusion to spectral reflectance
  real<lower=0> r1[max_plots]; //random effect for growth rate (r)
  real<lower=0> B_int[max_plots]; //random effect for initial conditions

for(m in 1:Nmis){
ndvi[missers_row[m],missers_col[m]] = ndvi_mis[m];//interpolate missing values

}

for(p in 1:Nobs){
ndvi[presents_row[p],presents_col[p]] = ndvi_obs[presents_row[p],presents_col[p]];//let observed values remain untouched
}

  for(jj in 1:max_plots){ //non-centered parameterization. this helped convergence a lot!
    r1[jj] = mu_r1 + tau_r1* r1_raw[jj]; //random effect for growth
     B_int[jj] = mu_int + tau_int* B_int_raw[jj];  //random effect for initial starting conditions
  }

  

  for(pp in 1:n) {
for(qq in 1:Nyear) {    
     ndvi1[pp,qq] = (log(ndvi_true[pp,qq])-HT_intercept)/HT_NDVIslope; //transforming predicted HT to NDVI
        }
  } 



}

model { //for explanation of parameters, see parameter block
  HT_NDVIslope~normal(0,10); 
  b_LiDARmm~normal(0,10);
  HT_intercept~normal(0,5);
  slope_mm~normal(year_mm_m,year_mm_sd);
  sigma_mm~normal(year_sigma_m,year_sigma_sd);
  K ~ gamma(a_HT,b_HT);
  sigma~exponential(200);
  mu_r1 ~ normal(1,1);
  tau_r1 ~ exponential(5);
  tau_int~exponential(50);
  mu_int~normal(0,5);
  gamma_shape~normal(0,10);

for(p in 1:Nlidar){ //landsat-lidar fusion model--basically, a gamma glm

HT[p]~gamma(gamma_shape,(gamma_shape/exp(HT_intercept + HT_NDVIslope*ndvi_lidar[p] + b_LiDARmm*mmLidar[p]))); //log link to ensure positive values 

}
   
   for(kk in 1:max_plots){
   B_int_raw[kk]~normal(0,1);//drawing from standard normal for non-centered parameterization
   r1_raw[kk]~normal(0,1);
   }
   
for(k in 1:Nyear) {
  
  year_effect[k]~normal(mean_years[k],sd_years[k]); //effect of years drawn from informative priors
  
}

  for(j in 1:n) {
    
     ndvi_true[j,1] ~ normal(B_int[ID[j]],sigma);//initial condition of state variable
      
        ndvi[j,1] ~ normal(ndvi1[j,1]+ slope_mm*mm[j,1]+year_effect[1],sigma_mm); //initial condition of observed data
      
    for (i in 2:Nyear) {
   
     ndvi_true[j,i] ~ normal(ndvi_true[j,i-1]*(1+r1[ID[j]]*(1-ndvi_true[j,i-1]/K)),sigma); //process model
      
      ndvi[j,i]~normal(ndvi1[j,i] + slope_mm*mm[j,i]+dum_years[i]*year_effect[i],sigma_mm); //measurement error model
    }
  }
 
}
// 
// 
 generated quantities{ //posterior predictions for model validation and visualization

  matrix[n,Nyear] ndvi_pred; //predicted matrix of NDVI

  for(jj in 1:n){
    for(ii in 1:Nyear){

      ndvi_pred[jj,ii] = ndvi1[jj,ii] + slope_mm*mm[jj,ii]+dum_years[ii]*year_effect[ii];


    }

  }

}