//check on NDVI mm relationship? gauranteed to be positive?
data {
    int<lower = 0> Nlidar;//  pixels for L5 to LiDAR
      vector[Nlidar] mmLidar; //Lidar
  int<lower = 0> HT[Nlidar];// LiDAR pixels > ht 
 vector[Nlidar] ndvi_lidar;//Landsat greenness for LiDAR
}

parameters {
  real<upper=0> HT_intercept;
  real<lower=0> HT_NDVIslope;
  real<lower=0> b_Lidar;
  real<lower=0> b_LiDARmm;
   real<lower=0> gamma_parameter1;
    real<lower=0> gamma_scale2;
}

}

model {
exp(HT_intercept + HT_NDVIslope*ndvi_lidar + b_LiDARmm*mmLidar); //using the log link 

    HT_NDVIslope~cauchy(0,10);
    b_LiDARmm~cauchy(0,10);
   HT_intercept~cauchy(0,10);


 
//gamma regression
HT ~ gamma(gamma_parameter1,gamma_parameter2);


}
