#set number of chains
NCHAINS=8
NITER=3500
NWARMUP=3250

library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

recodeFACTOR<-function(FACTOR) {
  number<-as.numeric(as.factor(as.numeric(FACTOR)))
  should<-c(1:max(number))
  if(length(which(should %in% number==F))>0) {
    return("BROKE BROKE BROKE")
  }
  else {return(number)}
}

############BASIC MODEL 
############BASIC MODEL 

load("prior_list.Rdata")
load("TEN_squares.Rdata")

dm<-read.csv("landsat_lidar_ndvi.csv")

dms<-dm[sample(x=c(1:nrow(dm)),size=4000),]

ndvi=TEN_squares$ndvi

#this code is to process the missing values, so that there is an index for missing values
go=matrix(NA,nrow=nrow(ndvi),ncol=ncol(ndvi))

for(i in 1:nrow(ndvi)){
  for(j in 2:ncol(ndvi)){
    
    go[i,j]<-ifelse(is.nan(ndvi[i,j-1])==T & is.nan(ndvi[i,j])==T,"OUT","IN")
    
  }
}

consecutive_NAs=which(go=="OUT",arr.ind=T)[,1]
bad_NDVI=which(TEN_squares$ndvi<0,arr.ind=T)[1]

revised_BASE<-TEN_squares$BASE[-c(consecutive_NAs,bad_NDVI),]
revised_NDVI0<-TEN_squares$ndvi[-c(consecutive_NAs,bad_NDVI),]
revised_MM0<-TEN_squares$mm[-c(consecutive_NAs,bad_NDVI),]
revised_MM<-ifelse(is.nan(revised_MM0)==T,0,revised_MM0)/79
revised_NDVI<-ifelse(is.nan(revised_NDVI0)==T,0,revised_NDVI0)

basic_data_list<-prior_list

basic_data_list$n<-nrow(revised_NDVI)
basic_data_list$ID<-recodeFACTOR(as.factor(revised_BASE$SquareID))
basic_data_list$Nyear<-ncol(revised_NDVI)
basic_data_list$yr<-c(2001:2015)

#sigma priors
basic_data_list$base_sigma_m<-basic_data_list$base_sigma[1]
basic_data_list$base_sigma_sd<-basic_data_list$base_sigma[2]

#year sigma priors
basic_data_list$year_sigma_m<-basic_data_list$year_sigma[1]
basic_data_list$year_sigma_sd<-basic_data_list$year_sigma[2]

#mm sigma
basic_data_list$mm_sigma_m<-basic_data_list$mm_sigma[1]
basic_data_list$mm_sigma_sd<-basic_data_list$mm_sigma[2]

#mm slope priors
basic_data_list$mm_mm_m<-basic_data_list$mm_mm[1]
basic_data_list$mm_mm_sd<-basic_data_list$mm_mm[2]

#year slope priors
basic_data_list$year_mm_m<-basic_data_list$year_mm[1]
basic_data_list$year_mm_sd<-basic_data_list$year_mm[2]

#yearly effect
basic_data_list$mean_years=basic_data_list$year_year[,1]
basic_data_list$sd_years=basic_data_list$year_year[,2]
basic_data_list$dum_years=ifelse(basic_data_list$mean_years==0,0,1)

#missing data
basic_data_list$missers_row=which(is.nan(revised_NDVI0)==T,arr.ind=T)[,1]
basic_data_list$presents_row=which(is.nan(revised_NDVI0)==F,arr.ind=T)[,1]
basic_data_list$missers_col=which(is.nan(revised_NDVI0)==T,arr.ind=T)[,2]
basic_data_list$presents_col=which(is.nan(revised_NDVI0)==F,arr.ind=T)[,2]
basic_data_list$Nmis=length(basic_data_list$missers_col)
basic_data_list$Nobs=length(basic_data_list$presents_col)
basic_data_list$ndvi_obs<-revised_NDVI
basic_data_list$mm<-revised_MM


basic_data_list$mmLidar=(dms$mm-11)/79
basic_data_list$Nlidar=length(basic_data_list$mmLidar)
basic_data_list$ndvi_lidar=dms$ndvi
basic_data_list$HT=dms$HT
basic_data_list$a_HT<-basic_data_list$max_ht[1]
basic_data_list$b_HT<-basic_data_list$max_ht[2]
  
ndvi<-ifelse(basic_data_list$ndvi_obs==0,NA,basic_data_list$ndvi_obs)

#Run Stan code

BASIC<-stan(file="BASIC.stan",data =basic_data_list,
                        
               #chains=1,iter=25,warmup=5,
            chains=NCHAINS,iter=NITER,warmup=NWARMUP,
            control=list(adapt_delta=0.9975,max_treedepth=25),
                   pars = c("r1","K","sigma_mm","tau_int","mu_int","B_int",
                            "mu_r1","tau_r1","ndvi_true","ndvi_mis"))

#
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
############mm MODEL
MM<-stan(file="MM.stan",data =basic_data_list,
                                   
         #chains=1,iter=25,warmup=5,
         chains=NCHAINS,iter=NITER,warmup=NWARMUP,
                    control=list(adapt_delta=0.98,max_treedepth=15),
                                  pars = c("r1","K","sigma_mm","tau_int","mu_int","B_int","sigma",
                                   "mu_r1","tau_r1","ndvi_true","ndvi_mis","slope_mm","ndvi_pred"))

############YEAR
############YEAR
############YEAR
############YEAR
############YEAR

YEAR<-stan(file="YEAR.stan",data =basic_data_list,
           chains=1,iter=25,warmup=5,
           #chains=NCHAINS,iter=NITER,warmup=NWARMUP,
            control=list(adapt_delta=0.98,max_treedepth=15),
           pars = c("r1","K","sigma_mm","tau_int","mu_int","B_int","sigma",
                    "mu_r1","tau_r1","ndvi_true","ndvi_mis","slope_mm","year_effect","ndvi_pred"))

########HT
########HT
########HT
########HT
########HT
########HT

basic_HT<-stan(file="HT.stan",data=basic_data_list,
                #chains=1,iter=25,warmup=5,
               chains=NCHAINS,iter=NITER,warmup=NWARMUP,
                              control=list(adapt_delta=0.99,max_treedepth=15),
               pars = c("r1","K","sigma_mm","tau_int","mu_int","B_int","sigma",
                        "mu_r1","tau_r1","ndvi_true","ndvi_mis","slope_mm","year_effect","ndvi_pred",
                        "HT_intercept","b_LiDARmm","HT_NDVIslope"))

save(basic_HT,file="basic_HT.Rdata")