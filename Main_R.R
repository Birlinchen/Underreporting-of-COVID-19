library(rstan)
library(bayesplot)
library(rstanarm)
library(shinystan)
library(ngspatial)
library(spdep)
library(parallel)
library(tidyverse)
library(nimble)
library(ggplot2)
library(INLA)

options("scipen"=100, "digits"=4)

###################souring functions#############
source("Functions.R")
source("StanCode_Covid19.R")

#readr::write_csv(dat_joint,file="April2020Covid_49states.csv")
dat = readr::read_csv("April2020Covid_49states.csv")
#read shape file
shp2 <- read_sf("shp5m/cb_2018_us_state_5m.shp")
shp2_48andDC <- shp2[-c(40,48,49,50,52,56,28),]

#constructing neighborhood relationship
rook.nb = poly2nb(shp2_48andDC,queen = F)
A = nb2mat(rook.nb,style = "B")
n_adj=rowSums(A)
D=diag(n_adj)
adj=as.carAdjacency(A)$adj
nbs = mungeCARdata4stan(adj, n_adj)
N = nbs$N
node1 = nbs$node1
node2 = nbs$node2
N_edges = nbs$N_edges
#Moran's I test 
lw=nb2listw(rook.nb, style="W", zero.policy=TRUE)
spdep::moran.test(dat$positive,lw)
MC = spdep::moran.mc(dat$positive,lw, nsim=599)
plot(MC, main="", las=1)

##################compute scaling_factor#############
#Build the adjacency matrix using INLA library functions
adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
#The ICAR precision matrix (note! This is singular)
Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)
# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
#See the inla.qinv function help for further details.
Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))
#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor = exp(mean(log(diag(Q_inv))))

CoresToUse = parallel::detectCores()-3
RNGkind("L'Ecuyer-CMRG")
set.seed(2021)

###############  Fit naive model ##################
X = dat %>% select(
  #Popdensity,
  Uninsured,
  #Inactivity,
  Obesity,
  #Smoking,
  Excessive_Drinking,
  AirPol,
  Drug_death,
  MDI)
X_scaled = scale(X)

df= X_scaled %>% as.data.frame()
df$y = dat$positive
df$E = dat$Pop
m1 = glm(y ~ Uninsured + Obesity + Excessive_Drinking + AirPol +
           Drug_death + MDI + offset(log(E)), data=df, family = poisson(link = "log"))

naive_data = list(N=N,
                  z=dat$positive,
                  E=dat$Pop,
                  K=ncol(X),
                  X=X_scaled)
naive_samples = rstan::sampling(compiled_naive, 
                                data=naive_data, 
                                warmup=2000, 
                                iter=8000, 
                                chains=2, 
                                thin=1, 
                                cores = 1,
                                seed = 2021,
                                #pars = c('deltas','gammas','tau','tauh','y_rep'),
                                #pars = c('deltas','gammas','tau','tauh',
                                #         'beta0','betas','z_rep','y_rep','p_rep'),
                                #pars=c('beta0','betas','gamma0','gammas'),
                                control = list(adapt_delta = 0.85, 
                                               max_treedepth = 15))
summa_Naive_fits = rstan::summary(naive_samples,
                                  permuted = T,
                                  pars = c("gamma0", "gammas"
                                  ),
                                  probs=c(0.05, 0.95))$summary
summa_Naive_fits
readr::write_csv(summa_Naive_fits %>% as.data.frame(),"Naive.csv")
#############Under only######################################
testing <- (dat$totalTestResults)/dat$Pop*1000
#W <- cbind(testing, df %>% dplyr::select(Physician,Uninsured))
W = matrix(testing,ncol = 1)
J <- ncol(W)
UnderOnly_data = list(N=N,
                      z=dat$positive,
                      E=dat$Pop,
                      K=ncol(X),
                      X=X_scaled,
                      W =scale(W),
                      J=ncol(scale(W)))
UnderOnly_samples = rstan::sampling(compiled_Under_Non_Spatial, 
                                    data=UnderOnly_data, 
                                    warmup=1000, 
                                    iter=2000, 
                                    chains=2, 
                                    thin=1, 
                                    cores = 1,
                                    seed = 2021,
                                    #pars = c('deltas','gammas','tau','tauh','y_rep'),
                                    #pars = c('deltas','gammas','tau','tauh',
                                    #         'beta0','betas','z_rep','y_rep','p_rep'),
                                    #pars=c('beta0','betas','gamma0','gammas'),
                                    control = list(adapt_delta = 0.85, 
                                                   max_treedepth = 15))
summa_UnderOnly_fits = rstan::summary(UnderOnly_samples,
                                      permuted = T,
                                      pars = c("gamma0","gammas", "beta0","betas"
                                      ),
                                      probs=c(0.05, 0.95))$summary
summa_UnderOnly_fits
readr::write_csv(summa_UnderOnly_fits %>% as.data.frame(),"Underonly.csv")
##################Spatial only#################################
data_SpatialOnly = list(N=N,
                        z=dat$positive,
                        E=dat$Pop,
                        K=ncol(X),
                        X=X_scaled,
                        N_edges=N_edges,
                        node1=node1,
                        node2=node2,
                        scaling_factor=scaling_factor)         
SpatialOnly_sample = rstan::sampling(compiled_Spa_Non_Under, 
                                     data=data_SpatialOnly, 
                                     warmup=1000, 
                                     iter=2000, 
                                     chains=2, 
                                     thin=1, 
                                     cores = 1,
                                     #pars = c('deltas','gammas','tau','tauh','y_rep'),
                                     #pars = c('deltas','gammas','tau','tauh',
                                     #         'beta0','betas','z_rep','y_rep','p_rep'),
                                     #pars=c('beta0','betas','gamma0','gammas'),
                                     control = list(adapt_delta = 0.85, 
                                                    max_treedepth = 15))
summa_SpatialOnly_fits = rstan::summary(SpatialOnly_sample,
                                        permuted = T,
                                        pars = c("gamma0","gammas"
                                        ),
                                        probs=c(0.05, 0.95))$summary
summa_SpatialOnly_fits
readr::write_csv(summa_SpatialOnly_fits %>% as.data.frame(),"Spatialonly.csv")
###############Fit full_BYM2 model#################

#Construct scaled risk factors for Poisson part
X = dat %>% select(
  #Popdensity,
  Uninsured,
  #Inactivity,
  Obesity,
  #Smoking,
  Excessive_Drinking,
  AirPol,
  Drug_death,
  MDI)
res <- cor(X)
round(res, 2)
#X_scaled = scale(X%>% select(Obesity,Drug_death))
X_scaled = scale(X)
stan_data=list(N=N,
               N_edges=N_edges,
               node1=node1,
               node2=node2,
               J=J,
               z=dat$positive,
               K=ncol(X_scaled ),
               x=X_scaled,
               w = scale(W),
               E= dat$Pop,
               scaling_factor=scaling_factor
)
sampling2 <- rstan::sampling(compiled_full_BYM2, data=stan_data, warmup=1000, 
                      iter=2000, chains=2, thin=1, cores = 2,
                      #pars = c('deltas','gammas','tau','tauh','y_rep'),
                      #pars = c('deltas','gammas','tau','tauh',
                      #         'beta0','betas','z_rep','y_rep','p_rep'),
                      #pars=c('beta0','betas','gamma0','gammas'),
                      control = list(adapt_delta = 0.85, 
                                     max_treedepth = 15))
summa2 = rstan::summary(sampling2,
                        permuted = T,
                        pars=c("gamma0","gammas","beta0",'betas'),probs=c(0.05, 0.95))
summa2$summary
readr::write_csv(summa2$summary %>% as.data.frame(),"Full.csv")
###############################
summa4models = rbind(summa_Naive_fits,
                     summa_SpatialOnly_fits,
                     summa_UnderOnly_fits,
                     summa2$summary) %>% as.data.frame()
row.names(summa4models) <- c()
paras = c("Intercept_Pois","uninsured","Obesity","Excessive_Drinking","AirPol","Drug_death","MDI",
          "Intercept_Pois","uninsured","Obesity","Excessive_Drinking","AirPol","Drug_death","MDI",
          "Intercept_Pois","uninsured","Obesity","Excessive_Drinking","AirPol","Drug_death","MDI","Intercept_logit","Testing",
          "Intercept_Pois","uninsured","Obesity","Excessive_Drinking","AirPol","Drug_death","MDI","Intercept_logit","Testing")
models = c(rep("M1",7), rep("M3",7), rep("M2",9),rep("M4",9))
summa4models$paras = paras
summa4models$model = models

ggplot(summa4models, aes(x=model, y=mean,colour=model)) + 
  geom_errorbar(aes(ymin=`5%`, ymax=`95%`), width=.2,size=1) + 
  #geom_line() + 
  geom_hline(yintercept=0, linetype="dashed",colour="black") +
  geom_point(size=2) + facet_wrap(vars(paras),ncol=3,scales = "free") + 
  theme_bw() + theme(legend.position = "none") + ylab("value")

#####################get 5 most and 5 least states
FiveMostFiveLeast =rstan::summary(sampling2,
                                  permuted = T,
                                  pars=c("y_rep"),probs=c(0.025, 0.975))$summary
FiveMostFiveLeast = FiveMostFiveLeast %>% as.data.frame()
FiveMostFiveLeast$State = dat$State
FiveMostFiveLeast$Reported = dat$positive
FiveMostFiveLeast=FiveMostFiveLeast[order(FiveMostFiveLeast$mean),]

TenStates = FiveMostFiveLeast[c(1:5,45:49),]
readr::write_csv(TenStates,"TenStates.csv")
###############################
row.names(summa2$summary) = c("intercept1", "Testing", "intercept2",
                              "Uninsured",
                              "Obesity","Drinking","AirPol","Drug_death",
                              "MDI","logit_rho","sigma","lp__")
summ = summa2$summary %>% round(digits=4)

summa_pi =  rstan::summary(sampling2,
                           permuted = T,
                           pars=c("p"),probs=c(0.05, 0.95))$summary
summa_pi

##############################plot the fitting results#############
pi = summa_pi[,1] %>% round(digits = 2)
pi_samples <- rstan::extract(sampling2, 
                             pars = 'p', permuted = F)
color_scheme_set("blue")
pars1 = rep(NA,12)
for (i in 1:12) {
  pars1[i] = stringr::str_c("p[",i,"]")
}

q1=mcmc_intervals(pi_samples,pars = pars1 ,
                  prob = 0.5,
                  prob_outer = 0.95,
                  point_est = c("median"),
                  inner_size = 1,
                  point_size = 2)+ 
  ggplot2::scale_x_continuous(limits=c(0,1),breaks = seq(0,1,length.out = 5))

q1=q1+ggplot2::scale_y_discrete(
  labels = dat$Abbre[13:24]
)+
  theme(axis.title=element_text(size=10))
q1
#########
pars2 = rep(NA,12)
for (i in 1:12) {
  pars2[i] = stringr::str_c("p[",i+12,"]")
}
q2=mcmc_intervals(pi_samples,pars = pars2 ,
                  prob = 0.5,
                  prob_outer = 0.95,
                  point_est = c("median"),
                  inner_size = 1,
                  point_size = 2)+ 
  ggplot2::scale_x_continuous(limits=c(0,1),breaks = seq(0,1,length.out = 5))

q2=q2+ggplot2::scale_y_discrete(
  labels = dat$Abbre[13:24]
)+
  theme(axis.title=element_text(size=10))
q2
#########
pars3 = rep(NA,12)
for (i in 1:12) {
  pars3[i] = stringr::str_c("p[",i+24,"]")
}
q3=mcmc_intervals(pi_samples,pars = pars3,
                  prob = 0.5,
                  prob_outer = 0.95,
                  point_est = c("median"),
                  inner_size = 1,
                  point_size = 2)+ ggplot2::scale_x_continuous(limits=c(0,1),
                                                               breaks = seq(0,1,length.out = 5))

q3=q3+ggplot2::scale_y_discrete(
  labels = dat$Abbre[25:36]
)+#ggplot2::xlab("Point estimate and 95% CI for reporting rate")+
  theme(axis.title=element_text(size=6))
q3
#########
pars4 = rep(NA,13)
for (i in 1:13) {
  pars4[i] = stringr::str_c("p[",i+36,"]")
}
q4=mcmc_intervals(pi_samples,pars = pars4,
                  prob = 0.5,
                  prob_outer = 0.95,
                  point_est = c("median"),
                  inner_size = 1,
                  point_size = 2) + ggplot2::scale_x_continuous(limits=c(0,1),
                                                               breaks = seq(0,1,length.out = 5))

q4=q4+ggplot2::scale_y_discrete(
  labels = dat$Abbre[37:49]
)+#ggplot2::xlab("Point estimate and 95% CI for reporting rate")+
  theme(axis.title=element_text(size=6))
q4

q=gridExtra::grid.arrange(q1,q2,q3,q4,nrow=2,ncol=2)
