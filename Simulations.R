library(lme4)
source("Empi_Null_fast.R")

#Define parameter values
zeta=0
beta=1
sig_alpha=0.25
sig_gamma=0.25
sig_psi=0.25
n_center=10
n_OPO=10
center=1:n_center
OPO=1:n_OPO
n_combo=length(rep(rep(n_center,n_OPO)*(OPO-1),n_center)+rep(center-1,n_OPO))
n_pat=1000
n_sim=100
#Main parameter of interest (0=null):
psi_size=0 #connection effect for first connection in dataset

sim_flag=function(O,E,alpha,alpha_sd,gamma,gamma_sd,n_sim=1000) {
  alpha_sim=rnorm(n_sim,alpha,alpha_sd)
  gamma_sim=rnorm(n_sim,gamma,gamma_sd)
  O_sim=rpois(n_sim,E*exp(alpha_sim+gamma_sim))
  temp=mean(O_sim > O)+0.5*mean(O_sim==O)
  mid_P=2*min(temp,1-temp)
  flag=as.numeric(mid_P < 0.05)
  lower=O/E-1.96*sd(O_sim/E)
  upper=O/E+1.96*sd(O_sim/E)
  return(list(flag=flag,lower=lower,upper=upper))
}

flag_naive=vector("numeric")
flag_sim=vector("numeric")
flag_EN=vector("numeric")
coverage_naive=vector("numeric")
coverage_sim=vector("numeric")
coverage_EN=vector("numeric")
for(i in 1:n_sim) {
  show(i)
  
  alpha=rnorm(n_center,0,sig_alpha)
  gamma=rnorm(n_OPO,0,sig_gamma)
  psi=rnorm(n_center*n_OPO,0,sig_psi)
  psi[1]=psi_size
  
  #Randomly select providers combinations for each patient
  #Generate covariate
  center_pat=sample(center,n_pat,replace=T)
  OPO_pat=sample(OPO,n_pat,replace=T)
  alpha_pat=alpha[center_pat]
  gamma_pat=gamma[OPO_pat]
  combo_pat=n_center*(OPO_pat-1)+center_pat
  psi_pat=psi[combo_pat]
  X_pat=rnorm(n_pat)
  
  theta=zeta+beta*X_pat+alpha_pat+gamma_pat+psi_pat
  Y=rpois(n_pat,exp(theta))
  
  #Fit model
  fit=glmer(Y~X_pat+(1|center_pat)+(1|OPO_pat)+(1|combo_pat),family="poisson")
  zeta_hat=fixef(fit)[1]
  beta_hat=fixef(fit)[2]
  re=as.data.frame(ranef(fit))
  alpha_hats=re[re$grpvar=="center_pat",]
  gamma_hats=re[re$grpvar=="OPO_pat",]
  
  df=data.frame(Y,combo_pat,gamma_pat,alpha_pat,X_pat,center_pat,OPO_pat)
  df=df[order(df$combo_pat),]
  
  O_ij=sapply(split(df$Y,df$combo_pat),sum)
  E_ij_Sim=sapply(split(exp(zeta_hat+beta_hat*df$X_pat+
                             alpha_hats$condval[match(df$center_pat,alpha_hats$grp)]+
                             gamma_hats$condval[match(df$OPO_pat,gamma_hats$grp)]),df$combo_pat),sum)
  E_ij_EN=sapply(split(exp(zeta_hat+beta_hat*df$X_pat),df$combo_pat),sum)
  
  #Naive Flagging
  S=(O_ij-E_ij_EN)/sqrt(E_ij_EN)
  flag_naive[i]=as.numeric(S[1] > 1.96 | S[1] < -1.96) 
  lower_naive=O_ij/E_ij_EN-1.96*sqrt(1/E_ij_EN)
  upper_naive=O_ij/E_ij_EN+1.96*sqrt(1/E_ij_EN)
  coverage_naive[i]=as.numeric(lower_naive[1] < exp(psi_size) & upper_naive[1] > exp(psi_size))
  
  #Simulation-Based Flagging
  sim_flagging=sim_flag(O_ij[1],E_ij_EN[1],
                        alpha_hats$condval[alpha_hats$grp==df$center_pat[1]],
                        alpha_hats$condsd[alpha_hats$grp==df$center_pat[1]],
                        gamma_hats$condval[gamma_hats$grp==df$OPO_pat[1]],
                        gamma_hats$condsd[gamma_hats$grp==df$OPO_pat[1]])
  flag_sim[i]=sim_flagging$flag
  coverage_sim[i]=as.numeric(sim_flagging$lower < exp(psi_size) & sim_flagging$upper > exp(psi_size))
  
  #Empirical Null Flagging
  a_id=df$center_pat[match(names(O_ij),df$combo_pat)]
  g_id=df$OPO_pat[match(names(O_ij),df$combo_pat)]
  EN_est=empirical_null_ind_RQ_corr(S, E_ij_EN, a_id, g_id)
  gt_90 <- function(x){ifelse(x<ceiling(quantile(x,0.9)), x, ceiling(quantile(x,0.9)))}
  S_EN=S/sqrt(1+(EN_est$phi^2)*gt_90(E_ij_EN))
  flag_EN[i]=as.numeric(S_EN[1] > 1.96 | S_EN[1] < -1.96) 
  lower_EN=O_ij/E_ij_EN-1.96*sqrt((1+(EN_est$phi^2)*gt_90(E_ij_EN))/E_ij_EN)
  upper_EN=O_ij/E_ij_EN+1.96*sqrt((1+(EN_est$phi^2)*gt_90(E_ij_EN))/E_ij_EN)
  coverage_EN[i]=as.numeric(lower_EN[1] < exp(psi_size) & upper_EN[1] > exp(psi_size))
}

#If psi_size=0, the simulation and empirical null methods should be close to 5% 
mean(flag_naive,na.rm=T) 
mean(flag_sim,na.rm=T)
mean(flag_EN,na.rm=T)







