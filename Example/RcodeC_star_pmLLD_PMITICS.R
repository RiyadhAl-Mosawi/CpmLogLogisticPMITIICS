##--------------------------------------------------------------------------------
# The R and C++ code for the paper 
# Estimation of $\mathcal C*_pm$ under Progressive 
# Multiple Interval Type-I Censoring: Bayesian and Classical 
# Approaches with Log-Logistic Distribution
# Jointly with : Sanku Dey
# We compute MLE, MPS and Bayes MCMC  based on LK and PS functions
# Written by:   Dr. Riyadh Al-Mosawi
# Date 11-5-2026
##--------------------------------------------------------------------------------

rm(list=ls())
cat("\014")

# --- Load or Install Packages ---
load_or_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

for (pkg in c("Rcpp", "pracma", "roptim", "coda", "progress", "xtable", "knitr", "dplyr")) {
  load_or_install(pkg)
}

options(width=100,length=200,max.print = 10000,digits = 4, scipen=999)
#path="//MASTER-COMPUTER/UnderworkPapers/CsMITIICloglog/Rcode"
#path="C:/Users/RiyadhRustam/Desktop/CurrentFolder/UnderworkPapers/CsMITIICloglog/Rcode"
path="//Master-Computer/UnderworkPapers/Estimation of C∗pm under PMITICS with LogLogisticD/RcodePMITICS"
path="C:/Users/hp/Riyadhrm Dropbox/Riyadh/Riyadh_AlMosawi/CurrentFolder/UnderworkPapers/Estimation of C∗pm under PMITICS with LogLogisticD/RcodePMITICS"
setwd(path)

#setwd="C:/Users/hp/Riyadhrm Dropbox/Riyadh/Riyadh_AlMosawi/CurrentFolder/UnderworkPapers/CsMITIICloglog/Rcode/CsLLD"
#ls("package:CsLogLogistic")

#sourceCpp("UserFunLogLogistPMITICS_cpp.cpp")
library(CpmLogLogisticPMITICS)

#-- Defining the Macro variables
DataSet <- function(index_set = 1) {
  vals <- matrix(c(
    1.2,2,50,6 ,0.3,0.20,0.1,9,4,1,
    
    # 1.2,2,50,6 ,0.3,0.30,0.1,9,4,1,1,
    # 1.2,2,50,6 ,0.3,0.30,0.1,9,4,1,2,
    # 1.2,2,50,6 ,0.3,0.30,0.1,9,4,1,3,
    # 1.2,2,50,6 ,0.3,0.30,0.1,9,4,1,4,
    # 1.2,2,50,6 ,0.3,0.30,0.1,9,4,1,5,
    
    1.2,2,50,6 ,0.3,0.60,0.1,9,4,1,
    
    # 1.2,2,50,6 ,0.7,0.20,0.1,9,4,1,1,
    # 1.2,2,50,6 ,0.7,0.20,0.1,9,4,1,2,
    # 1.2,2,50,6 ,0.7,0.20,0.1,9,4,1,3,
    # 1.2,2,50,6 ,0.7,0.20,0.1,9,4,1,4,
    # 1.2,2,50,6 ,0.7,0.20,0.1,9,4,1,5,
    # 
    # 1.2,2,50,6 ,0.7,0.30,0.1,9,4,1,1,
    # 1.2,2,50,6 ,0.7,0.30,0.1,9,4,1,2,
    # 1.2,2,50,6 ,0.7,0.30,0.1,9,4,1,3,
    # 1.2,2,50,6 ,0.7,0.30,0.1,9,4,1,4,
    # 1.2,2,50,6 ,0.7,0.30,0.1,9,4,1,5,
    
    # 1.2,2,50,6 ,0.7,0.60,0.1,9,4,1,1,
    # 1.2,2,50,6 ,0.7,0.60,0.1,9,4,1,2,
    # 1.2,2,50,6 ,0.7,0.60,0.1,9,4,1,3,
    # 1.2,2,50,6 ,0.7,0.60,0.1,9,4,1,4,
    # 1.2,2,50,6 ,0.7,0.60,0.1,9,4,1,5,
    
    1.2,2,50,10,0.3,0.20,0.1,9,4,1,
    
    # 1.2,2,50,10,0.3,0.30,0.1,9,4,1,1,
    # 1.2,2,50,10,0.3,0.30,0.1,9,4,1,2,
    # 1.2,2,50,10,0.3,0.30,0.1,9,4,1,3,
    # 1.2,2,50,10,0.3,0.30,0.1,9,4,1,4,
    # 1.2,2,50,10,0.3,0.30,0.1,9,4,1,5,
    
    1.2,2,50,10,0.3,0.60,0.1,9,4,1,
 
    # 1.2,2,50,10,0.7,0.20,0.1,9,4,1,1,
    # 1.2,2,50,10,0.7,0.20,0.1,9,4,1,2,
    # 1.2,2,50,10,0.7,0.20,0.1,9,4,1,3,
    # 1.2,2,50,10,0.7,0.20,0.1,9,4,1,4,
    # 1.2,2,50,10,0.7,0.20,0.1,9,4,1,5,
    
    # 1.2,2,50,10,0.7,0.30,0.1,9,4,1,1,
    # 1.2,2,50,10,0.7,0.30,0.1,9,4,1,2,
    # 1.2,2,50,10,0.7,0.30,0.1,9,4,1,3,
    # 1.2,2,50,10,0.7,0.30,0.1,9,4,1,4,
    # 1.2,2,50,10,0.7,0.30,0.1,9,4,1,5,
    
    
    # 1.2,2,50,10,0.7,0.60,0.1,9,4,1,1,
    # 1.2,2,50,10,0.7,0.60,0.1,9,4,1,2,
    # 1.2,2,50,10,0.7,0.60,0.1,9,4,1,3,
    # 1.2,2,50,10,0.7,0.60,0.1,9,4,1,4,
    # 1.2,2,50,10,0.7,0.60,0.1,9,4,1,5,
    
    1.2,2,100,6 ,0.3,0.20,0.1,9,4,1,
 
    # 1.2,2,100,6 ,0.3,0.30,0.1,9,4,1,1,
    # 1.2,2,100,6 ,0.3,0.30,0.1,9,4,1,2,
    # 1.2,2,100,6 ,0.3,0.30,0.1,9,4,1,3,
    # 1.2,2,100,6 ,0.3,0.30,0.1,9,4,1,4,
    # 1.2,2,100,6 ,0.3,0.30,0.1,9,4,1,5,
    
    1.2,2,100,6 ,0.3,0.60,0.1,9,4,1,
    
    # 1.2,2,100,6 ,0.7,0.20,0.1,9,4,1,1,
    # 1.2,2,100,6 ,0.7,0.20,0.1,9,4,1,2,
    # 1.2,2,100,6 ,0.7,0.20,0.1,9,4,1,3,
    # 1.2,2,100,6 ,0.7,0.20,0.1,9,4,1,4,
    # 1.2,2,100,6 ,0.7,0.20,0.1,9,4,1,5,
    
    # 1.2,2,100,6 ,0.7,0.60,0.1,9,4,1,1,
    # 1.2,2,100,6 ,0.7,0.60,0.1,9,4,1,2,
    # 1.2,2,100,6 ,0.7,0.60,0.1,9,4,1,3,
    # 1.2,2,100,6 ,0.7,0.60,0.1,9,4,1,4,
    # 1.2,2,100,6 ,0.7,0.60,0.1,9,4,1,5,
    
    # 1.2,2,100,6 ,0.7,0.90,0.1,9,4,1,1,
    # 1.2,2,100,6 ,0.7,0.90,0.1,9,4,1,2,
    # 1.2,2,100,6 ,0.7,0.90,0.1,9,4,1,3,
    # 1.2,2,100,6 ,0.7,0.90,0.1,9,4,1,4,
    # 1.2,2,100,6 ,0.7,0.90,0.1,9,4,1,5,
    
    1.2,2,100,10,0.3,0.20,0.1,9,4,1,
    
    1.2,2,100,10,0.3,0.60,0.1,9,4,1
    
    # 1.2,2,100,10,0.3,0.90,0.1,9,4,1,1,
    # 1.2,2,100,10,0.3,0.90,0.1,9,4,1,2,
    # 1.2,2,100,10,0.3,0.90,0.1,9,4,1,3,
    # 1.2,2,100,10,0.3,0.90,0.1,9,4,1,4,
    # 1.2,2,100,10,0.3,0.90,0.1,9,4,1,5,
    
    # 1.2,2,100,10,0.7,0.20,0.1,9,4,1,1,
    # 1.2,2,100,10,0.7,0.20,0.1,9,4,1,2,
    # 1.2,2,100,10,0.7,0.20,0.1,9,4,1,3,
    # 1.2,2,100,10,0.7,0.20,0.1,9,4,1,4,
    # 1.2,2,100,10,0.7,0.20,0.1,9,4,1,5,
    
    # 1.2,2,100,10,0.7,0.30,0.1,9,4,1,1,
    # 1.2,2,100,10,0.7,0.30,0.1,9,4,1,2,
    # 1.2,2,100,10,0.7,0.30,0.1,9,4,1,3,
    # 1.2,2,100,10,0.7,0.30,0.1,9,4,1,4,
    # 1.2,2,100,10,0.7,0.30,0.1,9,4,1,5,
    
    # 1.2,2,100,10,0.7,0.60,0.1,9,4,1,1,
    # 1.2,2,100,10,0.7,0.60,0.1,9,4,1,2,
    # 1.2,2,100,10,0.7,0.60,0.1,9,4,1,3,
    # 1.2,2,100,10,0.7,0.60,0.1,9,4,1,4,
    # 1.2,2,100,10,0.7,0.60,0.1,9,4,1,5
    ), ncol = 10, byrow = TRUE)
  
  colnames(vals) <- c(
    "bet", "gam",
    "n", "m",
    "censor_prop","l_prop",
    "l", "u", "t","lambda"
  )
  
  if(index_set < 1 || index_set > nrow(vals)) {
    stop("index_set must be between 1 and ", nrow(vals), ".")
  }
  
  out <- vals[index_set, ]
  return(out)
}

scheme = function(index_scheme = 1) {
  vals <- matrix(c(rep(1,2*m),
                   rep(2,m),rep(0,m),
                   rep(0,m),rep(2,m),2*m,rep(0,2*m-1),rep(0,2*m-1),2*m),nrow=5,ncol=2*m,byrow=TRUE)
  
  rownames(vals) = c("R1","R2","R3","R4","R5")
  out <- vals[index_scheme, ]
  return(out)
}

save_bayes <- function(ARR, est, par_true, it, HPD = NULL) {
  
  for(j in 1:3) {
    
    if(is.null(HPD)) {
      ARR[j,it,] <- c(
        est[j],
        est[j] - par_true[j],
        (est[j] - par_true[j])^2,
        NA, NA, NA, NA, NA
      )
    } else {
      ARR[j,it,] <- c(
        est[j],
        est[j] - par_true[j],
        (est[j] - par_true[j])^2,
        NA,
        HPD[j,]
      )
    }
  }
  
  return(ARR)
}
#library(devtools)

#devtools::build("C:/Users/hp/Riyadhrm Dropbox/Riyadh/Riyadh_AlMosawi/CurrentFolder/UnderworkPapers/CsMITIICloglog/Rcode/CsLogLogistic")

#-- MCMC size and burn-in sample
MC_size = 20500  
MC_burn = 500

# LINEX parameter
c_val1 = -1.5
c_val2 =  1.5

Sim_no = 500
B = 1000 # bootstrap sample size

#--------------------------------------

index_set=8
para_set = DataSet(index_set)
bet = as.numeric(para_set["bet"])
gam = as.numeric(para_set["gam"])
n   = as.numeric(para_set["n"])
m   = as.numeric(para_set["m"])
censor_prop  = as.numeric(para_set["censor_prop"])
l_prop = as.numeric(para_set["l_prop"])
Cpm_par   = as.numeric(para_set[c("l","u","t","lambda")])

par = c(bet, gam)
par_true = c(bet, gam, Cpm_star_fun_cpp(c(bet, gam),Cpm_par))

lower_lim = c(0.01,0.01)
fact = 2
upper_lim = fact*par_true[1:2]
h_non = c(0.001,0.001,0.001,0.001)


for(index_scheme in 3:3){
  #Change the value of j from 1,2,3,...,18
  #--------------------------------------

  Rprog =scheme(index_scheme)
  
  randomseed=c(2021,6,30) 
  set.seed(randomseed)
  
  h_inf <- as.numeric(HyperPara_PMITICS_cpp(
    para     = par_true[1:2],
    n        = n,
    m        = m,
    censor_prop  = censor_prop,
    l_prop       = l_prop,
    Rprog        = Rprog,
    lower    = lower_lim,
    upper    = upper_lim,
    Hyper_It = 1000,
    method   = "Nelder-Mead"
  ))
  
  it=1
  
  while(it<=Sim_no){
    if(it==1){
      MLE=MPS=array(dim=c(3,Sim_no,8),NA) 
      # 1st dim=parameter(beta,gam,C_star_pm), 
      # 2nd dim=iteration, 
      # 3rd dim=(Est, AB, MSE, ESE, (L,U,AL,CP) ACI)
      BOOT_p_MLE = BOOT_p_MPS =array(dim=c(3,Sim_no,8),NA) 
      BOOT_t_MLE = BOOT_t_MPS =array(dim=c(3,Sim_no,8),NA) 
      BOOT_log_t_MLE = BOOT_log_t_MPS =array(dim=c(3,Sim_no,8),NA) 
      # 3rd dim=(NA,Na,Na,Na (L,U,AL,CP) ACI)
      
      MCMC_SEL_MLE_non=MCMC_SEL_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      MCMC_SEL_MPS_non=MCMC_SEL_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      # 3rd dim=(Est, AB, MSE, NA, (L,U,AL,CP) HPD)
      MCMC_LIN1_MLE_non=MCMC_LIN1_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      MCMC_LIN1_MPS_non=MCMC_LIN1_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      MCMC_LIN2_MLE_non=MCMC_LIN2_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      MCMC_LIN2_MPS_non=MCMC_LIN2_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      
      TK_SEL_MLE_non=TK_SEL_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      TK_SEL_MPS_non=TK_SEL_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      TK_LIN1_MLE_non=TK_LIN1_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      TK_LIN1_MPS_non=TK_LIN1_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      TK_LIN2_MLE_non=TK_LIN2_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      TK_LIN2_MPS_non=TK_LIN2_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      
      
      Lindley_SEL_MLE_non=Lindley_SEL_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      Lindley_SEL_MPS_non=Lindley_SEL_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      Lindley_LIN1_MLE_non=Lindley_LIN1_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      Lindley_LIN1_MPS_non=Lindley_LIN1_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      Lindley_LIN2_MLE_non=Lindley_LIN2_MLE_inf=array(dim=c(3,Sim_no,8),NA)
      Lindley_LIN2_MPS_non=Lindley_LIN2_MPS_inf=array(dim=c(3,Sim_no,8),NA)
      
      # 3rd dim=(Est, AB, MSE, NA, NA, NA,NA,NA)
      start.time=date()
    }  
    
    data = try(simulate_PMITICS_cpp(par_true[1:2],n,m,censor_prop,l_prop,Rprog),silent=T)
    #data1 = try(simulate_MITIICS_cpp(n,par_true[1:2],s,m,odd_len,even_len),silent=T)
    
    if(is.character(data)) {cat("Error...\n"); next}
    
    X = data$X
    L = data$L
    R = data$R
    T = data$T
    final_right_censored =data$final_right_censored
    
    cat("***************************************************************************************************************\n")
    cat("===== iterations ",it," out of ",Sim_no, " dataset =(",index_set, ") scheme set =(",index_scheme,") true parameters=(",par_true[1],",",par_true[2],",",par_true[3],")========\n")
    cat("***************************************************************************************************************\n")
    #TK_cpp(par_true, par,n, X,L,R,T,Rprog,final_right_censored,Cpm_par, h_non, c(c_val1,c_val2), lower_lim, upper_lim, type="LK")
    #################################################################
    #------------- MLE using NR method
    ############s#####################################################
    
    mle.res = suppressWarnings(try(Estim_profile_cpp(par_true,
                                                     par_true[1:2],
                                                     X,L,R,T, Rprog,
                                                     final_right_censored,
                                                     Cpm_par,
                                                     lower_lim,
                                                     upper_lim,
                                                     type = "LK"),
                                   silent=TRUE))
    # Check convergence
    if(is.character(mle.res)) {message("Error...");next}
    if(mle.res$Err==TRUE) next
    if(any(!is.finite(as.numeric(unlist(mle.res))))) next
    
    mle.est = as.numeric(mle.res$Est)
    mle.ese = as.numeric(mle.res$ESE)
    mle.aci = mle.res$ACI
    
    
    ########################################################################
    # Bootstrap method based on LK function
    ########################################################################
    
    cat("Bootstrap sample based on LK funtion...\n")
    
    blk.res = try(Boot_fun_cpp(par_true,
                               init=mle.est,
                               se=mle.ese,
                               n,m,censor_prop,l_prop,Rprog,
                               lower_lim,
                               upper_lim,
                               Cpm_par,
                               BOOT = B,
                               alpha = 0.05,
                               type = "LK"),
                  silent=TRUE)
    
    if(is.character(blk.res)) {message("Error...");next}
    
    
    
    #################################################################
    #------------- MPSE using NR method
    ############s#####################################################
    mps.res = suppressWarnings(try(Estim_profile_cpp(par_true,
                                                     par_true[1:2],
                                                     X,L,R,T, Rprog,
                                                     final_right_censored,
                                                     Cpm_par,
                                                     lower_lim,
                                                     upper_lim,
                                                     type = "PS"),
                                   silent=TRUE))
    #Check convergence
    if(is.character(mps.res)) {message("Error...");next}
    if(mps.res$Err==TRUE) next
    if(any(!is.finite(as.numeric(unlist(mps.res))))) next
    
    mps.est=as.numeric(mps.res$Est)
    mps.ese=as.numeric(mps.res$ESE)
    mps.aci=mps.res$ACI
    
    ########################################################################
    # Bootstrap method based on LK function
    ########################################################################
    cat("Bootstrap sample based on PS funtion...\n")
    
    bps.res = try(Boot_fun_cpp(par_true,
                               init=mps.est,
                               se=mps.ese,
                               n,m,censor_prop,l_prop,Rprog,
                               lower_lim,
                               upper_lim,
                               Cpm_par,
                               BOOT = B,
                               alpha = 0.05,
                               type = "PS"),
                  silent=TRUE)
    
    if(is.character(bps.res)) {message("Error...");next}
    
    MLE[1,it,]= c(mle.est[1],mle.est[1]-par_true[1],(mle.est[1]-par_true[1])^2,mle.ese[1],
                  pmax(0,mle.aci[1,]))
    MLE[2,it,]= c(mle.est[2],mle.est[2]-par_true[2],(mle.est[2]-par_true[2])^2,mle.ese[2],
                  pmax(0,mle.aci[2,]))
    MLE[3,it,]= c(mle.est[3],mle.est[3]-par_true[3],(mle.est[3]-par_true[3])^2,mle.ese[3],
                  pmax(0,mle.aci[3,]))
    
    BOOT_p_MLE[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.bet)
    BOOT_p_MLE[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.gam)
    BOOT_p_MLE[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.Cpm)
    
    BOOT_t_MLE[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.bet)
    BOOT_t_MLE[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.gam)
    BOOT_t_MLE[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.Cpm)
    
    BOOT_log_t_MLE[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.bet)
    BOOT_log_t_MLE[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.gam)
    BOOT_log_t_MLE[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.Cpm)
    
    MPS[1,it,]= c(mps.est[1],mps.est[1]-par_true[1],(mps.est[1]-par_true[1])^2,mps.ese[1],
                  pmax(0,mps.aci[1,]))
    MPS[2,it,]= c(mps.est[2],mps.est[2]-par_true[2],(mps.est[2]-par_true[2])^2,mps.ese[2],
                  pmax(0,mps.aci[2,]))
    MPS[3,it,]= c(mps.est[3],mps.est[3]-par_true[3],(mps.est[3]-par_true[3])^2,mps.ese[3],
                  pmax(0,mps.aci[3,]))
    
    BOOT_p_MPS[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.bet)
    BOOT_p_MPS[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.gam)
    BOOT_p_MPS[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_p.Cpm)
    
    BOOT_t_MPS[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.bet)
    BOOT_t_MPS[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.gam)
    BOOT_t_MPS[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_t.Cpm)
    
    BOOT_log_t_MPS[1,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.bet)
    BOOT_log_t_MPS[2,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.gam)
    BOOT_log_t_MPS[3,it,]= c(NA,NA,NA,NA,blk.res$Boot_log_t.Cpm)
    
    
    ########################################################################
    #------------- Lindley based on LK function and under non-informative prior and SEL
    ############s###########################################################
    lindley_lk_sel_non = try(Lindley_GEL_cpp(True_par=par_true,
                                             Init_par=par,
                                             n=n,
                                             X=X,
                                             L=L,
                                             R=R,
                                             T=T,
                                             Rprog=Rprog,
                                             final_right_censored=final_right_censored,
                                             Cpm_par=Cpm_par,
                                             h=h_non,
                                             q=-1,
                                             type="LK"),silent=TRUE)
    
    
    if(inherits(lindley_lk_sel_non,"try-error")) next
    if(any(!is.finite(c(lindley_lk_sel_non$est_GEL)))) next
    
    ########################################################################
    #------------- Lindley based on LK function and under informative prior and SEL
    ############s###########################################################
    lindley_lk_sel_inf = Lindley_GEL_cpp(True_par=par_true,
                                         Init_par=par,
                                         n=n,
                                         X=X,
                                         L=L,
                                         R=R,
                                         T=T,
                                         Rprog=Rprog,
                                         final_right_censored=final_right_censored,
                                         Cpm_par=Cpm_par,
                                         h=h_inf,
                                         q=-1,
                                         type="LK")
    
    if(inherits(lindley_lk_sel_inf,"try-error")) next
    if(any(!is.finite(c(lindley_lk_sel_inf$est_GEL)))) next
    
    ########################################################################
    #------------- Lindley based on LK function and under non-informative prior and LINEX with c1
    ############s###########################################################
    lindley_lk_linex1_non = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_non,
                                              c=c_val1,
                                              type="LK")
    
    if(inherits(lindley_lk_linex1_non,"try-error")) next
    if(any(!is.finite(c(lindley_lk_linex1_non$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on LK function and under non_informative prior and LINEX with c2
    ############s###########################################################
    lindley_lk_linex2_non = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_non,
                                              c=c_val2,
                                              type="LK")
    
    if(inherits(lindley_lk_linex2_non,"try-error")) next
    if(any(!is.finite(c(lindley_lk_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under informative prior and LINEX with c1
    ############s###########################################################
    lindley_lk_linex1_inf = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_inf,
                                              c=c_val1,
                                              type="LK")
    
    if(inherits(lindley_lk_linex2_non,"try-error")) next
    if(any(!is.finite(c(lindley_lk_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on LK function and under informative prior and LINEX with c2
    ############s###########################################################
    lindley_lk_linex2_inf = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_inf,
                                              c=c_val2,
                                              type="LK")
    
    if(inherits(lindley_lk_linex2_inf,"try-error")) next
    if(any(!is.finite(c(lindley_lk_linex2_inf$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under non-informative prior and SEL
    ############s###########################################################
    lindley_ps_sel_non = Lindley_GEL_cpp(True_par=par_true,
                                         Init_par=par,
                                         n=n,
                                         X=X,
                                         L=L,
                                         R=R,
                                         T=T,
                                         Rprog=Rprog,
                                         final_right_censored=final_right_censored,
                                         Cpm_par=Cpm_par,
                                         h=h_non,
                                         q=-1,
                                         type="PS")
    
    if(inherits(lindley_ps_sel_non,"try-error")) next
    if(any(!is.finite(c(lindley_ps_sel_non$est_GEL)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under informative prior and SEL
    ############s###########################################################
    lindley_ps_sel_inf = Lindley_GEL_cpp(True_par=par_true,
                                         Init_par=par,
                                         n=n,
                                         X=X,
                                         L=L,
                                         R=R,
                                         T=T,
                                         Rprog=Rprog,
                                         final_right_censored=final_right_censored,
                                         Cpm_par=Cpm_par,
                                         h=h_inf,
                                         q=-1,
                                         type="PS")
    
    if(inherits(lindley_ps_sel_inf,"try-error")) next
    if(any(!is.finite(c(lindley_ps_sel_inf$est_GEL)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under non-informative prior and LINEX with c1
    ############s###########################################################
    lindley_ps_linex1_non = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_non,
                                              c=c_val1,
                                              type="PS")
    
    if(inherits(lindley_ps_linex1_non,"try-error")) next
    if(any(!is.finite(c(lindley_ps_linex1_non$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under non_informative prior and LINEX with c2
    ############s###########################################################
    lindley_ps_linex2_non = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_non,
                                              c=c_val2,
                                              type="PS")
    
    if(inherits(lindley_ps_linex2_non,"try-error")) next
    if(any(!is.finite(c(lindley_ps_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under informative prior and LINEX with c1
    ############s###########################################################
    lindley_ps_linex1_inf = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_inf,
                                              c=c_val1,
                                              type="PS")
    
    if(inherits(lindley_ps_linex1_inf,"try-error")) next
    if(any(!is.finite(c(lindley_ps_linex1_inf$est_LINEX)))) next
    
    ########################################################################
    #------------- Lindley based on PS function and under informative prior and LINEX with c2
    ############s###########################################################
    lindley_ps_linex2_inf = Lindley_LINEX_cpp(True_par=par_true,
                                              Init_par=par,
                                              n=n,
                                              X=X,
                                              L=L,
                                              R=R,
                                              T=T,
                                              Rprog=Rprog,
                                              final_right_censored=final_right_censored,
                                              Cpm_par=Cpm_par,
                                              h=h_inf,
                                              c=c_val2,
                                              type="PS")
    
    if(inherits(lindley_ps_linex2_inf,"try-error")) next
    if(any(!is.finite(c(lindley_ps_linex2_inf$est_LINEX)))) next
    
    Lindley_SEL_MLE_non  <- save_bayes(Lindley_SEL_MLE_non,  lindley_lk_sel_non$est_GEL,  par_true, it)
    Lindley_LIN1_MLE_non <- save_bayes(Lindley_LIN1_MLE_non, lindley_lk_linex1_non$est_LINEX, par_true, it)
    Lindley_LIN2_MLE_non <- save_bayes(Lindley_LIN2_MLE_non, lindley_lk_linex2_non$est_LINEX, par_true, it)
    
    Lindley_SEL_MLE_inf  <- save_bayes(Lindley_SEL_MLE_inf,  lindley_lk_sel_inf$est_GEL,  par_true, it)
    Lindley_LIN1_MLE_inf <- save_bayes(Lindley_LIN1_MLE_inf, lindley_lk_linex1_inf$est_LINEX, par_true, it)
    Lindley_LIN2_MLE_inf <- save_bayes(Lindley_LIN2_MLE_inf, lindley_lk_linex2_inf$est_LINEX, par_true, it)
    
    Lindley_SEL_MPS_non  <- save_bayes(Lindley_SEL_MPS_non,  lindley_ps_sel_non$est_GEL,  par_true, it)
    Lindley_LIN1_MPS_non <- save_bayes(Lindley_LIN1_MPS_non, lindley_ps_linex1_non$est_LINEX, par_true, it)
    Lindley_LIN2_MPS_non <- save_bayes(Lindley_LIN2_MPS_non, lindley_ps_linex2_non$est_LINEX, par_true, it)
    
    Lindley_SEL_MPS_inf  <- save_bayes(Lindley_SEL_MPS_inf,  lindley_ps_sel_inf$est_GEL,  par_true, it)
    Lindley_LIN1_MPS_inf <- save_bayes(Lindley_LIN1_MPS_inf, lindley_ps_linex1_inf$est_LINEX, par_true, it)
    Lindley_LIN2_MPS_inf <- save_bayes(Lindley_LIN2_MPS_inf, lindley_ps_linex2_inf$est_LINEX, par_true, it)
    
    ########################################################################
    #------------- TK based on LK function and under non-informative prior and SEL
    ############s###########################################################
    tk_lk_sel_non = try(TK_GEL_cpp(True_par=par_true,
                                   Init_par=par,
                                   n=n,
                                   X=X,
                                   L=L,
                                   R=R,
                                   T=T,
                                   Rprog=Rprog,
                                   final_right_censored=final_right_censored,
                                   Cpm_par=Cpm_par,
                                   h=h_non,
                                   q=-1, 
                                   type="LK"),silent=TRUE)
    
    
    if(inherits(tk_lk_sel_non,"try-error")) next
    if(tk_lk_sel_non$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_sel_non$est_GEL)))) next
    
    ########################################################################
    #------------- TK based on LK function and under informative prior and SEL
    ############s###########################################################
    tk_lk_sel_inf = TK_GEL_cpp(True_par=par_true,
                               Init_par=par,
                               n=n,
                               X=X,
                               L=L,
                               R=R,
                               T=T,
                               Rprog=Rprog,
                               final_right_censored=final_right_censored,
                               Cpm_par=Cpm_par,
                               h=h_inf,
                               q=-1,
                               type="LK")
    
    if(inherits(tk_lk_sel_inf,"try-error")) next
    if(tk_lk_sel_inf$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_sel_inf$est_GEL)))) next
    
    ########################################################################
    #------------- TK based on LK function and under non-informative prior and LINEX with c1
    ############s###########################################################
    tk_lk_linex1_non = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_non,
                                    c=c_val1,
                                    type="LK")
    
    if(inherits(tk_lk_linex1_non,"try-error")) next
    if(tk_lk_linex1_non$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_linex1_non$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on LK function and under non_informative prior and LINEX with c2
    ############s###########################################################
    tk_lk_linex2_non = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_non,
                                    c=c_val2,
                                    type="LK")
    
    if(inherits(tk_lk_linex2_non,"try-error")) next
    if(tk_lk_linex2_non$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on PS function and under informative prior and LINEX with c1
    ############s###########################################################
    tk_lk_linex1_inf = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_inf,
                                    c=c_val1,
                                    type="LK")
    
    if(inherits(tk_lk_linex2_non,"try-error")) next
    if(tk_lk_linex2_non$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on LK function and under informative prior and LINEX with c2
    ############s###########################################################
    tk_lk_linex2_inf = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_inf,
                                    c=c_val2,
                                    type="LK")
    
    if(inherits(tk_lk_linex2_inf,"try-error")) next
    if(tk_lk_linex2_inf$Err==TRUE) next
    if(any(!is.finite(c(tk_lk_linex2_inf$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on PS function and under non-informative prior and SEL
    ############s###########################################################
    tk_ps_sel_non = TK_GEL_cpp(True_par=par_true,
                               Init_par=par,
                               n=n,
                               X=X,
                               L=L,
                               R=R,
                               T=T,
                               Rprog=Rprog,
                               final_right_censored=final_right_censored,
                               Cpm_par=Cpm_par,
                               h=h_non,
                               q=-1,
                               type="PS")
    
    if(inherits(tk_ps_sel_non,"try-error")) next
    if(tk_ps_sel_non$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_sel_non$est_GEL)))) next
    
    ########################################################################
    #------------- TK based on PS function and under informative prior and SEL
    ############s###########################################################
    tk_ps_sel_inf = TK_GEL_cpp(True_par=par_true,
                               Init_par=par,
                               n=n,
                               X=X,
                               L=L,
                               R=R,
                               T=T,
                               Rprog=Rprog,
                               final_right_censored=final_right_censored,
                               Cpm_par=Cpm_par,
                               h=h_inf,
                               q=-1,
                               type="PS")
    
    if(inherits(tk_ps_sel_inf,"try-error")) next
    if(tk_ps_sel_inf$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_sel_inf$est_GEL)))) next
    
    ########################################################################
    #------------- TK based on PS function and under non-informative prior and LINEX with c1
    ############s###########################################################
    tk_ps_linex1_non = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_non,
                                    c=c_val1,
                                    type="PS")
    
    if(inherits(tk_ps_linex1_non,"try-error")) next
    if(tk_ps_linex1_non$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_linex1_non$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on PS function and under non_informative prior and LINEX with c2
    ############s###########################################################
    tk_ps_linex2_non = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_non,
                                    c=c_val2,
                                    type="PS")
    
    if(inherits(tk_ps_linex2_non,"try-error")) next
    if(tk_ps_linex2_non$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_linex2_non$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on PS function and under informative prior and LINEX with c1
    ############s###########################################################
    tk_ps_linex1_inf = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_inf,
                                    c=c_val1,
                                    type="PS")
    
    if(inherits(tk_ps_linex1_inf,"try-error")) next
    if(tk_ps_linex1_inf$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_linex1_inf$est_LINEX)))) next
    
    ########################################################################
    #------------- TK based on PS function and under informative prior and LINEX with c2
    ############s###########################################################
    tk_ps_linex2_inf = TK_LINEX_cpp(True_par=par_true,
                                    Init_par=par,
                                    n=n,
                                    X=X,
                                    L=L,
                                    R=R,
                                    T=T,
                                    Rprog=Rprog,
                                    final_right_censored=final_right_censored,
                                    Cpm_par=Cpm_par,
                                    h=h_inf,
                                    c=c_val2,
                                    type="PS")
    
    if(inherits(tk_ps_linex2_inf,"try-error")) next
    if(tk_ps_linex2_inf$Err==TRUE) next
    if(any(!is.finite(c(tk_ps_linex2_inf$est_LINEX)))) next
    
    TK_SEL_MLE_non  <- save_bayes(TK_SEL_MLE_non,  tk_lk_sel_non$est_GEL,  par_true, it)
    TK_LIN1_MLE_non <- save_bayes(TK_LIN1_MLE_non, tk_lk_linex1_non$est_LINEX, par_true, it)
    TK_LIN2_MLE_non <- save_bayes(TK_LIN2_MLE_non, tk_lk_linex2_non$est_LINEX, par_true, it)
    
    TK_SEL_MLE_inf  <- save_bayes(TK_SEL_MLE_inf,  tk_lk_sel_inf$est_GEL,  par_true, it)
    TK_LIN1_MLE_inf <- save_bayes(TK_LIN1_MLE_inf, tk_lk_linex1_inf$est_LINEX, par_true, it)
    TK_LIN2_MLE_inf <- save_bayes(TK_LIN2_MLE_inf, tk_lk_linex2_inf$est_LINEX, par_true, it)
    
    TK_SEL_MPS_non  <- save_bayes(TK_SEL_MPS_non,  tk_ps_sel_non$est_GEL,  par_true, it)
    TK_LIN1_MPS_non <- save_bayes(TK_LIN1_MPS_non, tk_ps_linex1_non$est_LINEX, par_true, it)
    TK_LIN2_MPS_non <- save_bayes(TK_LIN2_MPS_non, tk_ps_linex2_non$est_LINEX, par_true, it)
    
    TK_SEL_MPS_inf  <- save_bayes(TK_SEL_MPS_inf,  tk_ps_sel_inf$est_GEL,  par_true, it)
    TK_LIN1_MPS_inf <- save_bayes(TK_LIN1_MPS_inf, tk_ps_linex1_inf$est_LINEX, par_true, it)
    TK_LIN2_MPS_inf <- save_bayes(TK_LIN2_MPS_inf, tk_ps_linex2_inf$est_LINEX, par_true, it)
    
    ########################################################################
    #------------- MCMC based on LK function and under non-informative prior
    ############s###########################################################
    cat("Generating MH samples based on likelihood function for non-informative prior...\n");
    
    MH_MLE_non=try(MCMC_cpp(par=par_true, 
                            type = "LK", 
                            init = mle.est[1:2], 
                            X =X, L=L, R=R, T=T, Rprog = Rprog, 
                            final_right_censored=final_right_censored,
                            se = mle.ese[1:2], 
                            h=h_non, 
                            Cpm_par=Cpm_par, 
                            MC_size, MC_burn, 
                            c=c(c_val1,c_val2),
                            lower=lower_lim, 
                            upper=upper_lim),silent=TRUE)
    
    if(inherits(MH_MLE_non,"try-error")) next
    
    if(any(!is.finite(c(MH_MLE_non$est_SEL,
                        MH_MLE_non$est_LIN1,
                        MH_MLE_non$est_LIN2)))) next
    
    
    ########################################################################
    #------------- MCMC based on LK function and under informative prior
    ############s###########################################################
    cat("Generating MH samples based on likelihood function for informative prior...\n");
    
    MH_MLE_inf=try(MCMC_cpp(par=par_true, 
                            type = "LK", 
                            init = mle.est[1:2], 
                            X =X, L=L, R=R, T=T, Rprog = Rprog, 
                            final_right_censored=final_right_censored,
                            se = mle.ese[1:2], 
                            h=h_inf, 
                            Cpm_par=Cpm_par, 
                            MC_size, MC_burn, 
                            c=c(c_val1,c_val2),
                            lower=lower_lim, 
                            upper=upper_lim),silent=TRUE)
    
    if(inherits(MH_MLE_inf,"try-error")) next
    
    if(any(!is.finite(c(MH_MLE_inf$est_SEL,
                        MH_MLE_inf$est_LIN1,
                        MH_MLE_inf$est_LIN2)))) next
    
    ########################################################################
    #------------- MCMC based on PS function and under non-informative prior
    ############s###########################################################
    cat("Generating MH samples based on product spacings function for non-informative prior...\n");
    MH_MPS_non=try(MCMC_cpp(par=par_true, 
                            type = "PS", 
                            init = mps.est[1:2], 
                            X =X, L=L, R=R, T=T, Rprog = Rprog, 
                            final_right_censored=final_right_censored,
                            se = mps.ese[1:2], 
                            h=h_non, 
                            Cpm_par=Cpm_par, 
                            MC_size, MC_burn, 
                            c=c(c_val1,c_val2),
                            lower=lower_lim, 
                            upper=upper_lim),silent=TRUE)
    
    if(inherits(MH_MPS_non,"try-error")) next
    
    if(any(!is.finite(c(MH_MPS_non$est_SEL,
                        MH_MPS_non$est_LIN1,
                        MH_MPS_non$est_LIN2)))) next
    
    ########################################################################
    #------------- MCMC based on PS function and under inf-informative prior
    ############s###########################################################
    cat("Generating MH samples based on product spacings function for inf-informative prior...\n");
    MH_MPS_inf=try(MCMC_cpp(par=par_true, 
                            type = "PS", 
                            init = mps.est[1:2], 
                            X =X, L=L, R=R, T=T, Rprog = Rprog, 
                            final_right_censored=final_right_censored,
                            se = mps.ese[1:2], 
                            h=h_inf, 
                            Cpm_par=Cpm_par, 
                            MC_size, MC_burn, 
                            c=c(c_val1,c_val2),
                            lower=lower_lim, 
                            upper=upper_lim),silent=TRUE)
    
    if(inherits(MH_MPS_inf,"try-error")) next
    
    if(any(!is.finite(c(MH_MPS_inf$est_SEL,
                        MH_MPS_inf$est_LIN1,
                        MH_MPS_inf$est_LIN2)))) next
    
    MCMC_SEL_MLE_non  <- save_bayes(MCMC_SEL_MLE_non,  MH_MLE_non$est_SEL,  par_true, it, MH_MLE_non$HPD)
    MCMC_LIN1_MLE_non <- save_bayes(MCMC_LIN1_MLE_non, MH_MLE_non$est_LIN1, par_true, it)
    MCMC_LIN2_MLE_non <- save_bayes(MCMC_LIN2_MLE_non, MH_MLE_non$est_LIN2, par_true, it)
    
    MCMC_SEL_MLE_inf  <- save_bayes(MCMC_SEL_MLE_inf,  MH_MLE_inf$est_SEL,  par_true, it, MH_MLE_inf$HPD)
    MCMC_LIN1_MLE_inf <- save_bayes(MCMC_LIN1_MLE_inf, MH_MLE_inf$est_LIN1, par_true, it)
    MCMC_LIN2_MLE_inf <- save_bayes(MCMC_LIN2_MLE_inf, MH_MLE_inf$est_LIN2, par_true, it)
    
    MCMC_SEL_MPS_non  <- save_bayes(MCMC_SEL_MPS_non,  MH_MPS_non$est_SEL,  par_true, it, MH_MPS_non$HPD)
    MCMC_LIN1_MPS_non <- save_bayes(MCMC_LIN1_MPS_non, MH_MPS_non$est_LIN1, par_true, it)
    MCMC_LIN2_MPS_non <- save_bayes(MCMC_LIN2_MPS_non, MH_MPS_non$est_LIN2, par_true, it)
    
    MCMC_SEL_MPS_inf  <- save_bayes(MCMC_SEL_MPS_inf,  MH_MPS_inf$est_SEL,  par_true, it, MH_MPS_inf$HPD)
    MCMC_LIN1_MPS_inf <- save_bayes(MCMC_LIN1_MPS_inf, MH_MPS_inf$est_LIN1, par_true, it)
    MCMC_LIN2_MPS_inf <- save_bayes(MCMC_LIN2_MPS_inf, MH_MPS_inf$est_LIN2, par_true, it)
    
    if(it %% 50 == 0) {
      save.image(
        file= paste("CpmLogLogisticPMITICS_Results_index_set_",index_set,"_index_scheme_",index_scheme,".RData",sep="")
      )
    }
    
    if(it==Sim_no){break} else {it=it+1}
  }
  
}



#--- Final results

Results = data.frame(rbind(
apply(MLE, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_p_MLE, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_t_MLE, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_log_t_MLE, c(1,3), mean, na.rm=T)[3,],
apply(MPS, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_p_MPS, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_t_MPS, c(1,3), mean, na.rm=T)[3,],
apply(BOOT_log_t_MPS, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_SEL_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN1_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN2_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_SEL_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN1_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN2_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_SEL_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN1_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN2_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_SEL_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN1_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(Lindley_LIN2_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_SEL_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN1_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN2_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_SEL_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN1_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN2_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_SEL_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN1_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN2_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(TK_SEL_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN1_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(TK_LIN2_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_SEL_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN1_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN2_MLE_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_SEL_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN1_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN2_MLE_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_SEL_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN1_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN2_MPS_inf, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_SEL_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN1_MPS_non, c(1,3), mean, na.rm=T)[3,],
apply(MCMC_LIN2_MPS_non, c(1,3), mean, na.rm=T)[3,]))

colnames(Results)=c("Est","AB","MSE","ESE","LACI","UACI","Len","CP")
rownames(Results) = c("MLE","BOOT_p_MLE","BOOT_t_MLE","BOOT_log_t_MLE",
                      "MPS","BOOT_p_MPS","BOOT_t_MPS","BOOT_log_t_MPS",
                      "Lindley_SEL_MLE_inf","Lindley_LIN1_MLE_inf","Lindley_LIN2_MLE_inf",
                      "Lindley_SEL_MLE_non","Lindley_LIN1_MLE_non","Lindley_LIN2_MLE_non",
                      "Lindley_SEL_MPS_inf","Lindley_LIN1_MPS_inf","Lindley_LIN2_MPS_inf",
                      "Lindley_SEL_MPS_non","Lindley_LIN1_MPS_non","Lindley_LIN2_MPS_non",
                      "TK_SEL_MLE_inf","TK_LIN1_MLE_inf","TK_LIN2_MLE_inf",
                      "TK_SEL_MLE_non","TK_LIN1_MLE_non","TK_LIN2_MLE_non",
                      "TK_SEL_MPS_inf","TK_LIN1_MPS_inf","TK_LIN2_MPS_inf",
                      "TK_SEL_MPS_non","TK_LIN1_MPS_non","TK_LIN2_MPS_non",
                      "MCMC_SEL_MLE_inf","MCMC_LIN1_MLE_inf","MCMC_LIN2_MLE_inf",
                      "MCMC_SEL_MLE_non","MCMC_LIN1_MLE_non","MCMC_LIN2_MLE_non",
                      "MCMC_SEL_MPS_inf","MCMC_LIN1_MPS_inf","MCMC_LIN2_MPS_inf",
                      "MCMC_SEL_MPS_non","MCMC_LIN1_MPS_non","MCMC_LIN2_MPS_non")
print(Results)
