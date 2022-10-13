rm(list=ls())
library(jagsUI)

# setwd('C:/Users/stcunnin/Box Sync/Manuscripts/PacificFlyway_IPM-Ibis/revision/fromQing')

# load('data/data ready.RData')
load("data/ipm_data_noBBL.RData")

head(dat[[10]])

# Specify model in JAGS language
sink('gwfg ipm.txt')
cat("
  model {
  #-----------------------------------------
  # 1. Define the priors for the parameters
  #-----------------------------------------
  # 1.1 Population size
  Nobs_tau ~ dgamma(.01, .01)         # Precision for observation error in population size
  Nobs_sd <- 1 / sqrt(Nobs_tau)       # SD for observation error in population size
    
  # 1.2 Survival & recovery
  for (i in 1:nburn) {
    logit_nhs_mu[i] ~ dnorm(0, .01)          # Logit mean non-hunting survival
    nhs_mean[i] <- ilogit(logit_nhs_mu[i])   # Mean non-hunting survival
    logit_nhs_rice[i] ~ dnorm(0, .01)I(-1,1) # Effect of rice on logit non-hunting survival
  } # i
  logit_nhs_ddp ~ dnorm(0, .01)I(-1,1)       # Density dependence on logit non-hunting survival
  logit_nhs_best ~ dnorm(0, .01)I(-1,1)    # Effect of El Nino on logit non-hunting survival
  # logit_nhs_svp ~ dnorm(0, .01)              # Effect of Alaska precipitation on log female age ratio
  # logit_nhs_svt ~ dnorm(0, .01)              # Effect of Alaska temperature on log female age ratio
  logit_nhs_tau ~ dgamma(.01, .01)           # Precision of logit non-hunting survival error
  logit_nhs_sd <- 1 / sqrt(logit_nhs_tau)    # SD of logit non-hunting survival error

  logit_kill_j_mu ~ dnorm(0, .01)           # Logit mean kill rate for juveniles
  kill_j_mean <- ilogit(logit_kill_j_mu)    # Mean kill rate for juveniles
  logit_kill_tau ~ dgamma(.01, .01)         # Precision of logit kill rate
  logit_kill_sd <- 1 / sqrt(logit_kill_tau) # SD of logit kill rate
  kill_ratio ~ dunif(0, 1)                  # Ratio of adult kill rate to juvenile kill rate
  kill_a_mean <- kill_j_mean * kill_ratio   # Mean kill rate for adults
  vul <- 1 / kill_ratio                     # Vulnerability

  lrep_tau ~ dgamma(.01, .01)            # Precision of logit report rate error
  lrep_sd <- 1 / sqrt(lrep_tau)          # SD of logit report rate error
  lrep[1] ~ dnorm(-.8, 100)              # First-year logit report rate
  for (t in 2:nyear) {
    lrep[t] ~ dnorm(lrep[t-1], lrep_tau) # Yearly logit report rate (random walk)
  } # t
  crp ~ dunif(0, 1)                      # Crippling loss rate

  for (t in 1:nyear) {
    # Non-hunting survival for each year
    logit_nhs_pred[t] <- 
      logit_nhs_mu[burn[t]] + 
      logit_nhs_ddp * (log(N[t]) - logy_mean) / logy_sd +
      logit_nhs_rice[burn[t]] * rice[t] + 
      logit_nhs_best * best[t]
      # logit_nhs_svt * sv.temp[t] +
      # logit_nhs_svp * sv.prate[t]
    logit_nhs[t] ~ dnorm(logit_nhs_pred[t], logit_nhs_tau)
    nhs[t] <- ilogit(logit_nhs[t]) # Adult non-hunting survival

    # Kill rate for each year
    logit_kill_j[t] ~ dnorm(logit_kill_j_mu, logit_kill_tau)
    kill_j[t] <- ilogit(logit_kill_j[t])
    kill_a[t] <- kill_j[t] * kill_ratio

    # Annual survival for each year
    sa[t] <- nhs[t] * (1 - kill_a[t])
    sj[t] <- nhs[t] * (1 - kill_j[t])

    # Report rate for each year
    rep[t] <- ilogit(lrep[t])

    # Recovery rate for each year
    rec_a[t] <- kill_a[t] * (1 - crp) * rep[t]
    rec_j[t] <- kill_j[t] * (1 - crp) * rep[t]
  } # t

  # 1.3 Age ratio
  log_ar_mu ~ dnorm(0, .01)    # Log mean female age ratio
  ar_mean <- exp(log_ar_mu) # Mean female age ratio
  log_ar_ddp ~ dnorm(0, .01)        # Density dependence on log female age ratio
  log_ar_prec ~ dnorm(0, .01)       # Effect of Alaska precipitation on log female age ratio
  log_ar_temp ~ dnorm(0, .01)       # Effect of Alaska temperature on log female age ratio
  log_ar_tau ~ dgamma(.01, .01)     # Precision of log female age ratio
  log_ar_sd <- 1 / sqrt(log_ar_tau) # SD of log female age ratio
    
  for (t in 1:nyear) { 
    # Yearly female age ratio
      log_ar_pred[t] <- 
      log_ar_mu + 
      log_ar_ddp * (log(N[t]) - logy_mean) / logy_sd +
      log_ar_prec * prec[t] + 
      log_ar_temp * temp[t]
    log_ar[t] ~ dnorm(log_ar_pred[t], log_ar_tau)  # Log femaleage ratio
    ar[t] <- exp(log_ar[t])                        # Female age ratio
  } # t

  #------------------------------------
  # 2. Derived parameters
  #------------------------------------
  # Population growth rate
  for (t in 1:nyear) {
    lambda[t] <- N[t+1] / N[t]
    logla[t] <- log(lambda[t])
  }

  # Geometric mean for population growth
  mlam <- exp((1/(nyear-1)) * sum(logla[]))
    
  #------------------
  # 3. Process model
  #------------------
  # 3.1 Population size
  N[1] ~ dnorm(90, 1/90)T(0,)
  for (t in 2:(nyear+1)) {
    N_mu[t-1] <- N[t-1] * sa[t-1] + N[t-1] * ar[t-1] * sj[t-1]
    N[t] ~ dnorm(N_mu[t-1], 1/N_mu[t-1])T(0,)
  } # t
 
  # 3.2 M-array
  # 3.2.1 Juveniles
  for (t in 1:nyear) {
    pr.j[t,t] <- (1 - sj[t]) * rec_j[t]
    # Further above main diagonal
    for (j in (t+2):nyear){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):(j-1)]) * (1 - sa[j]) * rec_a[t]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
    } #j
  } #t
    
  for (t in 1:(nyear-1)){
    # One above main diagonal
    pr.j[t,t+1] <- sj[t] * (1 - sa[t+1]) * rec_a[t]
  } #t
    
  # Last column: probability of non-recovery
  for (t in 1:nyear){
    pr.j[t,nyear+1] <- 1 - sum(pr.j[t,1:nyear])
  } #t
    
  # 3.2.2 Adults
  for (t in 1:nyear){
    pr.a[t,t] <- (1 - sa[t]) * rec_a[t]
    # Above main diagonal
    for (j in (t+1):nyear){
      pr.a[t,j] <- prod(sa[t:(j-1)]) * (1 - sa[j]) * rec_a[t]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.a[t,j] <- 0
    } #j
  } #t
    
  # Last column: probability of non-recovery
  for (t in 1:nyear){
    pr.a[t,nyear+1] <- 1-sum(pr.a[t,1:nyear])
  } #t
    
  #----------------------
  # 4. Observation model
  #----------------------
  # 4.1 Population counts
  for (t in 1:(nyear+1)){
    y[t] ~ dnorm(N[t], Nobs_tau)
  } # t
    
  # 4.2 Survival (M-array)
  for (t in 1:nyear){
    marr.j[t,1:(nyear+1)] ~ dmulti(pr.j[t,1:(nyear+1)], band.j[t])
    marr.a[t,1:(nyear+1)] ~ dmulti(pr.a[t,1:(nyear+1)], band.a[t])
  } # t

  # 4.3 Report rate & harvest
  for (t in 1:nyear){
    lrep.mean[t] ~ dnorm(lrep[t], lrep.tau[t])
    harvest[t] ~ dnorm(N[t]*kill_a[t], 1/(N[t]*kill_a[t]*(1-kill_a[t])))
  } # t
    
  # 4.4 Productivity (part collection)
  for (t in 1:nyear) {
    qobs[t] <- (ar[t] * vul) / (1 + ar[t] * vul)
    wing.j[t] ~ dbin(qobs[t], wing.t[t])
  } # t
 
  } # model
  ",fill = TRUE)
sink()

#==========
# Run Jags
#==========
# Bundle data
burn <- dat$burn2
nburn <- length(unique(burn))

jags.data <- list(
  nyear = dat$nyear, 
  y=dat$popcount, 
  logy_mean=mean(log(dat$popcount), na.rm=TRUE), 
  logy_sd=sd(log(dat$popcount), na.rm=TRUE), 
  marr.a=dat$marr.a, marr.j=dat$marr.j, 
  band.a=dat$band.a, band.j=dat$band.j, 
  wing.j=dat$wing.j, wing.t=dat$wing.j + dat$wing.a, 
  rice=(dat$rice - mean(dat$rice)) / sd(dat$rice), 
  burn=burn, nburn=nburn, 
  prec=(dat$prec - mean(dat$prec)) / sd(dat$prec), 
  temp=(dat$temp - mean(dat$temp)) / sd(dat$temp), 
  best=dat$best,
  harvest=dat$harvest, 
  # sv.temp=(dat$sv.temp - mean(dat$sv.temp)) / sd(dat$sv.temp),
  # sv.prate=(dat$sv.prate - mean(dat$sv.prate)) / sd(dat$sv.prate),
  lrep.mean=dat$lrep.mean, lrep.tau=1/(dat$lrep.sd^2))

# Initial values
Ni <- dat$popcount
Ni[which(is.na(Ni))] <- (Ni[which(is.na(Ni))-1] + Ni[which(is.na(Ni))+1]) / 2
Ni <- round(Ni)
inits <- function() {
  list(
    N=Ni, Nobs_tau=1, 
    logit_nhs_mu=rep(0,nburn), logit_nhs_rice=rep(0,nburn), 
    logit_nhs_ddp=0, logit_nhs_tau=1, logit_nhs_best=0, 
    # logit_nhs_svt=0, logit_nhs_svp=0,
    logit_nhs=rep(0, dat$nyear), 
    logit_kill_j_mu=0, logit_kill_tau=1, 
    logit_kill_j=rep(0, dat$nyear), 
    kill_ratio=.5, crp=.2, 
    lrep_tau=1, lrep=dat$lrep.mean, 
    log_ar_mu=0, #log_ar_rice=rep(0,nburn), 
    log_ar_ddp=0, log_ar_prec=0, log_ar_temp=0, log_ar_tau=1, 
    log_ar=rep(0, dat$nyear)
  )}

# Parameters monitored
parameters <- c(
  'N', 'lambda', 'Nobs_sd', 'mlam',
  'nhs_mean', 'logit_nhs_ddp', 'logit_nhs_rice', 'logit_nhs_best', 
  # 'logit_nhs_svt', 'logit_nhs_svp',
  'logit_nhs_sd', 'nhs', 'sa', 'sj', 
  'ar_mean', 'log_ar_ddp', 'log_ar_rice', 'log_ar_prec', 'log_ar_temp', 
  'log_ar_sd', 'ar', 
  'kill_a_mean', 'kill_j_mean', 'logit_kill_sd', 'kill_a', 'kill_j', 
  'vul', 'crp', 
  'lrep_sd', 'rep'
)

# Call JAGS from R
fit <- jags(jags.data, inits, parameters, model.file='gwfg ipm.txt', 
            #            n.chains=1, n.adapt=100, n.burnin=100, n.iter=200, n.thin=1, 
            # n.chains=3, n.adapt=2000, n.burnin=16000, n.iter=20000, n.thin=1, 
            n.chains=3, n.adapt=2000, n.burnin=240000, n.iter=300000, n.thin=1, 
            parallel=TRUE)

print(fit, digits=3)

save(fit, file='data/ipmfit_noBBL_noCSE.RData')
# save(fit, file="data/ipm_YKDBBL.RData")

fit.sum <- fit$summary
write.csv(fit.sum, "output/ipmfit_noBBL_noCSE.csv")

