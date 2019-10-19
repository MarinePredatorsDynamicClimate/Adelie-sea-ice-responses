rm(list = ls()) #Clear working memory

# Required packages
my.packs <- c(
  
  #Accessing data
  'RPostgreSQL', 
  
  #Data manipulation
  'dplyr', 'reshape2', 'raster',
  
  #Data analysis
  'jagsUI',
  
  #Plotting
  'ggplot2', 'viridis', 'ggrepel',"RColorBrewer")

# if any of them are not installed, install them
if (any(!my.packs %in% installed.packages()[, 'Package'])) {install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my.packs, require, character.only = TRUE)

####################################################################
# Study parameters
####################################################################
# Minimum number of years of counts per site to include in analysis
min_obs = 10

# Range of years to consider in analysis
min_year = 1984
max_year = 2017

##############################################################################################################
# LOAD AND ANALYZE DATA
##############################################################################################################

load(file = paste0("../generated_files/jags_data_multilevel_SIC_",min_obs,"_yrs.rda"))

years = as.numeric(colnames(jags.data$nest_train))

sink("../generated_files/ADPE_model_multilevel_SIC.jags")
cat("
    
    model {
    
    #------------------------------------------------------
    # Priors
    #------------------------------------------------------
    
    # Global quadratic effect (test of Smith et al. 1999)
    beta0 ~ dnorm(0,0.01)
    beta1 ~ dnorm(0,0.01)
    beta2 ~ dnorm(0,0.01)
    
    # Deviations from global quadratic effect
    r.global.sd ~ dunif(0,2)
    r.global.tau <- pow(r.global.sd, -2)
    
    # Distribution of process variances (helps constrain estimates for some colonies)
    v.hyper.mean ~ dunif(0,2)
    log.v.hyper.mean <- log(v.hyper.mean)
    log.v.hyper.sd ~ dunif(0,2)
    log.v.hyper.tau <- pow(log.v.hyper.sd, -2)
    
    alpha1.intercept ~ dnorm(0,0.01)
    alpha1.slope ~ dnorm(0,0.01)
    alpha1.sd ~ dunif(0,5)
    alpha1.tau <- pow(alpha1.sd,-2)

    alpha5.intercept ~ dnorm(0,0.01)
    alpha5.slope ~ dnorm(0,0.01)
    alpha5.sd ~ dunif(0,5)
    alpha5.tau <- pow(alpha5.sd,-2)

    alpha6.intercept ~ dnorm(0,0.01)
    alpha6.slope ~ dnorm(0,0.01)
    alpha6.sd ~ dunif(0,5)
    alpha6.tau <- pow(alpha6.sd,-2)

    for (s in 1:n_sites){

      alpha1.deviate[s] ~ dnorm(0,alpha1.tau)
      alpha5.deviate[s] ~ dnorm(0,alpha5.tau)
      alpha6.deviate[s] ~ dnorm(0,alpha6.tau)

      alpha1[s] <- alpha1.intercept + alpha1.slope * SIC_sitemean[s] + alpha1.deviate[s]
      alpha5[s] <- alpha5.intercept + alpha5.slope * SIC_sitemean[s] + alpha5.deviate[s]
      alpha6[s] <- alpha6.intercept + alpha6.slope * SIC_sitemean[s] + alpha6.deviate[s]

      # Temporal process variance (different at each site)
      v[s] ~ dlnorm(log.v.hyper.mean,log.v.hyper.tau)
      v.tau[s] <- pow(v[s],-2)
      
      # Priors for initial abundance at each site
      N_init[s] ~ dunif(1,max_count[s]*10)
      N[s,first_obs[s]] <- N_init[s]
      logN[s,first_obs[s]] <- log(N[s,first_obs[s]])
    }
    
    #------------------------------------------------------
    # Process Model
    #------------------------------------------------------
    
    for (s in 1:n_sites){
      
      # Mean growth rate at the site (quadratic function of mean SIC)
      r.mu.determ[s] <-  beta0 + beta1 * SIC_sitemean[s] + beta2 * pow(SIC_sitemean[s],2)
      r.mu.deviate[s] ~ dnorm(0, r.global.tau) # Deviations from the global quadratic effect
      r.mu[s] <- r.mu.determ[s] + r.mu.deviate[s]
  
      # Describe temporal dynamics of abundance at the site
      for (t in first_obs[s]:(last_obs[s]-1)) {
        
        r.determ[s,t] <-  r.mu[s] + alpha1[s] * SIC_anomaly_1yr[s,t] + alpha5[s] * SIC_anomaly_5yr[s,t] + alpha6[s] * SIC_anomaly_6yr[s,t]
        r[s,t] ~ dnorm(r.determ[s,t], v.tau[s]) #Add temporal process variance
        r.deviate[s,t] <- r[s,t] - r.determ[s,t]

        N[s,t+1] <- exp(logN[s,t+1])
        logN[s,t+1] <- logN[s,t] + r[s,t]
        
      }
    }
    
    # Convert chick counts to adult counts
    mean.p ~ dunif(0,1)
    logit.p <- log(mean.p/(1-mean.p))
    temporal.sigma.reprod ~ dunif(0,5) # Temporal variance in annual reproductive success
    temporal.precision.reprod <- pow(temporal.sigma.reprod,-2)
    
    for (s in 1:n_sites){
      for (t in 1:n_seasons) {
        lp[s,t] ~ dnorm(logit.p, temporal.precision.reprod)
        p[s,t] <- 1/(1+exp(-lp[s,t]))
      }
    }
    
    #------------------------------------------------------
    # Observation model
    #------------------------------------------------------
    for (i in 1:n_length) {
      
      #Observed counts 
      nest_train[site[i],season[i]] ~ dlnorm(log(N[site[i],season[i]]), nest_precision[site[i],season[i]])
      chick_train[site[i],season[i]] ~ dlnorm(log(2*p[site[i],season[i]]*N[site[i],season[i]]), chick_precision[site[i],season[i]])
      adult_train[site[i],season[i]] ~ dlnorm(log(N[site[i],season[i]]), adult_precision[site[i],season[i]])
      
    }
    
    #------------------------------------------------------
    # Posterior predictive check
    #------------------------------------------------------
    # for (i in 1:n_length) {
    # 
    #   #Simulate data under fitted model
    #   nest_ynew[site[i], season[i]] ~ dlnorm(log(N[site[i],season[i]]), nest_precision[site[i],season[i]])
    #   chick_ynew[site[i], season[i]] ~ dlnorm(log(2*p[site[i],season[i]]*N[site[i],season[i]]), chick_precision[site[i],season[i]])
    #   adult_ynew[site[i], season[i]] ~ dlnorm(log(N[site[i],season[i]]), adult_precision[site[i],season[i]])
    # 
    #   #Compare observed data to fitted model
    #   nest_ysqActual[i]  <- pow((nest_train[site[i], season[i]] - N[site[i],season[i]]), 2)
    #   chick_ysqActual[i]  <- pow((chick_train[site[i], season[i]] - 2*p[site[i],season[i]]*N[site[i],season[i]]), 2)
    #   adult_ysqActual[i]  <- pow((adult_train[site[i], season[i]] - N[site[i],season[i]]), 2)
    # 
    #   #Compare simulated data to fitted model
    #   nest_ysqNew[i]  <- pow((nest_ynew[site[i], season[i]] - N[site[i],season[i]]), 2)
    #   chick_ysqNew[i]  <- pow((chick_ynew[site[i], season[i]] - 2*p[site[i],season[i]]*N[site[i],season[i]]), 2)
    #   adult_ysqNew[i]  <- pow((adult_ynew[site[i], season[i]] - N[site[i],season[i]]), 2)
    # 
    # }
    # 
    # nest_zActual <- sum(nest_ysqActual[])
    # nest_zNew <- sum(nest_ysqNew[])
    # 
    # chick_zActual <- sum(chick_ysqActual[])
    # chick_zNew <- sum(chick_ysqNew[])
    # 
    # adult_zActual <- sum(adult_ysqActual[])
    # adult_zNew <- sum(adult_ysqNew[])

    
    }
    
    ",fill = TRUE)
sink()

# Generate some initial values for abundance
N_init = array(NA,dim=c(jags.data$n_sites,1))
jags.data$max_count = rep(1000000,jags.data$n_sites)

#Add placeholder for Terre Adelie (this is not necessary when the actual data is included)
jags.data$first_obs[38] = 1
jags.data$last_obs[38] = 2

for (i in 1:jags.data$n_sites){
  fo = jags.data$first_obs[i] #Year of first observation
  m1 = mean(jags.data$nest_train[i,],na.rm=TRUE); if (is.nan(m1)) m1 = 0
  m2 = mean(jags.data$adult_train[i,],na.rm=TRUE); if (is.nan(m2)) m2 = 0
  N_init[i,1] = (m1 + m2)/2
}

inits <- function(){list(N_init = N_init[,1],beta0 = 0,beta1 = 0,beta2 = 0, r.global.sd  = 0.01,mean.p = 0.5, temporal.sigma.reprod = 0.2)}

start = Sys.time()

# Fit model
out <- jags(data=jags.data,
            model.file="../generated_files/ADPE_model_multilevel_SIC.jags",
            parameters.to.save=c(
              
              #Hyper-parameters describing effect of long-term SIC on long-term population growth rates
              "beta0",
              "beta1",
              "beta2",
              "r.global.sd",
              
              #Hyper-parameters describing effect of annual SIC covariates on annual growth rates
              "alpha1.intercept",
              "alpha1.slope",
              "alpha1.sd",
              
              "alpha5.intercept",
              "alpha5.slope",
              "alpha5.sd",
              
              "alpha6.intercept",
              "alpha6.slope",
              "alpha6.sd",
              
              #Hyper-parameters describing process variances
              "v.hyper.mean",
              "log.v.hyper.sd",
              
              # Long-term population growth rates
              "r.mu",
              "r.mu.determ",
              "r.mu.deviate",
              
              # Annual population growth rates
              "r",
              "r.determ",
              "r.deviate",
              
              # Regression coefficients at each site (responses to annual SIC)
              "alpha1",
              "alpha5",
              "alpha6",
              
              # Reproductive success
              "mean.p",
              "temporal.sigma.reprod",
              
              #Estimated nest abundance
              "N",
              
              #Process variances
              "v",
              
              #Posterior predictive checks
              "nest_zActual",
              "nest_zNew",
              "chick_zActual",
              "chick_zNew",
              "adult_zActual",
              "adult_zNew"
            ),
            
            inits = inits,
            n.chains=3,
            n.thin = 150,
            n.iter=1750000,
            n.burnin= 250000)

stop = Sys.time()

out

print(stop-start) #Time taken to fit the model
