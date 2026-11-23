# Mon Nov 23 12:26:03 2026 ------------------------------

# examine same time windown for 2022 and 2023 observations, 
# 20 July to & 7 August

# Thu May  1 15:24:31 2025 ------------------------------

# bring in ch data and apply Riecke et al. 2022 protocol
# adapted for Atlantic Puffins on Great Island, NL, Canada, 2021-2023 spanning the 2022 HPAIV outbreak
# by G. J. Robertson, Environment and Climate Change Canada
# data collected by Pierre-Paul Bitton and Antoine Morel, Memorial University of Newfoundland

library(jagsUI)
library(rjags)

library(here)
source('cmrfunctions.R')

load(here('data', 'choutsamewindow.Rdata'))

##### set-up data
## scrub out banding year (1st column)
c <- ch[,-1]

n.years <- dim(c)[2]                                       # study duration
releases <- table(ch$byear)                                 # releases per year
f <- rep(seq(1,n.years-1), times = releases)          # first encounter/marking vector
I <- length(f)                                       # total number of individuals


# MCMC specifics
nc <- 3
nt <- 10
ni <- 50000
nb <- 10000


# obtain m - the standard encounter history array by reducing counts to zeros

m <- c
m[c[] >= 1] <- 1

# last seen vector [l], helps with efficiency and establishing known states in
# state matrices a and z
l <- apply(m, 1, get.last)

###
### known alive matrix [z] (provided to JAGS as data)
###
z.dat <- matrix(NA, I, n.years)

for (i in 1:nrow(z.dat)){
  z.dat[i,f[i]] <- 1
  if (l[i] > f[i]){
    z.dat[i, f[i]:l[i]] <- 1
  }
}

###
### known available matrix [a] (provided to JAGS as data)
###
a.dat <- m
a.dat[a.dat == 0] <- NA


###
### marray for CJS
###
marr <- marray(m)
marr
rel <- rowSums(marr)
rel

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1b) CJS latent state
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink(here('JAGS', "cjs_latent.jags"))
cat("
 model {

  phi ~ dbeta(1,1)
  p ~ dbeta(1,1)

  for (i in 1:I){

    for (t in (f[i]+1):n.years){
      z[i,t] ~ dbern(phi * z[i,t-1])
      m[i,t] ~ dbern(z[i,t] * p)
    }
    
  }

}
",fill = TRUE)
sink()

jags.data <- list(n.years = n.years, z = z.dat,
                  I = I, f = f, m = m)

parameters <- c('phi', 'p')

Sys.time()
cjs <- jags(jags.data, NULL, parameters, here("JAGS","cjs_latent.jags"), 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            parallel = T)
Sys.time()  
print(cjs, digits = 3)

# removed robust design script, not needed for this puffin example

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Poisson
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink(here("JAGS", "zip_cmr.jags"))
cat("
model {

  phi ~ dbeta(1,1)
  gamma[1] ~ dbeta(1,1)
  gamma[2] ~ dbeta(1,1)  
  epsilon ~ dgamma(1,1)

  for (i in 1:I){

    for (t in (f[i]+1):n.years){
      z[i,t] ~ dbern(phi * z[i,t-1])
      a[i,t] ~ dbern(z[i,t] * gamma[a[i,t-1]+1])
      c[i,t] ~ dpois(a[i,t] * epsilon)
    }

   }
  
  
}
",fill = TRUE)
sink()

# Data
jags.data <- list(z = z.dat, a = a.dat, n.years = n.years, 
                  I = I, f = f, c = c)

# Initial values: not needed in this case
inits <-  NULL

# Parameters
parameters <- c('phi', 'epsilon', 'gamma')


Sys.time()
zip <- jags(jags.data, inits, parameters, "JAGS/zip_cmr.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            parallel = T)
Sys.time()  
print(zip, digits = 3)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Poisson w/ individual heterogeneity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink(here("JAGS", "ziph_cmr.jags"))
cat("
model {

  phi[1] ~ dbeta(1,1)
  gamma[1] ~ dbeta(1,1)
  epsilon[1] ~ dgamma(1,1)
  theta[1] ~ dunif(0,250)
  
  for (logical in (1:time.theta)) {
    theta[logical + 1] ~ dunif(0,250)
  }

  for (logical in (1:time.phi)) {
    phi[logical + 1] ~ dbeta(1,1)
  }

  for (logical in (1:time.eps)) {
    epsilon[logical + 1] ~ dgamma(1,1)
  }

  for (logical in (1:different.gamma)) {
    gamma[2] ~ dbeta(1,1)
  }



  for (i in 1:I){

    for (t in (f[i]+1):n.years){
    
        h[i, t-1] ~ dgamma(theta[(t-2) * time.theta +1], theta[(t-2) * time.theta +1])

    
      z[i,t] ~ dbern(phi[(t-2) * time.phi + 1] * z[i,t-1])
      a[i,t] ~ dbern(z[i,t] * gamma[(a[i,t-1]) * different.gamma + 1])
      c[i,t] ~ dpois(a[i,t] * epsilon[(t-2) * time.eps + 1] * h[i, t-1])      
    }

   }
  
  
}
",fill = TRUE)
sink()

# set logicals for pooling variables across time
# or remove Markovian process from gamma

time.theta <- FALSE
time.phi <-  FALSE
different.gamma <- TRUE
time.eps <- TRUE


jags.data <- list(z = z.dat, 
                  a = a.dat, 
                  n.years = n.years, 
                  I = I, f = f, c = c,
                  time.theta = as.numeric(time.theta),
                  different.gamma = as.numeric(different.gamma),
                  time.phi = as.numeric(time.phi),
                  time.eps = as.numeric(time.eps)
)
inits <- function(){list(phi = phi, gamma = gamma)}  
inits <- NULL
parameters <- c('phi', 'epsilon', 'gamma', 'theta')

# base model

ziph <- jags(jags.data, inits, parameters, "JAGS/ziph_cmr.jags", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
             parallel = T)

print(ziph, digits = 3)

# non-Markovian transition to available state

time.theta <- FALSE
time.phi <-  FALSE
different.gamma <- FALSE
time.eps <- TRUE

jags.data <- list(z = z.dat, 
                  a = a.dat, 
                  n.years = n.years, 
                  I = I, f = f, c = c,
                  time.theta = as.numeric(time.theta),
                  different.gamma = as.numeric(different.gamma),
                  time.phi = as.numeric(time.phi),
                  time.eps = as.numeric(time.eps)
)

# Initial values: not needed in this case
inits <- NULL

parameters <- c('phi', 'epsilon', 'gamma', 'theta')

ziph1 <- jags(jags.data, inits, parameters, "JAGS/ziph_cmr.jags", 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              parallel = T)

print(ziph1, digits = 3)

# different theta (individual heterogeneity) across years

time.theta <- TRUE
time.phi <-  FALSE
different.gamma <- TRUE
time.eps <- TRUE

jags.data <- list(z = z.dat, 
                  a = a.dat, 
                  n.years = n.years, 
                  I = I, f = f, c = c,
                  time.theta = as.numeric(time.theta),
                  different.gamma = as.numeric(different.gamma),
                  time.phi = as.numeric(time.phi),
                  time.eps = as.numeric(time.eps)
)

# Initial values: not needed in this case
inits <- NULL

parameters <- c('phi', 'epsilon', 'gamma', 'theta')

ziph2 <- jags(jags.data, inits, parameters, "JAGS/ziph_cmr.jags", 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              parallel = T)

print(ziph2, digits = 3)

# different survival across years

time.theta <- TRUE
time.phi <-  TRUE
different.gamma <- TRUE

jags.data <- list(z = z.dat, 
                  a = a.dat, 
                  n.years = n.years, 
                  I = I, f = f, c = c,
                  time.theta = as.numeric(time.theta),
                  different.gamma = as.numeric(different.gamma),
                  time.phi = as.numeric(time.phi),
                  time.eps = as.numeric(time.eps)
)
inits <- function(){list(phi = phi, gamma = gamma)}  
inits <- NULL
parameters <- c('phi', 'epsilon', 'gamma', 'theta')

ziph3 <- jags(jags.data, inits, parameters, "JAGS/ziph_cmr.jags", 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              parallel = T)

print(ziph3, digits = 3)


## bring back no time in eps


time.theta <- TRUE
time.phi <-  TRUE
different.gamma <- TRUE
time.eps <- FALSE

jags.data <- list(z = z.dat, 
                  a = a.dat, 
                  n.years = n.years, 
                  I = I, f = f, c = c,
                  time.theta = as.numeric(time.theta),
                  different.gamma = as.numeric(different.gamma),
                  time.phi = as.numeric(time.phi),
                  time.eps = as.numeric(time.eps)
)
inits <- function(){list(phi = phi, gamma = gamma)}  
inits <- NULL
parameters <- c('phi', 'epsilon', 'gamma', 'theta')

ziph4 <- jags(jags.data, inits, parameters, "JAGS/ziph_cmr.jags", 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
              parallel = T)

print(ziph4, digits = 3)



# summary table

models.to.include <- c("ziph", "ziph1", "ziph2", "ziph3")

out <- NULL

for (i in 1:length(models.to.include)){
  g <- get(models.to.include[i])
  out[[i]] <- c(models.to.include[i], g$DIC, g$mean$deviance)
}

out <- as.data.frame(t(as.data.frame(out)))
rownames(out)<- NULL
colnames(out) <- c('Model', 'DIC', 'deviance')
out[,2:3] <- lapply(out[,2:3], as.numeric)

theta <- c('.', '.', 'time', 'time')
phi <- c('.', '.', '.', 'time')
gamma <- c('Markovian', '.', 'Markovian', 'Markovian')

out <- cbind(out, theta, phi, gamma)

out <- out[,c(1, 4:6, 2:3)]

kableExtra::kbl(out[,-1], digits = 2, format = 'html') |> 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)

# model ziph1 is the best, collapse phi, theta and gamma

# NOTE: did not test for time effects on epsilon, 
# clearly different across two years
# also leads to a clear forcing of annual variation into theta, which is not desirable

round(ziph1$summary[,1:7],3)
kableExtra::kbl(ziph1$summary[1:5,1:7],digits = 3, row.names = TRUE, format = "html") |> 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)


# for comparison
round(zip$summary[,1:7],3)
round(cjs$summary[,1:7],3)

# what does the theta look like
plot(density(rgamma(1000, ziph1$mean$theta, ziph1$mean$theta)))

# something meaningful to report
round(quantile(rgamma(1000, ziph1$mean$theta, ziph1$mean$theta), probs = c(0.025, 0.975)),3)


