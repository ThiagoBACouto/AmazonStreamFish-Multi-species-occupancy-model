## Model
sink("model_com_occ.txt")     ##Define model
cat("
model {

   ## Priors
# Set uninformative priors for the hyper-parameters  
u.mean1 ~ dunif(0, 1)
muu1 <- log(u.mean1) - log(1 - u.mean1)
u.mean2 ~ dunif(0, 1)
muu2 <- log(u.mean2) - log(1 - u.mean2)
u.mean3 ~ dunif(0, 1)
muu3 <- log(u.mean3) - log(1 - u.mean3)
u.mean4 ~ dunif(0, 1)
muu4 <- log(u.mean4) - log(1 - u.mean4)
v.mean1 ~ dunif(0, 1)                                                  
muv1 <- log(v.mean1) - log(1 - v.mean1)
v.mean2 ~ dunif(0, 1)
muv2 <- log(v.mean2) - log(1 - v.mean2)
v.mean3 ~ dunif(0, 1)
muv3 <- log(v.mean3) - log(1 - v.mean3)

tau.u1 ~ dgamma(0.1, 0.1)                                                            
tau.u2 ~ dgamma(0.1, 0.1)
tau.u3 ~ dgamma(0.1, 0.1)
tau.u4 ~ dgamma(0.1, 0.1)
tau.v1 ~ dgamma(0.1, 0.1)
tau.v2 ~ dgamma(0.1, 0.1)
tau.v3 ~ dgamma(0.1, 0.1)

   ## Define the relationships between species-specific parameters and hyper-parameters
# Loop over species   
for (i in 1: S) {

   u1[i] ~ dnorm(muu1, tau.u1)                                          
   u2[i] ~ dnorm(muu2, tau.u2)
   u3[i] ~ dnorm(muu3, tau.u3)
   u4[i] ~ dnorm(muu4, tau.u4)
   v1[i] ~ dnorm(muv1, tau.v1)
   v2[i] ~ dnorm(muv2, tau.v2)
   v3[i] ~ dnorm(muv3, tau.v3)

   ## Likelihood
   # Loop over sites   
   for (j in 1: J) {

      # Ecological level
  logit(psi[j,i]) <- u1[i] * Chan [j] + u2[i] * Pond[j] + u3[i] * Dist[j] + u4[i] * Chan[j] * Dist[j]   
      z[j,i] ~ dbern(psi[j,i])                                                            
      
      # Loop over visits
      for (k in 1: K[j]){

            # Sampling level
            logit(p[j,k,i]) <- v1[i] * Chan [j] + v2[i] * Pond[j] + v3[i] * Chan[j] * Dist[j]                
            mup[j,k,i] <- p[j,k,i] * z[j,i]                                                    
            Y[j,k,i] ~ dbern(mup[j,k,i])                                                     
         }
      }
  }
  
## Richness estimates - Sum all “true” occupancy state (z) for all species in each site

# Loop over sites
for (j in 1: J){
   rich[j] <- sum(z[j,])
} 

## Similarity estimates between channel and ponds of each segment
# Estimates are based on the “true” occupancy state (z) for all species in each site

# Loop over segments
for (j in 1: W){

# Sorensen similarity index between adjacent sites
    similarity[j] <- 2 * (sum(z[j, ] * z[j + W, ])) / (sum(z[j, ]) + sum(z[j + W, ]))
 }


}
", fill = TRUE)
sink()
