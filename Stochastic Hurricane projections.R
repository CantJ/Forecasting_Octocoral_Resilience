# This script is for running a series of simulations exploring the impacts of changing recurrent hurricane regimes on the dynamics and viability of
# 3 different gorgonia populations: Gorgonia, Eunicea and Antillogorgia.
# 
# data last modified: Mar 2022
# Author: James Cant

# Load required packages
library(areaplot)

###############################################
# STEP 1: Calculate final model projection parameters
###############################################
# reset plot window
set_graph_pars("panel1")

# determine bin widths for integrating the population
h.store = (U.list-L.list)/m # bin widths

# define simulation storage
# all of the storage vectors are ordered as 1 = Antillogorgia, 2 = Eunicea, 3 = Gorgonia.
dates <- 2018:2100 # Projection dates
tmax = length(dates) # iteration length
sim <- 1000 # number of simulations
# set up one storage matrix for each species under each simulation
ant.nt <- flex.nt <- gorg.nt <- matrix(NA, m, tmax)
ant.nt2 <- flex.nt2 <- gorg.nt2 <- matrix(NA, m, tmax)
ant.nt3 <- flex.nt3 <- gorg.nt3 <- matrix(NA, m, tmax)

# total population stores
ant.Nt <- gorg.Nt <- flex.Nt <- numeric(tmax)
ant.Nt2 <- gorg.Nt2 <- flex.Nt2 <- numeric(tmax)
ant.Nt3 <- gorg.Nt3 <- flex.Nt3 <- numeric(tmax)

# Add initial population state for each projection
n0 <- rep(1, m)
N0 <- sum(n0)
ant.nt[,1] <- ant.nt2[,1] <- ant.nt3[,1] <- flex.nt[,1] <- flex.nt2[,1] <- flex.nt3[,1] <- gorg.nt[,1] <- gorg.nt2[,1] <- flex.nt3[,1] <- n0
ant.Nt[1] <- ant.Nt2[1] <- ant.Nt3[1] <- flex.Nt[1] <- flex.Nt2[1] <- flex.Nt3[1] <- gorg.Nt[1] <- gorg.Nt2[1] <- flex.Nt3[1] <- N0

# Because recruitment also relies on total colony density this needs saving.
total.Nt <- total.Nt2 <- total.Nt3 <- numeric(tmax)
total.Nt[1] <- total.Nt2[1] <- total.Nt3[1] <- N0*3

# Now ready for the moment of truth.
###############################################
# STEP 2: Run the simulations
###############################################
# iterate t = 1 to tmax, one column per year - after each iteration the loop needs to stop to calculate the current populations state.

# Assign hurricane frequency
# Scenario 1: No Hurricanes
hurricane.prob <- rep(0,length.out = tmax)
# Scenario 2:
# Hurricanes Maria and Irma passed the US Virgin Islands as Category 2 Hurricanes with wind speeds of 106Mph. Wind speeds of this magnitude have estimated return times of 1.54 years (0.65 annual probability)
# However Knutson (2019) show that climate change in the North Atlantic are expected to reduce Hurricane frequencies by 20% (0.65*0.2 = 0.13). So by the year 2055 under a 2 degree climate scenario hurricane frequencies will be 0.52.
# Create a probability vector representing this.
mid.century <- seq(0.65, 0.52,length.out = length(2018:2055)) #this reflects a 20% decline by 2055
# but what annual change  does this equate to?
ann.change <- (max(mid.century)-min(mid.century))/length(mid.century)
# now ensure the probability of hurricanes changes annually at the same rate until the end of the century.
hurricane.prob2 <- c(mid.century, seq(min(mid.century)-ann.change, min(mid.century)-(ann.change*(tmax-length(mid.century))), length.out = tmax-length(mid.century))) # hurricane declines will continue at the same rate until the end of the century
Recovery.years2 <- matrix(0, sim, tmax) # define storage for recording if a hurricane occurred in the previous year
# Scenario 3: Increasing Hurricane occurrence
seq(from = 0.65, by = 0.01, length.out = tmax)

hurricane.prob3 <- rep(0,length.out = tmax)
Recovery.years3 <- matrix(0, sim, tmax)

# create overall storage dataframe
ant.store <- matrix(NA, sim, tmax)
flex.store <- matrix(NA, sim, tmax)
gorg.store <- matrix(NA, sim, tmax)
ant.store2 <- matrix(NA, sim, tmax)
flex.store2 <- matrix(NA, sim, tmax)
gorg.store2 <- matrix(NA, sim, tmax)
colnames(ant.store) <- colnames(flex.store) <- colnames(gorg.store) <- colnames(ant.store2) <- colnames(flex.store2) <- colnames(gorg.store2) <- dates
# somewhere to store lambda s
lambdas.store <- matrix(NA, nrow = sim, ncol = 6) #two columns for each species (one for sim1 the other for sim 2)
colnames(lambdas.store) <- c("Antillogorgia.1", "Eunicea.1", "Gorgonia.1", "Antillogorgia.2", "Eunicea.2", "Gorgonia.2")

#--------------------------------------------------------------------------------------------------------------------------------
# run the grand loop for simulation 1
for(s in 1:sim) {
  # progress read out
  print(s)
  
  # this will run a load of different simulations so I can take a mean
  for(t in 2:tmax) {
    # extract associated hurricane potential
    Hurricane <- hurricane.prob2[t] 
    
    # force recovery dynamics in first run through
    if(t == 2){ # the model will first subject the population to recovery dynamics (the first annual step is 2018/2019 and so was the recovery from Hurricane Maria)
      # Antillogorgia
      P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
      ant.nt[,t] <- Iterate.a(nt = ant.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
      ant.Nt[t] <- h.store[1] * sum(ant.nt[,t])
      
      # Eunicea
      P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
      flex.nt[,t] <- Iterate.e(nt = flex.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
      flex.Nt[t] <- h.store[2] * sum(flex.nt[,t])
      
      # Gorgonia
      P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
      gorg.nt[,t] <- Iterate.g(nt = gorg.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
      gorg.Nt[t] <- h.store[3] * sum(gorg.nt[,t])
    } else {
      # define if a hurricane is happening
      prob <- runif(1, min = 0, max = 1)
      if (prob < Hurricane){ # when the model is in a hurricane year - use the parameters measured during 2017/18
        # the model will run through each species separately
        # Antillogorgia
        P.a <- mk_P(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
        ant.nt[,t] <- Iterate.a(nt = ant.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,5]) 
        ant.Nt[t] <- h.store[1] * sum(ant.nt[,t-1])
        
        # Eunicea
        P.e <- mk_P(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
        flex.nt[,t] <- Iterate.e(nt = flex.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,4]) 
        flex.Nt[t] <- h.store[2] * sum(flex.nt[,t])
        
        # Gorgonia
        P.g <- mk_P(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
        gorg.nt[,t] <- Iterate.g(nt = gorg.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,5]) 
        gorg.Nt[t] <- h.store[3] * sum(gorg.nt[,t])
        
        # record the hurricane
        Recovery.years[s, t] <- 1
      } else {
        if (Recovery.years[s, t-1] == 1){ #the year following a hurricane is recovery dynamics
          # Antillogorgia
          P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
          ant.nt[,t] <- Iterate.a(nt = ant.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
          ant.Nt[t] <- h.store[1] * sum(ant.nt[,t])
          
          # Eunicea
          P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
          flex.nt[,t] <- Iterate.e(nt = flex.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
          flex.Nt[t] <- h.store[2] * sum(flex.nt[,t])
          
          # Gorgonia
          P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
          gorg.nt[,t] <- Iterate.g(nt = gorg.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
          gorg.Nt[t] <- h.store[3] * sum(gorg.nt[,t])
        } else { #if its not a hurricane year or a recovery year then the loop will either randomly select transitions from the "good" years or if its a recovery year then select parameters from 2018/2019
          m.par.ant <- ant.params[,sample(1:4,1)] #randomly selects a column from the good year columns of the m.par matrices
          m.par.flex <- flex.params[,sample(1:3,1)]
          m.par.gorg <- gorg.params[,sample(1:4,1)]
          
          # Antillogorgia
          P.a <- mk_P(m = m, m.par = m.par.ant, L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
          ant.nt[,t] <- Iterate.a(nt = ant.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = m.par.ant) 
          ant.Nt[t] <- h.store[1] * sum(ant.nt[,t])
          
          # Eunicea
          P.e <- mk_P(m = m, m.par = m.par.flex, L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
          flex.nt[,t] <- Iterate.e(nt = flex.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = m.par.flex) 
          flex.Nt[t] <- h.store[2] * sum(flex.nt[,t])
          
          # Gorgonia
          P.g <- mk_P(m = m, m.par = m.par.gorg, L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
          gorg.nt[,t] <- Iterate.g(nt = gorg.nt[,t-1], Nt = total.Nt[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = m.par.gorg) 
          gorg.Nt[t] <- h.store[3] * sum(gorg.nt[,t])
        }}}
    
    # at the end of each loop save the total populations size
    total.Nt[t] <- ant.Nt[t] + flex.Nt[t] + gorg.Nt[t]
  }
  # Store the first simulation run through
  ant.store[s,] <- ant.Nt
  flex.store[s,] <- flex.Nt
  gorg.store[s,] <- gorg.Nt
  # and calculate lambda s for each simulation
  # create blank storage  
  lambdas.ant <- numeric(tmax-1) 
  lambdas.flex <- numeric(tmax-1) 
  lambdas.gorg <- numeric(tmax-1) 
  # determine stepwise lambda
  for (x in 1: tmax-1) {
    lambdas.ant[x] <- ant.Nt[x+1]/ant.Nt[x]
    lambdas.flex[x] <- flex.Nt[x+1]/flex.Nt[x]
    lambdas.gorg[x] <- gorg.Nt[x+1]/gorg.Nt[x]
    # and store overall lambda s for the simulation
    lambdas.store[s,1] <- exp(mean(log(lambdas.ant), na.rm = T))
    lambdas.store[s,2] <- exp(mean(log(lambdas.flex), na.rm = T))
    lambdas.store[s,3] <- exp(mean(log(lambdas.gorg), na.rm = T)) 
  }
}

## -----------------------------------------------------------------------------------------------------------------------------
## and re-run the grand loop for simulation 2
#for(s in 1:sim) {
#  # progress read out
#  print(s)
#
#  # this will run a load of different simulations so I can take a mean
#  for(t in 2:tmax) {
#    # extract associated hurricane potential
#    Hurricane <- hurricane.prob2[t] 
#    
#    # force recovery dynamics in first run through
#    if(t == 2){ # the model will first subject the population to recovery dynamics (the first annual step is 2018/2019 and so was the recovery from Hurricane Maria)
#      # Antillogorgia
#      P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
#      ant.nt2[,t] <- Iterate.a(nt = ant.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
#      ant.Nt2[t] <- h.store[1] * sum(ant.nt2[,t])
#      
#      # Eunicea
#      P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
#      flex.nt2[,t] <- Iterate.e(nt = flex.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
#      flex.Nt2[t] <- h.store[2] * sum(flex.nt2[,t])
#      
#      # Gorgonia
#      P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
#      gorg.nt2[,t] <- Iterate.g(nt = gorg.nt2[,t-1], Nt = gorg.Nt2[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
#      gorg.Nt2[t] <- h.store[3] * sum(gorg.nt2[,t])
#    } else {
#      # define if a hurricane is happening
#      prob <- runif(1, min = 0, max = 1)
#      if (prob < Hurricane){ # when the model is in a hurricane year - use the parameters measured during 2017/18
#        # the model will run through each species seperately
#        # Antillogorgia
#        P.a <- mk_P(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
#        ant.nt2[,t] <- Iterate.a(nt = ant.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,5]) 
#        ant.Nt2[t] <- h.store[1] * sum(ant.nt2[,t])
#        
#        # Eunicea
#        P.e <- mk_P(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
#        flex.nt2[,t] <- Iterate.e(nt = flex.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,4]) 
#        flex.Nt2[t] <- h.store[2] * sum(flex.nt2[,t])
#        
#        # Gorgonia
#        P.g <- mk_P(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
#        gorg.nt2[,t] <- Iterate.g(nt = gorg.nt2[,t-1], Nt = gorg.Nt2[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,5]) 
#        gorg.Nt2[t] <- h.store[3] * sum(gorg.nt2[,t])
#        
#        # record the hurricane
#        Recovery.years2[s, t] <- 1
#      } else {
#        if (Recovery.years2[s, t-1] == 1){ #the year following a hurricane is recovery dynamics
#          # Antillogorgia
#          P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
#          ant.nt2[,t] <- Iterate.a(nt = ant.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
#          ant.Nt2[t] <- h.store[1] * sum(ant.nt2[,t])
#          
#          # Eunicea
#          P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
#          flex.nt2[,t] <- Iterate.e(nt = flex.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
#          flex.Nt2[t] <- h.store[2] * sum(flex.nt2[,t])
#          
#          # Gorgonia
#          P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
#          gorg.nt2[,t] <- Iterate.g(nt = gorg.nt2[,t-1], Nt = gorg.Nt2[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
#          gorg.Nt2[t] <- h.store[3] * sum(gorg.nt2[,t])
#        } else { #if its not a hurricane year or a recovery year then the loop will either randomly select transitions from the "good" years or if its a recovery year then select parameters from 2018/2019
#          m.par.ant <- ant.params[,sample(1:4,1)] #randomly slects a column from the good year columns of the m.par matrices
#          m.par.flex <- flex.params[,sample(1:3,1)]
#          m.par.gorg <- gorg.params[,sample(1:4,1)]
#          
#          # Antillogorgia
#          P.a <- mk_P(m = m, m.par = m.par.ant, L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
#          ant.nt2[,t] <- Iterate.a(nt = ant.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = m.par.ant) 
#          ant.Nt2[t] <- h.store[1] * sum(ant.nt2[,t])
#          
#          # Eunicea
#          P.e <- mk_P(m = m, m.par = m.par.flex, L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
#          flex.nt2[,t] <- Iterate.e(nt = flex.nt2[,t-1], Nt = total.Nt2[t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = m.par.flex) 
#          flex.Nt2[t] <- h.store[2] * sum(flex.nt2[,t])
#          
#          # Gorgonia
#          P.g <- mk_P(m = m, m.par = m.par.gorg, L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
#          gorg.nt2[,t] <- Iterate.g(nt = gorg.nt2[,t-1], Nt = gorg.Nt2[t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = m.par.gorg) 
#          gorg.Nt2[t] <- h.store[3] * sum(gorg.nt2[,t])
#        }}}
#    
#    # at the end of each loop save the total populations size
#    total.Nt2[t] <- ant.Nt2[t] + flex.Nt2[t] + gorg.Nt2[t]
#  }
#  # Store the first simulation run through
#  ant.store2[s,] <- ant.Nt2
#  flex.store2[s,] <- flex.Nt2
#  gorg.store2[s,] <- gorg.Nt2
#  # and calculate lambda s for each simulation 
#  # create blank storage  
#  lambdas.ant2 <- numeric(tmax-1) 
#  lambdas.flex2 <- numeric(tmax-1) 
#  lambdas.gorg2 <- numeric(tmax-1) 
#  # determine stepwise lambda
#  for (x in 1: tmax-1) {
#    lambdas.ant2[x] <- ant.Nt2[x+1]/ant.Nt2[x]
#    lambdas.flex2[x] <- flex.Nt2[x+1]/flex.Nt2[x]
#    lambdas.gorg2[x] <- gorg.Nt2[x+1]/gorg.Nt2[x]
#    # and store overall lambda s for the simulation
#    lambdas.store[s,4] <- exp(mean(log(lambdas.ant2), na.rm = T))
#    lambdas.store[s,5] <- exp(mean(log(lambdas.flex2), na.rm = T))
#    lambdas.store[s,6] <- exp(mean(log(lambdas.gorg2), na.rm = T)) 
#  }
#}

###############################################
# STEP 3: Plot simulations and variance and statistical differences
###############################################

# Plot the simulations -----------------------------------------------------------------------------------
#what is my max y value
max(apply(flex.store,2, CI)[1,])
#max(apply(flex.store2,2, CI)[1,])

# plot simulation one and two together
# Eunicea RPC2.6
plot(x = dates, apply(flex.store,2, CI)[2,], typ = "n", ylim = c(0,300), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
#set a transparent colour
bluetrans <- rgb(0,114,178, alpha = 100, maxColorValue = 255)
confplot(x = dates, y1 = apply(flex.store,2, CI)[3,], # lower CI values
         y2 = apply(flex.store,2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store,2, CI)[2,], typ = "l", lty = "dashed", col = "#0072B2", add = TRUE, lwd = 1)
# Eunicea RPC8.5
#confplot(x = dates, y1 = apply(flex.store2,2, CI)[3,], # lower CI values
# y2 = apply(flex.store2,2, CI)[1,], #upper CI values
# col = bluetrans, # polygon colour
# add = TRUE)
# add the sample mean
#points(x = dates, apply(flex.store2,2, CI)[2,], typ = "l", lty = "dotted", col = "blue", add = TRUE, lwd = 1)

# Antilogorgia RCP2.6
#set a transparent colour
redtrans <- rgb(213,94,0, alpha = 100, maxColorValue = 255)
confplot(x = dates, y1 = apply(ant.store,2, CI)[3,], # lower CI values
         y2 = apply(ant.store,2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store,2, CI)[2,], typ = "l", lty = "dashed", col = "#D55E00", add = TRUE, lwd = 1)
# Antilogorgia RCP8.5
#confplot(x = dates, y1 = apply(ant.store2,2, CI)[3,], # lower CI values
#         y2 = apply(ant.store2,2, CI)[1,], #upper CI values
#         col = goldtrans, # polygon colour
#         add = TRUE)
## add the sample mean
#points(x = dates, apply(ant.store2,2, CI)[2,], typ = "l", lty = "dotted", col = "goldenrod4", add = TRUE, lwd = 1)

# Gorgonia RCP2.6
#set a transparent colour
blacktrans <- rgb(0,0,0, alpha = 150, maxColorValue = 255)
confplot(x = dates, y1 = apply(gorg.store,2, CI)[3,], # lower CI values
         y2 = apply(gorg.store,2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store,2, CI)[2,], typ = "l", lty = "dashed", col = "#000000", add = TRUE, lwd = 1)
# Gorgonia RCP8.5
#confplot(x = dates, y1 = apply(gorg.store2,2, CI)[3,], # lower CI values
#         y2 = apply(gorg.store2,2, CI)[1,], #upper CI values
#         col = redtrans, # polygon colour
#         add = TRUE)
## add the sample mean
#points(x = dates, apply(gorg.store,2, CI)[2,], typ = "l", lty = "dotted", col = "firebrick4", add = TRUE, lwd = 1)
# add legend!!!!
legend(x = 2030, y = 300, 
       legend = c("Eunicea flexuosa",#"Eunicea flexuosa",
                  "Antilogorgia americana",#"Antilogorgia americana RCP8.5",
                  "Gorgonia sp."),#,"Gorgonia sp. RCP8.5"), 
       fill = c(bluetrans, #bluetrans, 
                redtrans, #goldtrans, 
                blacktrans), #redtrans),
       col = c("#0072B2", #"blue",
               "#D55E00", #"goldenrod4", 
               "#000000"), #"firebrick4"),
       lty = c("dashed", #"dotted",
               "dashed", #"dotted",
               "dashed"), #"dotted"),
       lwd = 3,
       cex = 2,
       bty = "n")

# Calculate lambda s & its variance --------------------------------------------------------------------
# First determine the mean lambda s and its variance for reporting in the manuscript
# Antilogorgia - sim1
mean(lambdas.store[,1], na.rm = T); sd(lambdas.store[,1], na.rm = T)
# Eunicea - sim1
mean(lambdas.store[,2], na.rm = T); sd(lambdas.store[,2], na.rm = T)
# Gorgonia - sim1
mean(lambdas.store[,3], na.rm = T); sd(lambdas.store[,3], na.rm = T)
# Antilogorgia - sim2
#mean(lambdas.store[,4], na.rm = T); sd(lambdas.store[,4], na.rm = T)
# Eunicea - sim2
#mean(lambdas.store[,5], na.rm = T); sd(lambdas.store[,5], na.rm = T)
# Gorgonia - sim2
#mean(lambdas.store[,6], na.rm = T); sd(lambdas.store[,6], na.rm = T)

# and run a one way anova to determine if there is a significant difference in the simulated trajectories for each population
lambdas.store1 <- melt(lambdas.store[,1:3])
lambdas.store1[,2] <- as.character(lambdas.store1[,2])
# rename the variables
# Species ID
lambdas.store1[grep("Antillogorgia", lambdas.store1[,2]),2] <- "Antillogorgia"
lambdas.store1[grep("Eunicea", lambdas.store1[,2]),2] <- "Eunicea"
lambdas.store1[grep("Gorgonia", lambdas.store1[,2]),2] <- "Gorgonia"
# remove the first column
lambdas.store1 <- lambdas.store1[,2:3]
colnames(lambdas.store1) <- c("Species", "lambdaS")

# run the model
lambdaS.model <- aov(lambdaS ~ Species, data = lambdas.store1)
# and determine the significance.
summary(lambdaS.model) #there is a signficant difference between species and Simulation.
TukeyHSD(lambdaS.model)
# Each species is significantly different from each other.
