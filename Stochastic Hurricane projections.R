# This script is for running a series of simulations exploring the impacts of changing recurrent hurricane regimes on the dynamics and viability of
# 3 different gorgonian populations: Gorgonia, Eunicea and Antillogorgia.

# Author: James Cant
# Contact: james.cant91@gmail.com
#-----------------------------------

# Load required packages
library(areaplot)

# redefine the random seed
set.seed(44736)

###############################################
# STEP 1: Calculate final model projection parameters
###############################################

# determine bin widths for integrating the population
h.store = (U.list-L.list)/m # bin widths

# Determine the initial population integrals (for initiating the simulations)
# Use the bin widths to determine the bin boundaries for fitting observed populations structure (in 2018) into integral format
bin.bounds <- list()

for(i in 1:3){
  bin.bounds[[i]] <- L.list[i] + c(0:m) * (U.list[i] - L.list[i])/m
}

# Now fit observed population vector into a format that matches the bins of a population with a larger size range.
# using the bin boundaries as histogram vectors will achieve this.
initial.pop.integral <- list()

for(i in 1:3) {
  n0 <- hist(initial.pop[[i]], breaks = bin.bounds[[i]], plot = F)$counts
  # then using this line of code converts our population vector into a integration kernel
  n0 = length(initial.pop[[i]]) * n0/(sum(h.store[i] * n0))
  initial.pop.integral[[i]] <- n0
}
# all of the storage vectors are ordered as 1 = Antillogorgia, 2 = Eunicea, 3 = Gorgonia.

# define simulation details
dates <- 2019:2100 # Projection dates
tmax = length(dates) # iteration length
sim <- 1000 # number of simulations

# set up storage matrices for each species under each simulation
# One for population integral
ant.nt <- flex.nt <- gorg.nt <- list(matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax))
# And one for population vector
ant.vec <- flex.vec <- gorg.vec <- list(matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax))
# and set up storage for total population size
ant.Nt <- flex.Nt <- gorg.Nt <- list(numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax))

# Add initial population state for each projection
# initial population structure (nt)
ant.nt[[1]][,1] <- ant.nt[[2]][,1] <- ant.nt[[3]][,1] <- ant.nt[[4]][,1] <- ant.nt[[5]][,1] <- ant.nt[[6]][,1] <- ant.nt[[7]][,1] <- initial.pop.integral[[1]]
flex.nt[[1]][,1] <- flex.nt[[2]][,1] <- flex.nt[[3]][,1] <- flex.nt[[4]][,1] <- flex.nt[[5]][,1] <- flex.nt[[6]][,1] <- flex.nt[[7]][,1] <- initial.pop.integral[[2]]
gorg.nt[[1]][,1] <- gorg.nt[[2]][,1] <- gorg.nt[[3]][,1] <- gorg.nt[[4]][,1] <- gorg.nt[[5]][,1] <- gorg.nt[[6]][,1] <- gorg.nt[[7]][,1] <- initial.pop.integral[[3]]
# initial population size (Nt)
ant.Nt[[1]][1] <- ant.Nt[[2]][1] <- ant.Nt[[3]][1] <- ant.Nt[[4]][1] <- ant.Nt[[5]][1] <- ant.Nt[[6]][1] <- ant.Nt[[7]][1] <- h.store[1] * sum(initial.pop.integral[[1]])
flex.Nt[[1]][1] <- flex.Nt[[2]][1] <- flex.Nt[[3]][1] <- flex.Nt[[4]][1] <- flex.Nt[[5]][1] <- flex.Nt[[6]][1] <- flex.Nt[[7]][1] <- h.store[2] * sum(initial.pop.integral[[2]])
gorg.Nt[[1]][1] <- gorg.Nt[[2]][1] <- gorg.Nt[[3]][1] <- gorg.Nt[[4]][1] <- gorg.Nt[[5]][1] <- gorg.Nt[[6]][1] <- gorg.Nt[[7]][1] <- h.store[3] * sum(initial.pop.integral[[3]]) 
# initial population vector
ant.vec[[1]][,1] <- ant.vec[[2]][,1] <- ant.vec[[3]][,1] <- ant.vec[[4]][,1] <- ant.vec[[5]][,1] <- ant.vec[[6]][,1] <- ant.vec[[7]][,1] <- initial.pop.integral[[1]]/sum(initial.pop.integral[[1]]/sum(h.store[1]*initial.pop.integral[[1]]))
flex.vec[[1]][,1] <- flex.vec[[2]][,1] <- flex.vec[[3]][,1] <- flex.vec[[4]][,1] <- flex.vec[[5]][,1] <- flex.vec[[6]][,1] <- flex.vec[[7]][,1] <- initial.pop.integral[[2]]/sum(initial.pop.integral[[2]]/sum(h.store[2]*initial.pop.integral[[2]]))
gorg.vec[[1]][,1] <- gorg.vec[[2]][,1] <- gorg.vec[[3]][,1] <- gorg.vec[[4]][,1] <- gorg.vec[[5]][,1] <- gorg.vec[[6]][,1] <- gorg.vec[[7]][,1] <- initial.pop.integral[[3]]/sum(initial.pop.integral[[3]]/sum(h.store[3]*initial.pop.integral[[3]]))

# Because recruitment also relies on total colony density this needs saving.
total.Nt <- list(numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax), numeric(tmax))
total.Nt[[1]][1] <- total.Nt[[2]][1] <- total.Nt[[3]][1] <- total.Nt[[4]][1] <- total.Nt[[5]][1] <- total.Nt[[6]][1] <- total.Nt[[7]][1] <- sum(ant.Nt[[1]][1], flex.Nt[[1]][1], gorg.Nt[[1]][1])

# Now ready for the moment of truth.
###############################################
# STEP 2: Run the simulations
###############################################
# iterate t = 1 to tmax, one column per year - after each iteration the loop needs to stop to calculate the current populations state.
# Define the various hurricane simulations that will be used.

# Assign hurricane frequency
# Scenario 1: No Hurricanes
hurricane.prob <- rep(0,length.out = tmax)

# Remaining scenarios:
# Hurricanes Maria and Irma passed the US Virgin Islands as Category 2 Hurricanes with wind speeds of 106Mph. 
# The return time of storm impacts of this magnitude is highly variable and difficult to predict (Keim et al 2011; Knutson et al 2019).
# However, the return time for damaging hurricane storms (Category 1+) with the Caribbean is considered to vary between:
Keim_data <- c(12,12,12,8,15,10,7,10,8,6,8,8,10,35,26,9,7,5)
min_return <- 1/min(Keim_data); max_return <- 1/max(Keim_data) # 0.029 and 0.2
# with an average of:
av_return <- 1/mean(Keim_data) # ~0.0865% likelihood.
# the response of hurricane frequency to climate change is widely debated.
# For instance Knutson (2019) declare that climate change in the North Atlantic will result in a 20% decline in the frequency of hurricane storms
# whilst Bhatia (2018) show that climate change will cause a 12.5% increase by 2035 and a 25% by 2100.
# The remaining scenarios will be used to reflect this impact.

# Scenarios 2,3 & 4: declining Hurricane frequency 
# Create probability vector representing the projected decline in Hurricane frequencies reported in Knutson 2019 for each of the above return time regimes.
# Minimum estimated return time regime:
min.mid.century <- seq(min_return, min_return-(min_return*0.2), length.out = length(2019:2055)) # this reflects a 20% decline by 2055
# but what annual change  does this equate to?
min.ann.change <- (max(min.mid.century)-min(min.mid.century))/length(min.mid.century)
# now ensure the probability of hurricanes changes annually at the same rate until the end of the century.
hurricane.prob2 <- c(min.mid.century, seq(min(min.mid.century)-min.ann.change, min(min.mid.century)-(min.ann.change*(tmax-length(min.mid.century))), length.out = tmax-length(min.mid.century))) # hurricane declines will continue at the same rate until the end of the time series
# and repeat for the other return time regimes:
# Maximum estimated return time regime:
max.mid.century <- seq(max_return, max_return-(max_return*0.2), length.out = length(2019:2055))
max.ann.change <- (max(max.mid.century)-min(max.mid.century))/length(max.mid.century)
hurricane.prob3 <- c(max.mid.century, seq(min(max.mid.century)-max.ann.change, min(max.mid.century)-(max.ann.change*(tmax-length(max.mid.century))), length.out = tmax-length(max.mid.century))) 
# Mean estimated return time regime:
av.mid.century <- seq(av_return, av_return-(av_return*0.2), length.out = length(2019:2055))
av.ann.change <- (max(av.mid.century)-min(av.mid.century))/length(av.mid.century)
hurricane.prob4 <- c(av.mid.century, seq(min(av.mid.century)-av.ann.change, min(av.mid.century)-(av.ann.change*(tmax-length(av.mid.century))), length.out = tmax-length(av.mid.century))) 

# Scenarios 5,6, & 7: Increasing Hurricane occurrence
# Create probability vector representing the projected increase in Hurricane frequencies reported in Bhatia 2018 for each of the above return time regimes.
# Minimum estimated return time regime:
min.2035 <- seq(min_return, min_return+(min_return*0.125), length.out = length(2019:2035)) # this reflects a 12.5% increase by 2035
min.2100 <- seq(max(min.2035), min_return+(min_return*0.25), length.out = length(2035:2100)) # this continues the increase reaching 25% by 2100
# but what annual change  does this equate to?
min.ann.change <- (max(min.2100)-min(min.2100))/length(min.2100)
# now ensure the probability of hurricanes changes annually at the same rate until the end of the time series.
hurricane.prob5 <- c(min.2035, min.2100[-1], seq(max(min.2100)+min.ann.change, max(min.2100)+(min.ann.change*(tmax-(length(min.2035)+length(min.2100[-1])))), length.out = tmax-(length(min.2035)+length(min.2100[-1])))) # hurricane increases will continue at the same rate until the end of the time series
# and repeat for the other return time regimes:
# Maximum estimated return time regime:
max.2035 <- seq(max_return, max_return+(max_return*0.125), length.out = length(2019:2035)) 
max.2100 <- seq(max(max.2035), max_return+(max_return*0.25), length.out = length(2035:2100))
max.ann.change <- (max(max.2100)-min(max.2100))/length(max.2100)
hurricane.prob6 <- c(max.2035, max.2100[-1], seq(max(max.2100)+max.ann.change, max(max.2100)+(max.ann.change*(tmax-(length(max.2035)+length(max.2100[-1])))), length.out = tmax-(length(max.2035)+length(max.2100[-1])))) 
# Mean estimated return time regime:
av.2035 <- seq(av_return, av_return+(av_return*0.125), length.out = length(2019:2035)) 
av.2100 <- seq(max(av.2035), av_return+(av_return*0.25), length.out = length(2035:2100)) 
av.ann.change <- (max(av.2100)-min(av.2100))/length(av.2100)
hurricane.prob7 <- c(av.2035, av.2100[-1], seq(max(av.2100)+av.ann.change, max(av.2100)+(av.ann.change*(tmax-(length(av.2035)+length(av.2100[-1])))), length.out = tmax-(length(av.2035)+length(av.2100[-1]))))

# create overall storage dataframes
ant.store <- flex.store <- gorg.store <- replicate(7, matrix(NA, sim, tmax), simplify = FALSE) 
colnames(ant.store[[1]]) <- colnames(flex.store[[1]]) <- colnames(gorg.store[[1]]) <- colnames(ant.store[[2]]) <- colnames(flex.store[[2]]) <- colnames(gorg.store[[2]]) <-colnames(ant.store[[3]]) <- colnames(flex.store[[3]]) <- colnames(gorg.store[[3]]) <- dates
colnames(ant.store[[4]]) <- colnames(flex.store[[4]]) <- colnames(gorg.store[[4]]) <- colnames(ant.store[[5]]) <- colnames(flex.store[[5]]) <- colnames(gorg.store[[5]]) <-colnames(ant.store[[6]]) <- colnames(flex.store[[6]]) <- colnames(gorg.store[[6]]) <- dates
colnames(ant.store[[7]]) <- colnames(flex.store[[7]]) <- colnames(gorg.store[[7]]) <- dates
# somewhere to store recruits
recruit.store <- replicate(3*7, matrix(NA, sim, tmax), simplify = FALSE)
# somewhere to store population vectors
ant.vec.store <- replicate(7, matrix(NA, m, tmax), simplify = FALSE)
flex.vec.store <- replicate(7, matrix(NA, m, tmax), simplify = FALSE)
gorg.vec.store <- replicate(7, matrix(NA, m, tmax), simplify = FALSE)

# somewhere to store lambda s
lambdas.store <- matrix(NA, nrow = sim, ncol = 3*7) # seven columns for each species (one for each scenario)
colnames(lambdas.store) <- c("Antillogorgia.1", "Eunicea.1", "Gorgonia.1", "Antillogorgia.2", "Eunicea.2", "Gorgonia.2", "Antillogorgia.3", "Eunicea.3", "Gorgonia.3",
                             "Antillogorgia.4", "Eunicea.4", "Gorgonia.4", "Antillogorgia.5", "Eunicea.5", "Gorgonia.5", "Antillogorgia.6", "Eunicea.6", "Gorgonia.6",
                             "Antillogorgia.7", "Eunicea.7", "Gorgonia.7")

#--------------------------------------------------------------------------------------------------------------------------------

# run a grand loop covering all three scenarios for each of the species
for(ii in 1:7){
  # Progress read out
  cat("Intiating Scenario", ii, "\n")
  # Extract associated hurricane potential for the selected scenario
  if(ii == 1){ Hurricane <- hurricane.prob }
  if(ii == 2){ Hurricane <- hurricane.prob2 }
  if(ii == 3){ Hurricane <- hurricane.prob3 }
  if(ii == 4){ Hurricane <- hurricane.prob4 }
  if(ii == 5){ Hurricane <- hurricane.prob5 }
  if(ii == 6){ Hurricane <- hurricane.prob6 }
  if(ii == 7){ Hurricane <- hurricane.prob7 }
  
  # Initiate iterative simulations for the selected scenario (used for determining confidence bounds)
  for(s in 1:sim) {
    # progress read out
    print(s)
    
    # define a vector for storing whether a hurricane occurred in the previous year during each forecast
    Recovery.years <- numeric(tmax)
    
    # Initiate forecast for current simulation
    for(t in 2:tmax) {
      
      # for each step of the forecast first determining if a hurricane should happen
      prob <- runif(1, min = 0, max = 1)
      
      # when the model is in a hurricane year - use the parameters measured during 2017/18
      if (prob < Hurricane[t]){ 
        # the model will run through each species separately
        # Antillogorgia
        P.a <- mk_P_ceiling(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1], U1 = ant.params[15,5]) # calculate the IPM matrix for the singular time step.
        Temp_a <- Iterate(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,5])
        ant.nt[[ii]][,t] <- Temp_a$main 
        ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
        recruit.store[[ii]][s,t] <- Temp_a$R
        ant.vec[[ii]][,t] <- ant.nt[[ii]][,t]/sum(ant.nt[[ii]][,t]/sum(h.store[1]*ant.nt[[ii]][,t]))
          
        # Eunicea
        P.e <- mk_P_ceiling(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2], U1 = flex.params[15,4]) # calculate the IPM matrix for the singular time step.
        Temp_e <- Iterate(nt = flex.nt[[ii]][,t-1], Nt = flex.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,4])
        flex.nt[[ii]][,t] <- Temp_e$main 
        flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
        recruit.store[[ii+7]][s,t] <- Temp_e$R
        flex.vec[[ii]][,t] <- flex.nt[[ii]][,t]/sum(flex.nt[[ii]][,t]/sum(h.store[2]*flex.nt[[ii]][,t]))
        
        # Gorgonia
        P.g <- mk_P_ceiling(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3], U1 = gorg.params[15,5]) # calculate the IPM matrix for the singular time step.
        Temp_g <- Iterate(nt = gorg.nt[[ii]][,t-1], Nt = gorg.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,5])
        gorg.nt[[ii]][,t] <- Temp_g$main 
        gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
        recruit.store[[ii+14]][s,t] <- Temp_g$R
        gorg.vec[[ii]][,t] <- gorg.nt[[ii]][,t]/sum(gorg.nt[[ii]][,t]/sum(h.store[3]*gorg.nt[[ii]][,t]))
        
        # record that a hurricane occurred
        Recovery.years[t] <- 1
        
        # if a hurricane does not occur the model needs to determine if it is in a recovery state or not
      } else {
        if (Recovery.years[t-1] == 1){ #the year following a hurricane is recovery dynamics
          # Antillogorgia
          P.a <- mk_P_ceiling(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1], U1 = ant.params[15,6]) # calculate the IPM matrix for the singular time step.
          Temp_a <- Iterate(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6])
          ant.nt[[ii]][,t] <- Temp_a$main 
          ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
          recruit.store[[ii]][s,t] <- Temp_a$R
          ant.vec[[ii]][,t] <- ant.nt[[ii]][,t]/sum(ant.nt[[ii]][,t]/sum(h.store[1]*ant.nt[[ii]][,t]))
          
          # Eunicea
          P.e <- mk_P_ceiling(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2], U1 = flex.params[15,5]) # calculate the IPM matrix for the singular time step.
          Temp_e <- Iterate(nt = flex.nt[[ii]][,t-1], Nt = flex.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5])
          flex.nt[[ii]][,t] <- Temp_e$main 
          flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
          recruit.store[[ii+7]][s,t] <- Temp_e$R
          flex.vec[[ii]][,t] <- flex.nt[[ii]][,t]/sum(flex.nt[[ii]][,t]/sum(h.store[2]*flex.nt[[ii]][,t]))
          
          # Gorgonia
          P.g <- mk_P_ceiling(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3], U1 = gorg.params[15,6]) # calculate the IPM matrix for the singular time step.
          Temp_g <- Iterate(nt = gorg.nt[[ii]][,t-1], Nt = gorg.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6])
          gorg.nt[[ii]][,t] <- Temp_g$main 
          gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
          recruit.store[[ii+14]][s,t] <- Temp_g$R
          gorg.vec[[ii]][,t] <- gorg.nt[[ii]][,t]/sum(gorg.nt[[ii]][,t]/sum(h.store[3]*gorg.nt[[ii]][,t]))
        
        } else { # if the simulation is not in a hurricane year, nor a recovery year, then the loop will randomly select transitions from the "good" years
          m.par.ant <- ant.params[,sample(1:4,1)] #randomly selects a column from the good year columns of the m.par matrices
          m.par.flex <- flex.params[,sample(1:3,1)]
          m.par.gorg <- gorg.params[,sample(1:4,1)]
          
          # Antillogorgia
          P.a <- mk_P_ceiling(m = m, m.par = m.par.ant, L = L.list[1], U = U.list[1], U1 = m.par.ant[15]) # calculate the IPM matrix for the singular time step.
          Temp_a <- Iterate(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = m.par.ant)
          ant.nt[[ii]][,t] <- Temp_a$main 
          ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
          recruit.store[[ii]][s,t] <- Temp_a$R
          ant.vec[[ii]][,t] <- ant.nt[[ii]][,t]/sum(ant.nt[[ii]][,t]/sum(h.store[1]*ant.nt[[ii]][,t]))
          
          # Eunicea
          P.e <- mk_P_ceiling(m = m, m.par = m.par.flex, L = L.list[2], U = U.list[2], U1 = m.par.flex[15]) # calculate the IPM matrix for the singular time step.
          Temp_e <- Iterate(nt = flex.nt[[ii]][,t-1], Nt = flex.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = m.par.flex)
          flex.nt[[ii]][,t] <- Temp_e$main 
          flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
          recruit.store[[ii+7]][s,t] <- Temp_e$R
          flex.vec[[ii]][,t] <- flex.nt[[ii]][,t]/sum(flex.nt[[ii]][,t]/sum(h.store[2]*flex.nt[[ii]][,t]))
          
          # Gorgonia
          P.g <- mk_P_ceiling(m = m, m.par = m.par.gorg, L = L.list[3], U = U.list[3], U1 = m.par.gorg[15]) # calculate the IPM matrix for the singular time step.
          Temp_g <- Iterate(nt = gorg.nt[[ii]][,t-1], Nt = gorg.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = m.par.gorg)
          gorg.nt[[ii]][,t] <- Temp_g$main 
          gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
          recruit.store[[ii+14]][s,t] <- Temp_g$R
          gorg.vec[[ii]][,t] <- gorg.nt[[ii]][,t]/sum(gorg.nt[[ii]][,t]/sum(h.store[3]*gorg.nt[[ii]][,t]))
        }
      }
      # at the end of each loop save the total populations size
      total.Nt[[ii]][t] <- ant.Nt[[ii]][t] + flex.Nt[[ii]][t] + gorg.Nt[[ii]][t]
    }
    # Store the first simulation run through
    ant.store[[ii]][s,] <- ant.Nt[[ii]]
    flex.store[[ii]][s,] <- flex.Nt[[ii]]
    gorg.store[[ii]][s,] <- gorg.Nt[[ii]]
    
    # Re-estimate and store the mean population vector for each time step
    Ant.array <- array(c(ant.vec.store[[ii]], ant.vec[[ii]]), c(m, tmax, 2))
    ant.vec.store[[ii]] <- rowMeans(Ant.array, dims = 2, na.rm = T)
    Flex.array <- array(c(flex.vec.store[[ii]], flex.vec[[ii]]), c(m, tmax, 2))
    flex.vec.store[[ii]] <- rowMeans(Flex.array, dims = 2, na.rm = T)
    Gorg.array <- array(c(gorg.vec.store[[ii]], gorg.vec[[ii]]), c(m, tmax, 2))
    gorg.vec.store[[ii]] <- rowMeans(Gorg.array, dims = 2, na.rm = T)
    
    # and calculate lambda s for each simulation
    # create blank storage  
    lambdas.ant <- numeric(tmax-1) 
    lambdas.flex <- numeric(tmax-1) 
    lambdas.gorg <- numeric(tmax-1) 
    # determine stepwise lambda
    for (x in 1: tmax-1) {
      lambdas.ant[x] <- ant.Nt[[ii]][x+1]/ant.Nt[[ii]][x]
      lambdas.flex[x] <- flex.Nt[[ii]][x+1]/flex.Nt[[ii]][x]
      lambdas.gorg[x] <- gorg.Nt[[ii]][x+1]/gorg.Nt[[ii]][x]
    }
    
    # and store overall lambda s for the simulation
    if(ii == 1){
      lambdas.store[s,1] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,2] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,3] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 2){
      lambdas.store[s,4] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,5] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,6] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 3){
      lambdas.store[s,7] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,8] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,9] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 4){
      lambdas.store[s,10] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,11] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,12] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 5){
      lambdas.store[s,13] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,14] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,15] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 6){
      lambdas.store[s,16] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,17] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,18] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
    if(ii == 7){
      lambdas.store[s,19] <- exp(mean(log(lambdas.ant), na.rm = T))
      lambdas.store[s,20] <- exp(mean(log(lambdas.flex), na.rm = T))
      lambdas.store[s,21] <- exp(mean(log(lambdas.gorg), na.rm = T)) }
  }
}

###############################################
# STEP 3: Plot simulations and variance and statistical differences
###############################################

# Plot the simulations -----------------------------------------------------------------------------------
#what is my max y value across each scenario
ymax <- max(round((max(apply(flex.store[[1]], 2, CI)[1,])+50), digits = -2), # Scenario 1. The plus 50 ensures that values are always rounded up.
            round((max(apply(flex.store[[2]], 2, CI)[1,])+50), digits = -2), # Scenario 2
            round((max(apply(flex.store[[3]], 2, CI)[1,])+50), digits = -2), # Scenario 3 
            round((max(apply(flex.store[[4]], 2, CI)[1,])+50), digits = -2), # Scenario 4
            round((max(apply(flex.store[[5]], 2, CI)[1,])+50), digits = -2), # Scenario 5 
            round((max(apply(flex.store[[6]], 2, CI)[1,])+50), digits = -2), # Scenario 6
            round((max(apply(flex.store[[7]], 2, CI)[1,])+50), digits = -2)) # Scenario 7

# Manually defining ymax can help to expand some of the visuals for the manuscript
ymax <- 420

# plot each scenario
##### Scenario 1: No Hurricane
# Eunicea
plot(x = dates, apply(flex.store[[1]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
#set a transparent colour
bluetrans <- rgb(0,114,178, alpha = 50, maxColorValue = 255)
confplot(x = dates, y1 = apply(flex.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[1]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[1]], 2, CI)[2,], typ = "l", lty = "solid", col = "#0072B2", lwd = 1)

# Antilogorgia 
#set a transparent colour
redtrans <- rgb(213,94,0, alpha = 50, maxColorValue = 255)
confplot(x = dates, y1 = apply(ant.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[1]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[1]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)

# Gorgonia 
#set a transparent colour
blacktrans <- rgb(0,0,0, alpha = 50, maxColorValue = 255)
confplot(x = dates, y1 = apply(gorg.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[1]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[1]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Add references lines for the initial population sizes used for each species
segments(x0 = dates[1], y0 = h.store[2] * sum(initial.pop.integral[[2]]), x1 = dates[30], lty = "solid", col = "#0072B2", lwd = 2) # Flexuosa
segments(x0 = dates[1], y0 = h.store[1] * sum(initial.pop.integral[[1]]), x1 = dates[30], lty = "solid", col = "#D55E00", lwd = 2) # Antillogorgia
segments(x0 = dates[1], y0 = h.store[3] * sum(initial.pop.integral[[3]]), x1 = dates[30], lty = "solid", col = "#000000", lwd = 2) # Gorgonia

# add legend!!!!
legend(x = 2030, y = ymax, 
       legend = c("Eunicea flexuosa",#"Eunicea flexuosa",
                  "Antilogorgia americana",#"Antilogorgia americana RCP8.5",
                  "Gorgonia sp."),#,"Gorgonia sp. RCP8.5"), 
       fill = c("#0072B2", #"blue",
               "#D55E00", #"goldenrod4", 
               "#000000"), #"firebrick4"),
      # lty = c("solid", #"dotted",
      #         "solid", #"dotted",
      #         "solid"), #"dotted"),
      # lwd = 3,
       cex = 1.5,
       bty = "n")

# For the remaining plots combine together scenarios dealing with the same initial conditions (e.g. min, mean and max return time)
##### Minimum return time
# Eunicea
# Scenario 2: Decreasing hurricane likelihood
plot(x = dates, apply(flex.store[[2]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Colours remain consistent across plots
confplot(x = dates, y1 = apply(flex.store[[2]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[2]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[2]], 2, CI)[2,], typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
bluetrans2 <- rgb(0,0,255, alpha = 40, maxColorValue = 255)
confplot(x = dates, y1 = apply(flex.store[[5]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[5]], 2, CI)[1,], #upper CI values
         col = bluetrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[5]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)

# Antilogorgia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(ant.store[[2]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[2]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[2]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
redtrans2 <- rgb(255,0,0, alpha = 80, maxColorValue = 255)
confplot(x = dates, y1 = apply(ant.store[[5]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[5]], 2, CI)[1,], #upper CI values
         col = redtrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[5]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)

# Gorgonia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(gorg.store[[2]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[2]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[2]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
blacktrans2 <- rgb(0,0,0, alpha = 60, maxColorValue = 255)
confplot(x = dates, y1 = apply(gorg.store[[5]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[5]], 2, CI)[1,], #upper CI values
         col = blacktrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[5]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)
# for publication legend will be copied across using illustrator 


##### Maximum return time
# Eunicea
# Scenario 3: Decreasing hurricane likelihood
plot(x = dates, apply(flex.store[[3]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Colours remain consistent across plots
confplot(x = dates, y1 = apply(flex.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[3]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[3]], 2, CI)[2,], typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Scenario 6: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(flex.store[[6]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[6]], 2, CI)[1,], #upper CI values
         col = bluetrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[6]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)

# Antilogorgia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(ant.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[3]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[3]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(ant.store[[6]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[6]], 2, CI)[1,], #upper CI values
         col = redtrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[6]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)

# Gorgonia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(gorg.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[3]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[3]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(gorg.store[[6]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[6]], 2, CI)[1,], #upper CI values
         col = blacktrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[6]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)
# for publication legend will be copied across using illustrator 


##### Mean return time
# Eunicea
# Scenario 2: Decreasing hurricane likelihood
plot(x = dates, apply(flex.store[[4]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Colours remain consistent across plots
confplot(x = dates, y1 = apply(flex.store[[4]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[4]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[4]], 2, CI)[2,], typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(flex.store[[7]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[7]], 2, CI)[1,], #upper CI values
         col = bluetrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[7]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)

# Antilogorgia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(ant.store[[4]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[4]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[4]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(ant.store[[7]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[7]], 2, CI)[1,], #upper CI values
         col = redtrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[7]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)

# Gorgonia 
# Scenario 2: Decreasing likelihood
confplot(x = dates, y1 = apply(gorg.store[[4]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[4]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[4]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
confplot(x = dates, y1 = apply(gorg.store[[7]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[7]], 2, CI)[1,], #upper CI values
         col = blacktrans2, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[7]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)
# for publication legend will be copied across using illustrator 

# Calculate lambda s & its variance --------------------------------------------------------------------
# First determine the mean lambda s and its variance for reporting in the manuscript

# Scenario 1
# Antilogorgia
CI(lambdas.store[,1])
# Eunicea
CI(lambdas.store[,2])
# Gorgonia
CI(lambdas.store[,3])

# Scenario 2
# Antilogorgia
CI(lambdas.store[,4])
# Eunicea
CI(lambdas.store[,5])
# Gorgonia
CI(lambdas.store[,6])

# Scenario 3
# Antilogorgia
CI(lambdas.store[,7])
# Eunicea
CI(lambdas.store[,8])
# Gorgonia
CI(lambdas.store[,9])

# Scenario 4
# Antilogorgia
CI(lambdas.store[,10])
# Eunicea
CI(lambdas.store[,11])
# Gorgonia
CI(lambdas.store[,12])

# Scenario 5
# Antilogorgia
CI(lambdas.store[,13])
# Eunicea
CI(lambdas.store[,14])
# Gorgonia 
CI(lambdas.store[,15])

# Scenario 6
# Antilogorgia
CI(lambdas.store[,16])
# Eunicea
CI(lambdas.store[,17])
# Gorgonia
CI(lambdas.store[,18])

# Scenario 7
# Antilogorgia
CI(lambdas.store[,19])
# Eunicea
CI(lambdas.store[,20])
# Gorgonia
CI(lambdas.store[,21])

# and finally run a two-way anova to determine if there is a significant difference in the simulated trajectories for each population
lambdas.ANOVA.data <- melt(lambdas.store)
# reassign nessecary variable names
lambdas.ANOVA.data[grep("1", lambdas.ANOVA.data$Var2),]$Var1 <- 1
lambdas.ANOVA.data[grep("2", lambdas.ANOVA.data$Var2),]$Var1 <- 2
lambdas.ANOVA.data[grep("3", lambdas.ANOVA.data$Var2),]$Var1 <- 3
lambdas.ANOVA.data[grep("4", lambdas.ANOVA.data$Var2),]$Var1 <- 4
lambdas.ANOVA.data[grep("5", lambdas.ANOVA.data$Var2),]$Var1 <- 5
lambdas.ANOVA.data[grep("6", lambdas.ANOVA.data$Var2),]$Var1 <- 6
lambdas.ANOVA.data[grep("7", lambdas.ANOVA.data$Var2),]$Var1 <- 7
lambdas.ANOVA.data$Var2 <- as.character(lambdas.ANOVA.data$Var2)
lambdas.ANOVA.data[grep("Antillogorgia", lambdas.ANOVA.data$Var2),]$Var2 <- "Antillogorgia"
lambdas.ANOVA.data[grep("Eunicea", lambdas.ANOVA.data$Var2),]$Var2 <- "Eunicea"
lambdas.ANOVA.data[grep("Gorgonia", lambdas.ANOVA.data$Var2),]$Var2 <- "Gorgonia"

# rename and reformat (if nessecary) the columns
colnames(lambdas.ANOVA.data) <- c("Scenario", "Species", "lambdaS")
lambdas.ANOVA.data$Scenario <- as.factor(lambdas.ANOVA.data$Scenario)
lambdas.ANOVA.data$Species <- as.factor(lambdas.ANOVA.data$Species)

# run the model
lambdaS.model <- aov(lambdaS ~ Species + Scenario + Species:Scenario, data = lambdas.ANOVA.data)
# and determine the significance.
car::Anova(lambdaS.model, type = "II")
# and identify where this significance lies
TukeyHSD(lambdaS.model)

### Visualize recruitment variability across years and scenarios -----------------------------

# what is my max y value across each scenario
ymax2 <- ceiling(max(apply(recruit.store[[8]], 2, mean), # Scenario 1
            apply(recruit.store[[9]], 2, mean), # Scenario 2
            apply(recruit.store[[10]], 2, mean), # Scenario 3 
            apply(recruit.store[[11]], 2, mean), # Scenario 4
            apply(recruit.store[[12]], 2, mean), # Scenario 5 
            apply(recruit.store[[13]], 2, mean), # Scenario 6
            apply(recruit.store[[14]], 2, mean), na.rm = TRUE)) # Scenario 7

# Manually defining ymax can help to expand some of the visuals for the manuscript
ymax2 <- 200

# plot recruitment across each scenario
##### Scenario 1: No Hurricane
# Eunicea
plot(x = dates, apply(recruit.store[[8]], 2, max), typ = "n", ylim = c(0, ymax2), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
points(x = dates, apply(recruit.store[[8]], 2, mean), typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Antilogorgia 
points(x = dates, apply(recruit.store[[1]], 2, mean), typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[15]], 2, mean), typ = "l", lty = "solid", col = "#000000", lwd = 1)

##### Minimum return time
plot(x = dates, apply(recruit.store[[8]], 2, mean), typ = "n", ylim = c(0, ymax2), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Scenario 2: Decreasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[9]], 2, mean), typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Antilogorgia 
points(x = dates, apply(recruit.store[[2]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[16]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[12]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)
# Antilogorgia
points(x = dates, apply(recruit.store[[5]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[19]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)

##### Maximum return time
plot(x = dates, apply(recruit.store[[8]], 2, mean), typ = "n", ylim = c(0, ymax2), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Scenario 2: Decreasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[11]], 2, mean), typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Antilogorgia 
points(x = dates, apply(recruit.store[[4]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[18]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[14]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)
# Antilogorgia
points(x = dates, apply(recruit.store[[7]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[21]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)

##### Mean return time
plot(x = dates, apply(recruit.store[[8]], 2, mean), typ = "n", ylim = c(0, ymax2), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
# Scenario 2: Decreasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[10]], 2, mean), typ = "l", lty = "solid", col = "#0072B2", lwd = 1)
# Antilogorgia 
points(x = dates, apply(recruit.store[[3]], 2, CI)[2,], typ = "l", lty = "solid", col = "#D55E00", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[17]], 2, CI)[2,], typ = "l", lty = "solid", col = "#000000", lwd = 1)
# Scenario 5: Increasing hurricane likelihood
# Eunicea
points(x = dates, apply(recruit.store[[13]], 2, CI)[2,], typ = "l", lty = "dashed", col = "blue", lwd = 1)
# Antilogorgia
points(x = dates, apply(recruit.store[[6]], 2, CI)[2,], typ = "l", lty = "dashed", col = "red", lwd = 1)
# Gorgonia 
points(x = dates, apply(recruit.store[[20]], 2, CI)[2,], typ = "l", lty = "dashed", col = "black", lwd = 1)

### Visualize population state distributions across years and scenarios -----------------------------

# First extract a series of colony sizes from bin bounds vectors for use as x-axis labels.
Xaxis_ant <- round(bin.bounds[[1]][c(1,50,100,150,200)])
Xaxis_flex <- round(bin.bounds[[2]][c(1,50,100,150,200)])
Xaxis_gorg <- round(bin.bounds[[3]][c(1,50,100,150,200)])

# First plot size distribution projections for the non-hurricane scenario for each species
# These will go into the main manuscript.

# Antillogorgia
matplot(ant.vec.store[[1]], typ = "n", ylim = c(0, 10), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5,
     xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_ant, cex.axis = 1.5)
# Add size distributions
# 2020
points(ant.vec.store[[1]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(ant.vec.store[[1]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(ant.vec.store[[1]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(ant.vec.store[[1]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(ant.vec.store[[1]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(ant.vec.store[[1]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(ant.vec.store[[1]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Eunicea
matplot(flex.vec.store[[1]], typ = "n", ylim = c(0, 12), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_flex, cex.axis = 1.5)
# Add size distributions
# 2020
points(flex.vec.store[[1]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(flex.vec.store[[1]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(flex.vec.store[[1]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(flex.vec.store[[1]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(flex.vec.store[[1]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(flex.vec.store[[1]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(flex.vec.store[[1]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Gorgonia
matplot(gorg.vec.store[[1]], typ = "n", ylim = c(0, 10), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_gorg, cex.axis = 1.5)
# Add size distributions
# 2020
points(gorg.vec.store[[1]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(gorg.vec.store[[1]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(gorg.vec.store[[1]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(gorg.vec.store[[1]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(gorg.vec.store[[1]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(gorg.vec.store[[1]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(gorg.vec.store[[1]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Next plot all remaining scenarios for Supplementary
# Scenario 2
# Antillogorgia
matplot(ant.vec.store[[2]], typ = "n", ylim = c(0, 11), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_ant, cex.axis = 1.5)
# Add size distributions
# 2020
points(ant.vec.store[[2]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(ant.vec.store[[2]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(ant.vec.store[[2]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(ant.vec.store[[2]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(ant.vec.store[[2]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(ant.vec.store[[2]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(ant.vec.store[[2]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Eunicea
matplot(flex.vec.store[[2]], typ = "n", ylim = c(0, 12), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_flex, cex.axis = 1.5)
# Add size distributions
# 2020
points(flex.vec.store[[2]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(flex.vec.store[[2]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(flex.vec.store[[2]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(flex.vec.store[[2]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(flex.vec.store[[2]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(flex.vec.store[[2]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(flex.vec.store[[2]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Gorgonia
matplot(gorg.vec.store[[2]], typ = "n", ylim = c(0, 10), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_gorg, cex.axis = 1.5)
# Add size distributions
# 2020
points(gorg.vec.store[[2]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(gorg.vec.store[[2]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(gorg.vec.store[[2]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(gorg.vec.store[[2]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(gorg.vec.store[[2]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(gorg.vec.store[[2]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(gorg.vec.store[[2]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Scenarios 3 - 7
# Antillogorgia
matplot(ant.vec.store[[7]], typ = "n", ylim = c(0, 8), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_ant, cex.axis = 1.5)
# Add size distributions
# 2020
points(ant.vec.store[[7]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(ant.vec.store[[7]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(ant.vec.store[[7]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(ant.vec.store[[7]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(ant.vec.store[[7]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(ant.vec.store[[7]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(ant.vec.store[[7]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Eunicea
matplot(flex.vec.store[[7]], typ = "n", ylim = c(0, 9), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_flex, cex.axis = 1.5)
# Add size distributions
# 2020
points(flex.vec.store[[7]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(flex.vec.store[[7]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(flex.vec.store[[7]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(flex.vec.store[[7]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(flex.vec.store[[7]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(flex.vec.store[[7]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(flex.vec.store[[7]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)

# Gorgonia
matplot(gorg.vec.store[[7]], typ = "n", ylim = c(0, 14), yaxs = "i", xaxs = "i",
        xlab = "",
        ylab = "",
        cex.axis = 1.5,
        xaxt = 'n')
# Re-label x axis tick labels.
axis(1, at = c(1, 50, 100, 150, 200), labels = Xaxis_gorg, cex.axis = 1.5)
# Add size distributions
# 2020
points(gorg.vec.store[[7]][,2], typ = "l", lty = "solid", col = "black", lwd = 3)
# 2023
points(gorg.vec.store[[7]][,5], typ = "l", lty = "solid", col = "dark grey", lwd = 3)
# 2039
points(gorg.vec.store[[7]][,21], typ = "l", lty = "solid", col = "light grey", lwd = 3)
#2055
points(gorg.vec.store[[7]][,37], typ = "l", lty = "solid", col = "blue", lwd = 3)
# 2071
points(gorg.vec.store[[7]][,53], typ = "l", lty = "solid", col = "red", lwd = 3)
# 2087
points(gorg.vec.store[[7]][,69], typ = "l", lty = "solid", col = "purple", lwd = 3)
# 2100
points(gorg.vec.store[[7]][,82], typ = "l", lty = "solid", col = "orange", lwd = 3)
legend(150, 10, legend = c("2020","2023", "2039","2055","2071","2087","2100"), fill = c("Black","dark grey","light grey","blue","red","purple",'orange'), bg = NULL, bty = "n", cex = 1.2)



##### ----------------------------------------------------- End of Code & Analysis ---------------------------------------------
