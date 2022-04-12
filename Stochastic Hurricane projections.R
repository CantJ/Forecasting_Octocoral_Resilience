# This script is for running a series of simulations exploring the impacts of changing recurrent hurricane regimes on the dynamics and viability of
# 3 different gorgonia populations: Gorgonia, Eunicea and Antillogorgia.
# 
# data last modified: Mar 2022
# Author: James Cant

# Load required packages
library(areaplot)

# redefine the random seed
set.seed(44736)

###############################################
# STEP 1: Calculate final model projection parameters
###############################################
# reset plot window
set_graph_pars("panel1")

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
dates <- 2018:2150 # Projection dates
tmax = length(dates) # iteration length
sim <- 1000 # number of simulations

# set up one storage matrix for each species under each simulation
ant.nt <- flex.nt <- gorg.nt <- list(matrix(NA, m, tmax), matrix(NA, m, tmax), matrix(NA, m, tmax))
# and set up storage for total population size
ant.Nt <- flex.Nt <- gorg.Nt <- list(numeric(tmax), numeric(tmax), numeric(tmax))

# Add initial population state for each projection
# initial population structure (nt)
ant.nt[[1]][,1] <- ant.nt[[2]][,1] <- ant.nt[[3]][,1] <- initial.pop.integral[[1]]
flex.nt[[1]][,1] <- flex.nt[[2]][,1] <- flex.nt[[3]][,1] <- initial.pop.integral[[2]]
gorg.nt[[1]][,1] <- gorg.nt[[2]][,1] <- gorg.nt[[3]][,1] <- initial.pop.integral[[3]]
# initial populatioe size (Nt)
ant.Nt[[1]][1] <- ant.Nt[[2]][1] <- ant.Nt[[3]][1] <- h.store[1] * sum(initial.pop.integral[[1]])
flex.Nt[[1]][1] <- flex.Nt[[2]][1] <- flex.Nt[[3]][1] <- h.store[2] * sum(initial.pop.integral[[2]])
gorg.Nt[[1]][1] <- gorg.Nt[[2]][1] <- gorg.Nt[[3]][1] <- h.store[3] * sum(initial.pop.integral[[3]]) 

# Because recruitment also relies on total colony density this needs saving.
total.Nt <- list(numeric(tmax), numeric(tmax), numeric(tmax))
total.Nt[[1]][1] <- total.Nt[[2]][1] <- total.Nt[[3]][1] <- sum(ant.Nt[[1]][1], flex.Nt[[1]][1], gorg.Nt[[1]][1])

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
# A little fix to prevent impossible probabilities
hurricane.prob2[which(hurricane.prob2 >= 1)] <- 0.99

# Scenario 3: Increasing Hurricane occurrence
# Artificially define hurricane probability to increase by same rate as it is expected to decrease.
hurricane.prob3 <- seq(from = 0.65, by = ann.change, length.out = tmax)
# Again prevent improbable hurricane scenarios
hurricane.prob3[which(hurricane.prob3 >= 1)] <- 0.99

# create overall storage dataframe
ant.store <- flex.store <- gorg.store <- list(matrix(NA, sim, tmax), matrix(NA, sim, tmax), matrix(NA, sim, tmax))
colnames(ant.store[[1]]) <- colnames(flex.store[[1]]) <- colnames(gorg.store[[1]]) <- colnames(ant.store[[2]]) <- colnames(flex.store[[2]]) <- colnames(gorg.store[[2]]) <-colnames(ant.store[[3]]) <- colnames(flex.store[[3]]) <- colnames(gorg.store[[3]]) <- dates
# somewhere to store lambda s
lambdas.store <- matrix(NA, nrow = sim, ncol = 9) # three columns for each species (one for sim1, sim2, and sim3)
colnames(lambdas.store) <- c("Antillogorgia.1", "Eunicea.1", "Gorgonia.1", "Antillogorgia.2", "Eunicea.2", "Gorgonia.2", "Antillogorgia.3", "Eunicea.3", "Gorgonia.3")

#--------------------------------------------------------------------------------------------------------------------------------

# run a grand loop covering all three scenarios for each of the species
for(ii in 1:3){
  # Progress read out
  cat("Intiating Scenario", ii, "\n")
  # Extract associated hurricane potential for the selected scenario
  if(ii == 1){ Hurricane <- hurricane.prob }
  if(ii == 2){ Hurricane <- hurricane.prob2 }
  if(ii == 3){ Hurricane <- hurricane.prob3 }
  
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
      
      # The population data being used to define initial population structure was collected in 2018 (not 2019, from which we have demographic data) 
      # therefore the simulation is forced to follow the dynamics calculated for the 2018/19 period in order to better reflect reality.
      if(t == 2) {
        # Here all IPMs will be constructed using 2018/19 coefficients.
        # Antillogorgia
        P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
        ant.nt[[ii]][,t] <- Iterate.a(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
        ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
        
        # Eunicea
        P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
        flex.nt[[ii]][,t] <- Iterate.e(nt = flex.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
        flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
        
        # Gorgonia
        P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
        gorg.nt[[ii]][,t] <- Iterate.g(nt = gorg.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
        gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
        
      } else { 
        # otherwise the simulation uses the preassigned probabilities to determine hurricane liklihood.
        # when the model is in a hurricane year - use the parameters measured during 2017/18
        if (prob < Hurricane[t]){ 
          # the model will run through each species separately
          # Antillogorgia
          P.a <- mk_P(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
          ant.nt[[ii]][,t] <- Iterate.a(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,5]) 
          ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t-1])
            
          # Eunicea
          P.e <- mk_P(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
          flex.nt[[ii]][,t] <- Iterate.e(nt = flex.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,4]) 
          flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
            
          # Gorgonia
          P.g <- mk_P(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
          gorg.nt[[ii]][,t] <- Iterate.g(nt = gorg.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,5]) 
          gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
            
          # record that a hurricane occurred
          Recovery.years[t] <- 1
          
        # if a hurricane does not occur the model needs to determine if it is in a recovery state or not
        } else {
          if (Recovery.years[t-1] == 1){ #the year following a hurricane is recovery dynamics
            # Antillogorgia
            P.a <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
            ant.nt[[ii]][,t] <- Iterate.a(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = ant.params[,6]) 
            ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
              
            # Eunicea
            P.e <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
            flex.nt[[ii]][,t] <- Iterate.e(nt = flex.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = flex.params[,5]) 
            flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
              
            # Gorgonia
            P.g <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
            gorg.nt[[ii]][,t] <- Iterate.g(nt = gorg.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = gorg.params[,6]) 
            gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
          
          } else { # if the simulation is not in a hurricane year, nor a recovery year, then the loop will randomly select transitions from the "good" years
            m.par.ant <- ant.params[,sample(1:4,1)] #randomly selects a column from the good year columns of the m.par matrices
            m.par.flex <- flex.params[,sample(1:3,1)]
            m.par.gorg <- gorg.params[,sample(1:4,1)]
              
            # Antillogorgia
            P.a <- mk_P(m = m, m.par = m.par.ant, L = L.list[1], U = U.list[1]) # calculate the IPM matrix for the singular time step.
            ant.nt[[ii]][,t] <- Iterate.a(nt = ant.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.a$meshpts, P = P.a$P, m.par = m.par.ant) 
            ant.Nt[[ii]][t] <- h.store[1] * sum(ant.nt[[ii]][,t])
              
            # Eunicea
            P.e <- mk_P(m = m, m.par = m.par.flex, L = L.list[2], U = U.list[2]) # calculate the IPM matrix for the singular time step.
            flex.nt[[ii]][,t] <- Iterate.e(nt = flex.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.e$meshpts, P = P.e$P, m.par = m.par.flex) 
            flex.Nt[[ii]][t] <- h.store[2] * sum(flex.nt[[ii]][,t])
              
            # Gorgonia
            P.g <- mk_P(m = m, m.par = m.par.gorg, L = L.list[3], U = U.list[3]) # calculate the IPM matrix for the singular time step.
            gorg.nt[[ii]][,t] <- Iterate.g(nt = gorg.nt[[ii]][,t-1], Nt = total.Nt[[ii]][t-1], meshpts = P.g$meshpts, P = P.g$P, m.par = m.par.gorg) 
            gorg.Nt[[ii]][t] <- h.store[3] * sum(gorg.nt[[ii]][,t])
          }}
      }
      # at the end of each loop save the total populations size
      total.Nt[[ii]][t] <- ant.Nt[[ii]][t] + flex.Nt[[ii]][t] + gorg.Nt[[ii]][t]
    }
    # Store the first simulation run through
    ant.store[[ii]][s,] <- ant.Nt[[ii]]
    flex.store[[ii]][s,] <- flex.Nt[[ii]]
    gorg.store[[ii]][s,] <- gorg.Nt[[ii]]
    
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
  }
}

###############################################
# STEP 3: Plot simulations and variance and statistical differences
###############################################

# Plot the simulations -----------------------------------------------------------------------------------
#what is my max y value across each scenario
ymax <- max(round((max(apply(flex.store[[1]], 2, CI)[1,])+50), digits = -2), # Scenario 1. The plus 50 ensures that values are always rounded up.
            round((max(apply(flex.store[[2]], 2, CI)[1,])+50), digits = -2),  # Scenario 2
            round((max(apply(flex.store[[3]], 2, CI)[1,])+50), digits = -2)) # Scenario 3

# plot each scenario
##### Scenario 1: No Hurricane
# Eunicea
plot(x = dates, apply(flex.store[[1]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
#set a transparent colour
bluetrans <- rgb(0,114,178, alpha = 100, maxColorValue = 255)
confplot(x = dates, y1 = apply(flex.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[1]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[1]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#0072B2", lwd = 1)

# Antilogorgia 
#set a transparent colour
redtrans <- rgb(213,94,0, alpha = 100, maxColorValue = 255)
confplot(x = dates, y1 = apply(ant.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[1]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[1]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#D55E00", lwd = 1)

# Gorgonia 
#set a transparent colour
blacktrans <- rgb(0,0,0, alpha = 150, maxColorValue = 255)
confplot(x = dates, y1 = apply(gorg.store[[1]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[1]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[1]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#000000", lwd = 1)
# add legend!!!!
legend(x = 2030, y = ymax, 
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
       cex = 1.5,
       bty = "n")

##### Scenario 2: Decreasing hurricane likelihood
# Eunicea
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
points(x = dates, apply(flex.store[[2]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#0072B2", lwd = 1)

# Antilogorgia 
confplot(x = dates, y1 = apply(ant.store[[2]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[2]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[2]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#D55E00", lwd = 1)

# Gorgonia 
confplot(x = dates, y1 = apply(gorg.store[[2]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[2]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[2]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#000000", lwd = 1)
# for publication legend will be copied across using illustrator 

##### Scenario 1: No Hurricane
# Eunicea
plot(x = dates, apply(flex.store[[3]], 2, CI)[2,], typ = "n", ylim = c(0, ymax), yaxs = "i", xaxs = "i",
     xlab = "",
     ylab = "",
     cex.axis = 1.5)
confplot(x = dates, y1 = apply(flex.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(flex.store[[3]], 2, CI)[1,], #upper CI values
         col = bluetrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(flex.store[[3]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#0072B2", lwd = 1)

# Antilogorgia 
confplot(x = dates, y1 = apply(ant.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(ant.store[[3]], 2, CI)[1,], #upper CI values
         col = redtrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(ant.store[[3]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#D55E00", lwd = 1)

# Gorgonia 
confplot(x = dates, y1 = apply(gorg.store[[3]], 2, CI)[3,], # lower CI values
         y2 = apply(gorg.store[[3]], 2, CI)[1,], #upper CI values
         col = blacktrans, # polygon colour
         add = TRUE)
# add the sample mean
points(x = dates, apply(gorg.store[[3]], 2, CI)[2,], typ = "l", lty = "dashed", col = "#000000", lwd = 1)
# Again legend will be manually copied across for publication

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

# and finally run a two-way anova to determine if there is a significant difference in the simulated trajectories for each population
lambdas.ANOVA.data <- melt(lambdas.store)
# reassign nessecary variable names
lambdas.ANOVA.data[grep("1", lambdas.ANOVA.data$Var2),]$Var1 <- 1
lambdas.ANOVA.data[grep("2", lambdas.ANOVA.data$Var2),]$Var1 <- 2
lambdas.ANOVA.data[grep("3", lambdas.ANOVA.data$Var2),]$Var1 <- 3
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

##### ----------------------------------------------------- End of Code & Analysis ---------------------------------------------
