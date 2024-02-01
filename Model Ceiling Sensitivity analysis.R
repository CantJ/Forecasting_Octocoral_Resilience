# This script is for running a sensitivity analysis on the position of the P kernel ceiling and its effect on LambdaP across the 
# 3 different gorgonian populations: Gorgonia, Eunicea and Antillogorgia.

# Author: James Cant
# Contact: james.cant91@gmail.com
#-----------------------------------

# Load required packages
library(tidyr) 
library(dplyr)

# redefine the random seed
set.seed(500123)

#####################################
# STEP 1: Define variance in ceiling value
#####################################

# For this sensitivity analysis all model parameters will be held constant except for the species specific value of U1.
# over a series of iterations, U1 will be altered for each species, with the yearly values of lambda P calculated and stored.

# Firstly define a vector of U1 values for each species.
# This vector should be defined such that the 'baseline' U1 value for each species is the median value.
# For each species values within this vector also cannot exceed U.

# View baseline U1 values.
ant.params['U1',][1]; flex.params['U1',][1]; gorg.params['U1',][1]

# define vectors (adding in also the 'true' U1 value)
ant.vec <- c(seq((ant.params['U1',][1])-0.6, (ant.params['U1',][1])+0.6, by = 0.01), ant.params['U1',][1])
flex.vec <- c(seq((flex.params['U1',][1])-0.6, (flex.params['U1',][1])+0.6, by = 0.01), flex.params['U1',][1])
gorg.vec <- c(seq((gorg.params['U1',][1])-0.6, (gorg.params['U1',][1])+0.6, by = 0.01), gorg.params['U1',][1])

# Ensure vectors are in numerical order for plotting
ant.vec <- sort(ant.vec, decreasing = FALSE)
flex.vec <- sort(flex.vec, decreasing = FALSE)
gorg.vec <- sort(gorg.vec, decreasing = FALSE)

#####################################
# STEP 2: Run interative sensitivity loop
#####################################
# This loop will simply generate a series of IPMs each constructed with differing ceiling values and estimate and store their lambda values.

# Generate storage
# Antillogorgia
ant.sens <- data.frame(U1 = ant.vec,
                       Lambda_2013 = rep(NA, length.out = length(ant.vec)),
                       Lambda_2014 = rep(NA, length.out = length(ant.vec)),
                       Lambda_2015 = rep(NA, length.out = length(ant.vec)),
                       Lambda_2016 = rep(NA, length.out = length(ant.vec)),
                       Lambda_2017 = rep(NA, length.out = length(ant.vec)),
                       Lambda_2018 = rep(NA, length.out = length(ant.vec)))
# Eunicea
flex.sens <- data.frame(U1 = flex.vec,
                        Lambda_2014 = rep(NA, length.out = length(flex.vec)),
                        Lambda_2015 = rep(NA, length.out = length(flex.vec)),
                        Lambda_2016 = rep(NA, length.out = length(flex.vec)),
                        Lambda_2017 = rep(NA, length.out = length(flex.vec)),
                        Lambda_2018 = rep(NA, length.out = length(flex.vec)))
# Gorgonia
gorg.sens <- data.frame(U1 = gorg.vec,
                        Lambda_2013 = rep(NA, length.out = length(gorg.vec)),
                        Lambda_2014 = rep(NA, length.out = length(gorg.vec)),
                        Lambda_2015 = rep(NA, length.out = length(gorg.vec)),
                        Lambda_2016 = rep(NA, length.out = length(gorg.vec)),
                        Lambda_2017 = rep(NA, length.out = length(gorg.vec)),
                        Lambda_2018 = rep(NA, length.out = length(gorg.vec)))



# Work through each U1 value for each species generating new lambdaP values to explore the sensitivity of lambdaP to the ceiling threshold.
for (ii in 1:length(ant.vec)) {
  # Progress read out 
  print(ii)
  
  # During each iteration the only parameter to be changed is U1. Everything else will be kept constant.
  
  # store species specific parameters to prevent overwriting the originals.
  m.par.ant <- ant.params
  m.par.flex <- flex.params
  m.par.gorg <- gorg.params
  
  # reassign U1
  m.par.ant['U1',] <- ant.sens$U1[ii]
  m.par.flex['U1',] <- flex.sens$U1[ii]
  m.par.gorg['U1',] <- gorg.sens$U1[ii]
  
  # Generate species-specific IPMs
  ## Antillogorgia
  Ant1 <- mk_P_ceiling(m = m, m.par = m.par.ant[,1], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',1])
  Ant2 <- mk_P_ceiling(m = m, m.par = m.par.ant[,2], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',2])
  Ant3 <- mk_P_ceiling(m = m, m.par = m.par.ant[,3], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',3])
  Ant4 <- mk_P_ceiling(m = m, m.par = m.par.ant[,4], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',4])
  Ant5 <- mk_P_ceiling(m = m, m.par = m.par.ant[,5], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',5])
  Ant6 <- mk_P_ceiling(m = m, m.par = m.par.ant[,6], L = L.list[1], U = U.list[1], U1 = m.par.ant['U1',6])
  ## Eunicea
  flex1 <- mk_P_ceiling(m = m, m.par = m.par.flex[,1], L = L.list[2], U = U.list[2], U1 = m.par.flex['U1',1])
  flex2 <- mk_P_ceiling(m = m, m.par = m.par.flex[,2], L = L.list[2], U = U.list[2], U1 = m.par.flex['U1',2])
  flex3 <- mk_P_ceiling(m = m, m.par = m.par.flex[,3], L = L.list[2], U = U.list[2], U1 = m.par.flex['U1',3])
  flex4 <- mk_P_ceiling(m = m, m.par = m.par.flex[,4], L = L.list[2], U = U.list[2], U1 = m.par.flex['U1',4])
  flex5 <- mk_P_ceiling(m = m, m.par = m.par.flex[,5], L = L.list[2], U = U.list[2], U1 = m.par.flex['U1',5])
  ## Gorgonia
  gorg1 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,1], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',1])
  gorg2 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,2], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',2])
  gorg3 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,3], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',3])
  gorg4 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,4], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',4])
  gorg5 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,5], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',5])
  gorg6 <- mk_P_ceiling(m = m, m.par = m.par.gorg[,6], L = L.list[3], U = U.list[3], U1 = m.par.gorg['U1',6])

  # Calculate and store lambda
  # Antillogorgia
  ant.sens$Lambda_2013[ii] <- Re(eigen(Ant1$P)$values[1])
  ant.sens$Lambda_2014[ii] <- Re(eigen(Ant2$P)$values[1])
  ant.sens$Lambda_2015[ii] <- Re(eigen(Ant3$P)$values[1])
  ant.sens$Lambda_2016[ii] <- Re(eigen(Ant4$P)$values[1])
  ant.sens$Lambda_2017[ii] <- Re(eigen(Ant5$P)$values[1])
  ant.sens$Lambda_2018[ii] <- Re(eigen(Ant6$P)$values[1])
  # Eunicea
  flex.sens$Lambda_2014[ii] <- Re(eigen(flex1$P)$values[1])
  flex.sens$Lambda_2015[ii] <- Re(eigen(flex2$P)$values[1])
  flex.sens$Lambda_2016[ii] <- Re(eigen(flex3$P)$values[1])
  flex.sens$Lambda_2017[ii] <- Re(eigen(flex4$P)$values[1])
  flex.sens$Lambda_2018[ii] <- Re(eigen(flex5$P)$values[1]) 
  # Gorgonia
  gorg.sens$Lambda_2013[ii] <- Re(eigen(gorg1$P)$values[1])
  gorg.sens$Lambda_2014[ii] <- Re(eigen(gorg2$P)$values[1])
  gorg.sens$Lambda_2015[ii] <- Re(eigen(gorg3$P)$values[1])
  gorg.sens$Lambda_2016[ii] <- Re(eigen(gorg4$P)$values[1])
  gorg.sens$Lambda_2017[ii] <- Re(eigen(gorg5$P)$values[1])
  gorg.sens$Lambda_2018[ii] <- Re(eigen(gorg6$P)$values[1]) 
  
}

#####################################
# STEP 3: Quantify and display sensitivities
#####################################

# Reformat data from processing
ant.data <- ant.sens %>%
  gather(Year, Lambda, -U1)
flex.data <- flex.sens %>%
  gather(Year, Lambda, -U1)
gorg.data <- gorg.sens %>%
  gather(Year, Lambda, -U1)

# use regression analysis to determine the relationship between lambda and the positioning of the ceiling size thresholds.
ant.mod <- glm(Lambda ~ U1 * Year, data = ant.data)
flex.mod <- glm(Lambda ~ U1 * Year, data = flex.data)
gorg.mod <- glm(Lambda ~ U1 * Year, data = gorg.data)
# Extract predicted lambda estimates
ant.mod.lines <- ggpredict(ant.mod, terms = c("U1[10.26278:11.46278]", "Year"))
ant.mod.lines <- ant.mod.lines %>% group_split(group)
flex.mod.lines <- ggpredict(flex.mod, terms = c("U1[9.74408:10.94408]", "Year"))
flex.mod.lines <- flex.mod.lines %>% group_split(group)
gorg.mod.lines <- ggpredict(gorg.mod, terms = c("U1[9.09536:10.29536]", "Year"))
gorg.mod.lines <- gorg.mod.lines %>% group_split(group)

# Plot iterative lambda estimates, and the regression predictions
ggplot(data = ant.data) +
  #geom_point(aes(x = U1, y = Lambda, col = Year)) +
  geom_line(aes(x = U1, y = Lambda, col = Year), size = 0.9) +
  scale_colour_manual(values = c("#0000FF", "#708090", "#B0C4DE", "#000000", "#A9A9A9", '#DCDCDC')) +
  geom_vline(aes(xintercept = ant.params['U1',][1]), lty = 'dotted', size = 1) +
  coord_cartesian(ylim = c(0.2, 1), 
                  xlim = c(10.26, 11.46),
                  expand = FALSE) +
  xlab('') + 
  ylab('') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.position='none',
                     legend.box.background = element_rect(colour = "black"))
ggplot(data = flex.data) +
  #geom_point(aes(x = U1, y = Lambda, col = Year)) +
  geom_line(aes(x = U1, y = Lambda, col = Year), size = 0.9) +
  scale_colour_manual(values = c("#708090", "#B0C4DE", "#000000", "#A9A9A9", '#DCDCDC')) +
  geom_vline(aes(xintercept = flex.params['U1',][1]), lty = 'dotted', size = 1) +
  coord_cartesian(ylim = c(0.55, 0.90), 
                  xlim = c(9.80,  10.94),
                  expand = FALSE) +
  xlab('') + 
  ylab('') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.position='none',
                     legend.box.background = element_rect(colour = "black"))
ggplot(data = gorg.data) +
  #geom_point(aes(x = U1, y = Lambda, col = Year)) +
  geom_line(aes(x = U1, y = Lambda, col = Year), size = 0.9) +
  scale_colour_manual(values = c("#0000FF", "#708090", "#B0C4DE", "#000000", "#A9A9A9", '#DCDCDC')) +
  geom_vline(aes(xintercept = gorg.params['U1',][1]), lty = 'dotted', size = 1) +
  coord_cartesian(ylim = c(0.55, 0.95), 
                  xlim = c(9.09, 10.29),
                  expand = FALSE) +
  xlab('') + 
  ylab('') +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                     legend.position='none',
                     legend.box.background = element_rect(colour = "black"))


# Calculate the numerical sensitivity of the models to the positioning of the size threshold.
# Sensitivity = delta(lambda)/delta(U1)
# Antillogorgia
(ant.sens[122,'Lambda_2013']-ant.sens[1,'Lambda_2013'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
(ant.sens[122,'Lambda_2014']-ant.sens[1,'Lambda_2014'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
(ant.sens[122,'Lambda_2015']-ant.sens[1,'Lambda_2015'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
(ant.sens[122,'Lambda_2016']-ant.sens[1,'Lambda_2016'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
(ant.sens[122,'Lambda_2017']-ant.sens[1,'Lambda_2017'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
(ant.sens[122,'Lambda_2018']-ant.sens[1,'Lambda_2018'])/(ant.sens[122,'U1']-ant.sens[1,'U1'])
# Eunicea
(flex.sens[122,'Lambda_2014']-flex.sens[1,'Lambda_2014'])/(flex.sens[122,'U1']-flex.sens[1,'U1'])
(flex.sens[122,'Lambda_2015']-flex.sens[1,'Lambda_2015'])/(flex.sens[122,'U1']-flex.sens[1,'U1'])
(flex.sens[122,'Lambda_2016']-flex.sens[1,'Lambda_2016'])/(flex.sens[122,'U1']-flex.sens[1,'U1'])
(flex.sens[122,'Lambda_2017']-flex.sens[1,'Lambda_2017'])/(flex.sens[122,'U1']-flex.sens[1,'U1'])
(flex.sens[122,'Lambda_2018']-flex.sens[1,'Lambda_2018'])/(flex.sens[122,'U1']-flex.sens[1,'U1'])
# Gorgonia
(gorg.sens[122,'Lambda_2013']-gorg.sens[1,'Lambda_2013'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])
(gorg.sens[122,'Lambda_2014']-gorg.sens[1,'Lambda_2014'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])
(gorg.sens[122,'Lambda_2015']-gorg.sens[1,'Lambda_2015'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])
(gorg.sens[122,'Lambda_2016']-gorg.sens[1,'Lambda_2016'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])
(gorg.sens[122,'Lambda_2017']-gorg.sens[1,'Lambda_2017'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])
(gorg.sens[122,'Lambda_2018']-gorg.sens[1,'Lambda_2018'])/(gorg.sens[122,'U1']-gorg.sens[1,'U1'])

######################### End of code #############################