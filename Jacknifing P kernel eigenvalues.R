# This script follows on from the construction of the individual species IPM's and will Jacknife through each population to determine variance for the 
# calculated P kernel eigenvalues.

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
#-----------------------------------

# load required packages
library(ggplot2)
library(Rmisc)
library(reshape2)

# set the random seed
set.seed(46789)

###############################################
# STEP 1: Organise required data and storage locations.
###############################################

# re-store species dataframe so that they are not affected during Jackknifing.
G.data <- Gorg_sp 
A.data <- A_american
E.data <- E_flexuosa
Recruits.ant <- Recruits[which(Recruits$Genus == "Antillogorgia"),]
Recruits.flex <- Recruits[which(Recruits$Genus == "Eunicea"),]
Recruits.gorg <- Recruits[which(Recruits$Genus == "Gorgonia"),]
Adults.ant <- Adults[which(Adults$Genus == "Antillogorgia"),]
Adults.flex <- Adults[which(Adults$Genus == "Eunicea"),]
Adults.gorg <- Adults[which(Adults$Genus == "Gorgonia"),]

# lambda storage
lambdas.sc <- matrix(NA, 1000, 17) #there is a list of lambdas for each species in each year
# labelling each column will make it easier later
colnames(lambdas.sc) <- c("Ant.13","Ant.14","Ant.15","Ant.16","Ant.17","Ant.18",
                          "Flex.14","Flex.15","Flex.16","Flex.17","Flex.18",
                          "Gorg.13","Gorg.14","Gorg.15","Gorg.16","Gorg.17","Gorg.18")

###########################################
# STEP 2: Run Jacknife loop
###########################################

for (x in 1:1000) {
  # progress read out
  print(x)
        
  # determine the new Jackknife sample
  # the use of try should mean the loop will return NA's rather than breaking for errors but means that = signs can't be used.
  try(G.data.1 <- G.data[sample(nrow(G.data), 0.95*dim(G.data)[1], replace = F),]) #this randomly takes a 95% subset of each species dataset.
  try(A.data.1 <- A.data[sample(nrow(A.data), 0.95*dim(A.data)[1], replace = F),])
  try(E.data.1 <- E.data[sample(nrow(E.data), 0.95*dim(E.data)[1], replace = F),])
  # and again for the full community data set (this is only used here to determine the max and min size boundaries of the populations)
  try(Recruit.ant.1 <- Recruits.ant[sample(nrow(Recruits.ant), 0.95*dim(Recruits.ant)[1], replace = F),]) 
  try(Recruit.flex.1 <- Recruits.flex[sample(nrow(Recruits.flex), 0.95*dim(Recruits.flex)[1], replace = F),])
  try(Recruit.gorg.1 <- Recruits.gorg[sample(nrow(Recruits.gorg), 0.95*dim(Recruits.gorg)[1], replace = F),])
  try(Adult.ant.1 <- Adults.ant[sample(nrow(Adults.ant), 0.95*dim(Adults.ant)[1], replace = F),]) 
  try(Adult.flex.1 <- Adults.flex[sample(nrow(Adults.flex), 0.95*dim(Adults.flex)[1], replace = F),])
  try(Adult.gorg.1 <- Adults.gorg[sample(nrow(Adults.gorg), 0.95*dim(Adults.gorg)[1], replace = F),])
  
  # Re-calculate each of the vital rates of growth and survival forcing the data to follow the same regression formats as the orginal models
  #### Growth ####
      ## A.american
      try(jAmer_grow1 <-lm(Size.t1 ~ Size.t * Year.t, data = A.data.1)) 
      ## E. flexuosa
      try(jflex_grow1 <-lm(Size.t1 ~ Size.t * Year.t, data = E.data.1)) 
      ## Gorgonia spp.
      try(jgorg_grow1 <-lm(Size.t1 ~ Size.t * Year.t, data = G.data.1)) 
        
  #### Variance in Growth ####
      ## A.american
      try(jAmer.grow.sd <- glm(abs(resid(jAmer_grow1))~jAmer_grow1$model$Size.t, family = Gamma(link = "log"))) 
      ## E.flexuosa
      try(jflex.grow.sd <- glm(abs(resid(jflex_grow1))~jflex_grow1$model$Size.t, family = Gamma(link = "log"))) 
      ## Gorgonia
      try(jgorg.grow.sd <- glm(abs(resid(jgorg_grow1))~jgorg_grow1$model$Size.t, family = Gamma(link = "log"))) 
      
  #### Survival ####
      ## A. american
      try(jant_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = A.data.1))
      ## E. flexuosa
      try(jflex_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = E.data.1))
      ## Gorgonia
      try(jgorg_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = G.data.1))
       
  #### Store model coefficients ####
      # blank matrices to store the yearly parameters for each species separately
      try(jgorg.params <- matrix(NA, 8, 6))
      try(jflex.params <- matrix(NA, 8, 5))
      try(jant.params <- matrix(NA, 8, 6))
      
      # Fill each matrix.
      # survival
      try(for(i in 1:6){
        if (i == 1){
          jant.params[1,i] <- coefficients(jant_surv1)[1]
          jant.params[2,i] <- coefficients(jant_surv1)[2]
          jant.params[3,i] <- sigma(jant_surv1)
        } else {
          jant.params[1,i] <- coefficients(jant_surv1)[1] + coefficients(jant_surv1)[i+1] 
          jant.params[2,i] <- coefficients(jant_surv1)[2] + coefficients(jant_surv1)[i+6]
          jant.params[3,i] <- sigma(jant_surv1)
        }  
      })
      
      try(for(i in 1:6){
        if (i == 1){
          jgorg.params[1,i] <- coefficients(jgorg_surv1)[1]
          jgorg.params[2,i] <- coefficients(jgorg_surv1)[2]
          jgorg.params[3,i] <- sigma(jgorg_surv1)
        } else {
          jgorg.params[1,i] <- coefficients(jgorg_surv1)[1] + coefficients(jgorg_surv1)[i+1] 
          jgorg.params[2,i] <- coefficients(jgorg_surv1)[2] + coefficients(jgorg_surv1)[i+6]
          jgorg.params[3,i] <- sigma(jgorg_surv1)
        }  
      })
      
      try(for(i in 1:5){
        if (i == 1){
          jflex.params[1,i] <- coefficients(jflex_surv1)[1]
          jflex.params[2,i] <- coefficients(jflex_surv1)[2]
          jflex.params[3,i] <- sigma(jflex_surv1)
        } else {
          jflex.params[1,i] <- coefficients(jflex_surv1)[1] + coefficients(jflex_surv1)[i+1] 
          jflex.params[2,i] <- coefficients(jflex_surv1)[2] + coefficients(jflex_surv1)[i+5]
          jflex.params[3,i] <- sigma(jflex_surv1)
        }  
      })
      
      # growth
      try(for(i in 1:6){
        if (i == 1){
          jant.params[4,i] <- coefficients(jAmer_grow1)[1]
          jant.params[5,i] <- coefficients(jAmer_grow1)[2]
        } else {
          jant.params[4,i] <- coefficients(jAmer_grow1)[1] + coefficients(jAmer_grow1)[i+1] 
          jant.params[5,i] <- coefficients(jAmer_grow1)[2] + coefficients(jAmer_grow1)[i+6]
        }  
      })
      
      try(for(i in 1:6){
        if (i == 1){
          jgorg.params[4,i] <- coefficients(jgorg_grow1)[1]
          jgorg.params[5,i] <- coefficients(jgorg_grow1)[2]
        } else {
          jgorg.params[4,i] <- coefficients(jgorg_grow1)[1] + coefficients(jgorg_grow1)[i+1] 
          jgorg.params[5,i] <- coefficients(jgorg_grow1)[2] + coefficients(jgorg_grow1)[i+6]
        }  
      })
      
      try(for(i in 1:5){
        if (i == 1){
          jflex.params[4,i] <- coefficients(jflex_grow1)[1]
          jflex.params[5,i] <- coefficients(jflex_grow1)[2]
        } else {
          jflex.params[4,i] <- coefficients(jflex_grow1)[1] + coefficients(jflex_grow1)[i+1] 
          jflex.params[5,i] <- coefficients(jflex_grow1)[2] + coefficients(jflex_grow1)[i+5]
        }  
      })
      
      # growth variability 
      try(for(i in 1:6){
        jant.params[6,i] <- coefficients(jAmer.grow.sd)[1]
        jant.params[7,i] <- coefficients(jAmer.grow.sd)[2]
        jgorg.params[6,i] <- coefficients(jgorg.grow.sd)[1]
        jgorg.params[7,i] <- coefficients(jgorg.grow.sd)[2]
      }) 
      
      try(for(i in 1:5){
        jflex.params[6,i] <- coefficients(jflex.grow.sd)[1]
        jflex.params[7,i] <- coefficients(jflex.grow.sd)[2]
      }) 
      
      # Now if tidy up the parameter matrices
      try(rownames(jant.params) <- rownames(jgorg.params) <- rownames(jflex.params) <- c(
        #survival
        "surv.int","surv.slope", "surv.sd",
        #growth
        "grow.int", "grow.slope",
        #growth variability
        "grow.sd.int","grow.sd.slope",
        # Size ceiling
        "U1"))
      
      try(colnames(jant.params) <- colnames(jgorg.params) <- c("2013","2014","2015","2016","2017","2018"))
      try(colnames(jflex.params) <- c("2014","2015","2016","2017","2018"))
      
   #### These jackknife models will use the same model functions generated for the original models.   
   #### Define final model details
      
      # define model parameters
      try(jL.list <- c(0.9*min(c(A.data.1$Size.t, A.data.1$Size.t1, Recruit.ant.1$Height..cm.), na.rm = T), 
                  0.9*min(c(E.data.1$Size.t, E.data.1$Size.t1, Recruit.flex.1$Height..cm.), na.rm = T),
                  0.9*min(c(G.data.1$Size.t, G.data.1$Size.t1, Recruit.gorg.1$Height..cm.), na.rm = T)))
      try(jU.list <- c(1.1*max(c(A.data.1$Size.t, A.data.1$Size.t1, Adult.ant.1$Height..cm.), na.rm = T), 
                  1.1*max(c(E.data.1$Size.t, E.data.1$Size.t1, Adult.flex.1$Height..cm.), na.rm = T),
                  1.1*max(c(G.data.1$Size.t, G.data.1$Size.t1, Adult.gorg.1$Height..cm.), na.rm = T)))
  
      # define ceiling parameters
      try(jant.params[8,] <- max(c(A.data.1$Size.t, A.data.1$Size.t1, Adult.ant.1$Height..cm.), na.rm = T))
      try(jflex.params[8,] <- max(c(E.data.1$Size.t, E.data.1$Size.t1, Adult.flex.1$Height..cm.), na.rm = T))
      try(jgorg.params[8,] <- max(c(G.data.1$Size.t, G.data.1$Size.t1, Adult.gorg.1$Height..cm.), na.rm = T))
      
   #### The value of m used will stay the same, as will the initial Nt values (as Im only interested here in the P values)
      
   #### Build the IPMs ####
      ## Antillogorgia
      try(jAnt.2013 <- mk_P_ceiling(m = m, m.par = jant.params[,1], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,1]))
      try(jAnt.2014 <- mk_P_ceiling(m = m, m.par = jant.params[,2], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,2]))
      try(jAnt.2015 <- mk_P_ceiling(m = m, m.par = jant.params[,3], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,3]))
      try(jAnt.2016 <- mk_P_ceiling(m = m, m.par = jant.params[,4], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,4]))
      try(jAnt.2017 <- mk_P_ceiling(m = m, m.par = jant.params[,5], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,5]))
      try(jAnt.2018 <- mk_P_ceiling(m = m, m.par = jant.params[,6], L = jL.list[1], U = jU.list[1], U1 = jant.params[8,6]))
      
      ## Eunicea
      try(jflex.2014 <- mk_P_ceiling(m = m, m.par = jflex.params[,1], L = jL.list[2], U = jU.list[2], U1 = jflex.params[8,1]))
      try(jflex.2015 <- mk_P_ceiling(m = m, m.par = jflex.params[,2], L = jL.list[2], U = jU.list[2], U1 = jflex.params[8,2]))
      try(jflex.2016 <- mk_P_ceiling(m = m, m.par = jflex.params[,3], L = jL.list[2], U = jU.list[2], U1 = jflex.params[8,3]))
      try(jflex.2017 <- mk_P_ceiling(m = m, m.par = jflex.params[,4], L = jL.list[2], U = jU.list[2], U1 = jflex.params[8,4]))
      try(jflex.2018 <- mk_P_ceiling(m = m, m.par = jflex.params[,5], L = jL.list[2], U = jU.list[2], U1 = jflex.params[8,5]))
      
      ## Gorgonia
      try(jgorg.2013 <- mk_P_ceiling(m = m, m.par = jgorg.params[,1], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,1]))
      try(jgorg.2014 <- mk_P_ceiling(m = m, m.par = jgorg.params[,2], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,2]))
      try(jgorg.2015 <- mk_P_ceiling(m = m, m.par = jgorg.params[,3], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,3]))
      try(jgorg.2016 <- mk_P_ceiling(m = m, m.par = jgorg.params[,4], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,4]))
      try(jgorg.2017 <- mk_P_ceiling(m = m, m.par = jgorg.params[,5], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,5]))
      try(jgorg.2018 <- mk_P_ceiling(m = m, m.par = jgorg.params[,6], L = jL.list[3], U = jU.list[3], U1 = jgorg.params[8,6]))
 
   #### Calculate and store lambda ####     
      ## Antilogorgia
      try(lambdas.sc[x,1] <- Re(eigen(jAnt.2013$P)$values[1]))
      try(lambdas.sc[x,2] <- Re(eigen(jAnt.2014$P)$values[1]))
      try(lambdas.sc[x,3] <- Re(eigen(jAnt.2015$P)$values[1]))
      try(lambdas.sc[x,4] <- Re(eigen(jAnt.2016$P)$values[1]))
      try(lambdas.sc[x,5] <- Re(eigen(jAnt.2017$P)$values[1])) 
      try(lambdas.sc[x,6] <- Re(eigen(jAnt.2018$P)$values[1]))
      
      ## Eunicea
      try(lambdas.sc[x,7] <- Re(eigen(jflex.2014$P)$values[1]))
      try(lambdas.sc[x,8] <- Re(eigen(jflex.2015$P)$values[1]))
      try(lambdas.sc[x,9] <- Re(eigen(jflex.2016$P)$values[1]))
      try(lambdas.sc[x,10] <- Re(eigen(jflex.2017$P)$values[1]))
      try(lambdas.sc[x,11] <- Re(eigen(jflex.2018$P)$values[1]))
        
      ## Gorgonia
      try(lambdas.sc[x,12] <- Re(eigen(jgorg.2013$P)$values[1]))
      try(lambdas.sc[x,13] <- Re(eigen(jgorg.2014$P)$values[1]))
      try(lambdas.sc[x,14] <- Re(eigen(jgorg.2015$P)$values[1]))
      try(lambdas.sc[x,15] <- Re(eigen(jgorg.2016$P)$values[1]))
      try(lambdas.sc[x,16] <- Re(eigen(jgorg.2017$P)$values[1]))
      try(lambdas.sc[x,17] <- Re(eigen(jgorg.2018$P)$values[1]))
  
} # End of Jackknife loop

###########################################
# STEP 3: Determine variance and test for significant differences
###########################################

# Now I have a dataframe of different lambda values I can determine the variance of each sample
# Standard deviation
sd(lambdas.sc[,1], na.rm = T)
sd(lambdas.sc[,2], na.rm = T)
sd(lambdas.sc[,3], na.rm = T)
sd(lambdas.sc[,4], na.rm = T)
sd(lambdas.sc[,5], na.rm = T)
sd(lambdas.sc[,6], na.rm = T)
sd(lambdas.sc[,7], na.rm = T)
sd(lambdas.sc[,8], na.rm = T)
sd(lambdas.sc[,9], na.rm = T)
sd(lambdas.sc[,10], na.rm = T)
sd(lambdas.sc[,11], na.rm = T)
sd(lambdas.sc[,12], na.rm = T)
sd(lambdas.sc[,13], na.rm = T)
sd(lambdas.sc[,14], na.rm = T)
sd(lambdas.sc[,15], na.rm = T)
sd(lambdas.sc[,16], na.rm = T)
sd(lambdas.sc[,17], na.rm = T)
# Confidence intervals
CI(na.omit(lambdas.sc[,1]))
CI(na.omit(lambdas.sc[,2]))
CI(na.omit(lambdas.sc[,3]))
CI(na.omit(lambdas.sc[,4]))
CI(na.omit(lambdas.sc[,5]))
CI(na.omit(lambdas.sc[,6]))
CI(na.omit(lambdas.sc[,7]))
CI(na.omit(lambdas.sc[,8]))
CI(na.omit(lambdas.sc[,9]))
CI(na.omit(lambdas.sc[,10]))
CI(na.omit(lambdas.sc[,11]))
CI(na.omit(lambdas.sc[,12]))
CI(na.omit(lambdas.sc[,13]))
CI(na.omit(lambdas.sc[,14]))
CI(na.omit(lambdas.sc[,15]))
CI(na.omit(lambdas.sc[,16]))
CI(na.omit(lambdas.sc[,17]))

# Now I can run a two-way ANOVA to find if there is a difference between species
# for this I need to rearrange the lambda dataframe
# first stack all the lambdas together
lambdas.sc.1 <- melt(lambdas.sc)
# rename variables
colnames(lambdas.sc.1) <- c("Species", "Year", "LambdaP")
lambdas.sc.1$Species <- c(rep("A_Americana", 6000),rep("E_flexuosa",5000),rep("Gorgonia_sp", 6000))
lambdas.sc.1$Year <- c(rep(2013, 1000),rep(2014, 1000),rep(2015, 1000),rep(2016, 1000),rep(2017, 1000),rep(2018, 1000),
                     rep(2014, 1000),rep(2015, 1000),rep(2016, 1000),rep(2017, 1000),rep(2018, 1000),
                     rep(2013, 1000),rep(2014, 1000),rep(2015, 1000),rep(2016, 1000),rep(2017, 1000),rep(2018, 1000))

# Now run the ANOVA
model <- lm(LambdaP ~ Species + Year + Species:Year, data = lambdas.sc.1)
# and determine the significance.
car::Anova(model, type = "II")
# since there is a significant interaction we have to evaluate everything at the same time.
# this can be done using an interaction plot
set_graph_pars("panel1")
# store the mean and sd of each species
interaction.data <- data.frame(Year = c(2013,2014,2015,2016,2017,2018,
                                        2014,2015,2016,2017,2018,
                                        2013,2014,2015,2016,2017,2018),
                               Species = c(rep("Antilogorgia", length.out = 6),
                                           rep("Eunicea", length.out = 5),
                                           rep("Gorgornia", length.out = 6)),
                               Mean = c(mean(lambdas.sc[,1], na.rm = T), mean(lambdas.sc[,2], na.rm = T), mean(lambdas.sc[,3], na.rm = T), mean(lambdas.sc[,4], na.rm = T), mean(lambdas.sc[,5], na.rm = T), mean(lambdas.sc[,6], na.rm = T),
                                        mean(lambdas.sc[,7], na.rm = T), mean(lambdas.sc[,8], na.rm = T), mean(lambdas.sc[,9], na.rm = T), mean(lambdas.sc[,10], na.rm = T), mean(lambdas.sc[,11], na.rm = T),
                                        mean(lambdas.sc[,12], na.rm = T), mean(lambdas.sc[,13], na.rm = T), mean(lambdas.sc[,14], na.rm = T), mean(lambdas.sc[,15], na.rm = T), mean(lambdas.sc[,16], na.rm = T), mean(lambdas.sc[,17], na.rm = T)),
                               SD = c(sd(lambdas.sc[,1], na.rm = T), sd(lambdas.sc[,2], na.rm = T), sd(lambdas.sc[,3], na.rm = T), sd(lambdas.sc[,4], na.rm = T), sd(lambdas.sc[,5], na.rm = T), sd(lambdas.sc[,6], na.rm = T),
                                      sd(lambdas.sc[,7], na.rm = T), sd(lambdas.sc[,8], na.rm = T), sd(lambdas.sc[,9], na.rm = T), sd(lambdas.sc[,10], na.rm = T), sd(lambdas.sc[,11], na.rm = T),
                                      sd(lambdas.sc[,12], na.rm = T), sd(lambdas.sc[,13], na.rm = T), sd(lambdas.sc[,14], na.rm = T), sd(lambdas.sc[,15], na.rm = T), sd(lambdas.sc[,16], na.rm = T), sd(lambdas.sc[,17], na.rm = T)),
                               CIlow = c(CI(na.omit(lambdas.sc[,1]))[3], CI(na.omit(lambdas.sc[,2]))[3], CI(na.omit(lambdas.sc[,3]))[3], CI(na.omit(lambdas.sc[,4]))[3], CI(na.omit(lambdas.sc[,5]))[3], CI(na.omit(lambdas.sc[,6]))[3], 
                                         CI(na.omit(lambdas.sc[,7]))[3], CI(na.omit(lambdas.sc[,8]))[3], CI(na.omit(lambdas.sc[,9]))[3], CI(na.omit(lambdas.sc[,10]))[3], CI(na.omit(lambdas.sc[,11]))[3], CI(na.omit(lambdas.sc[,12]))[3], 
                                         CI(na.omit(lambdas.sc[,13]))[3], CI(na.omit(lambdas.sc[,14]))[3], CI(na.omit(lambdas.sc[,15]))[3], CI(na.omit(lambdas.sc[,16]))[3], CI(na.omit(lambdas.sc[,17]))[3]),
                               CIhigh = c(CI(na.omit(lambdas.sc[,1]))[1], CI(na.omit(lambdas.sc[,2]))[1], CI(na.omit(lambdas.sc[,3]))[1], CI(na.omit(lambdas.sc[,4]))[1], CI(na.omit(lambdas.sc[,5]))[1], CI(na.omit(lambdas.sc[,6]))[1], 
                                          CI(na.omit(lambdas.sc[,7]))[1], CI(na.omit(lambdas.sc[,8]))[1], CI(na.omit(lambdas.sc[,9]))[1], CI(na.omit(lambdas.sc[,10]))[1], CI(na.omit(lambdas.sc[,11]))[1], CI(na.omit(lambdas.sc[,12]))[1], 
                                          CI(na.omit(lambdas.sc[,13]))[1], CI(na.omit(lambdas.sc[,14]))[1], CI(na.omit(lambdas.sc[,15]))[1], CI(na.omit(lambdas.sc[,16]))[1], CI(na.omit(lambdas.sc[,17]))[1]))

# create the interaction plot
interaction <- ggplot(interaction.data, aes(x = Year, y = Mean, col = Species, group = Species), xlab(NULL), ylab(NULL)) + 
  labs(y = "", # add axis labels
       x = "") + 
  coord_cartesian(ylim = c(0, 1)) 
interaction.plot <- interaction + geom_line(aes(linetype = Species), size = 0.6) + 
                                  scale_linetype_manual(values = c("dashed", "dashed", "dashed")) +
                                  geom_point(aes(col = Species), size = 3) + 
                                  scale_colour_manual(values = c("#D55E00","#0072B2","#000000")) +
                                  geom_errorbar(aes(ymax = CIhigh, ymin = CIlow), width = 0.2) + 
                                  theme_bw() + theme(panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                                                     legend.justification=c(0,0), legend.position=c(0.1,0.1), #put legend in plot
                                                     legend.background = element_blank(),
                                                     legend.box.background = element_rect(colour = "black"))

######################### End of code #############################