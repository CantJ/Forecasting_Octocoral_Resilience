# This script is for working with the Caribbean octocoral demography data to construct density dependent IPMs.
# This script goes through the calculation of vital rates (growth survival and recruitment), model parameterisation and
# determination of adult viability (This is done simultaneously for all species).

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
#-----------------------------------

# Clear work space memory
rm(list=ls(all=TRUE))

#set working directory and load in required packages.
library(lme4)
library(ggeffects)
library(emmeans)
library(splines)
library(fitdistrplus)
library(FSA)
library(energy)
library(gamlss)
library(gamlss.dist)
library(fields)
library(dplyr)

# set working directory
setwd("DATA DIRECTORY PATHWAY")

#load in the data set
soft_corals <- read.csv("IPM_DATA_FINAL", stringsAsFactors = TRUE)
# this is the size and survival data collected from tagged adult colonies between 2013 and 2019 annually.

###############################################
# STEP 1: Sort the data.
###############################################

# Define variable formats
soft_corals$ID <- as.factor(soft_corals$ID)
soft_corals$Year.t <- as.factor(soft_corals$Year.t)
soft_corals$Surv <- as.numeric(soft_corals$Surv)

# Determine nessecary size data transformations.
hist(soft_corals$Size.t) #little bit skewed to the left
hist(soft_corals$Size.t1) #same here - the data does need transforming
hist(log(soft_corals$Size.t))
hist(log(soft_corals$Size.t1)) #this skews far to much to the right
hist(log(soft_corals$Size.t+1))
hist(log(soft_corals$Size.t1+1)) #less so but still not normal
hist(sqrt(soft_corals$Size.t))
hist(sqrt(soft_corals$Size.t1)) #square root transformation works the best.

#and implement transformation
soft_corals$Size.t <- sqrt(soft_corals$Size.t)
soft_corals$Size.t1 <- sqrt(soft_corals$Size.t1) # size data (lots of little, little of large) is easier to analyse on a transformed scale.

#and lastly split the data into the component species
Gorg_sp <- subset(soft_corals, Species == "GV")
A_american <- subset(soft_corals, Species == "AA")
E_flexuosa <- subset(soft_corals, Species == "EFL")
# Across the analyses sites will remained pooled and only considered as a random effect

# The data is now ready for analysis!!

###############################################
# STEP 2: Determine annual variance in survival and growth.
###############################################

########### Growth ###########
## A.american
plot(Size.t1 ~ Size.t, col = Year.t, data = A_american)

# for these models time will be used as a fixed effect. But I will run multiple regression formats to determine whether the model is non-linear, polynomial or contains random effects.
Amer_grow1 <- lm(Size.t1 ~ Size.t * Year.t, data = A_american) #fixed effect model
Amer_grow2 <- lmer(Size.t1 ~ Size.t * Year.t + (Size.t|Site), data = A_american) #includes the random effect of site (we are pooling colonies tagged in East Cabritte and Europa)
# this model has a singular fit
# which model fits best?
AIC(Amer_grow1,Amer_grow2) #the model excluding the random effect of site is the best fit.
BIC(Amer_grow1,Amer_grow2) #this is agreed by the BIC selection (the BIC penalizes model complexity more than an AIC)

# but is the relationship linear?
Amer_grow3 <- glm(Size.t1 ~ Size.t + I(Size.t^2) * Year.t, data = A_american)
AIC(Amer_grow1, Amer_grow3) # linear fits best.
BIC(Amer_grow1, Amer_grow3) # the BIC score prefers the less complex linear regression.

Amer1 <- ggpredict(Amer_grow1, terms = c("Size.t[0:11]", "Year.t")) #predict values for a 'line of best fit'
Amer1 <- split(Amer1, Amer1$group) #split the predicted values by year so they can all be plotted separately.

# Now plot these lines on the growth plot
lines(Amer1$'2013'$x, Amer1$'2013'$predicted, lty = "dashed", col = "blue")
lines(Amer1$'2014'$x, Amer1$'2014'$predicted, lty = "dashed", col = "purple")
lines(Amer1$'2015'$x, Amer1$'2015'$predicted, lty = "dashed", col = "green")
lines(Amer1$'2016'$x, Amer1$'2016'$predicted, lty = "dashed", col = "gray")
lines(Amer1$'2017'$x, Amer1$'2017'$predicted, lty = "solid", col = "red", lwd = 2) # Hurricane year
lines(Amer1$'2018'$x, Amer1$'2018'$predicted, lty = "dashed", col = "black")
abline(0,1,lwd = 3) #this confirms that there is no weird growth occurring!
# visually there is no real difference between model types so I will follow the AICs and so I'll go linear. Interestingly the growth in both 2017 and 2018 are just as bad...
#..perhaps suggests an argument for incorporating a recovery year?

#check the regression outputs for the model
summary(Amer_grow1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(Amer_grow1, test = "Chisq")
emmeans(Amer_grow1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

## E. flexuosa
plot(Size.t1 ~ Size.t, col = Year.t, data = E_flexuosa)

# for these models time will be used as a fixed effect. But I will run multiple regression formats to determine whether the model is non-linear, polynomial or contains random effects.
flex_grow1 <-lm(Size.t1 ~ Size.t * Year.t, data = E_flexuosa) #fixed effect model
flex_grow2 <- lmer(Size.t1 ~ Size.t * Year.t + (Size.t|Site), data = E_flexuosa) #includes the random effect of site (we are pooling colonies tagged in East Cabritte and Europa)
# Another singular fit (i.e not enough repetition of random effects)
# which model fits best?
AIC(flex_grow1, flex_grow2) #the model excluding the random effect of site is the best fit.
BIC(flex_grow1, flex_grow2)

# but is the relationship linear?
flex_grow3 <- glm(Size.t1 ~ Size.t + I(Size.t^2) * Year.t, data = E_flexuosa)
AIC(flex_grow1, flex_grow3) # non-linear fits better according to AIC.
BIC(flex_grow1, flex_grow3) # but linear fits better according to BIC so lets check visually.

flex1 <- ggpredict(flex_grow1, terms = c("Size.t[0:11]", "Year.t")) #predict values for a 'line of best fit'
flex1 <- split(flex1, flex1$group) #split the predicted values by year so they can all be plotted seperately.
flex2 <- ggpredict(flex_grow3, terms = c("Size.t[0:11]", "Year.t"))
flex2 <- split(flex2, flex2$group)

# Now plot these lines on the growth plot
lines(flex1$'2014'$x, flex1$'2014'$predicted, lty = "dashed")
lines(flex1$'2015'$x, flex1$'2015'$predicted, lty = "dashed")
lines(flex1$'2016'$x, flex1$'2016'$predicted, lty = "dashed")
lines(flex1$'2017'$x, flex1$'2017'$predicted, lty = "dashed", col = "red", lwd = 2)
lines(flex1$'2018'$x, flex1$'2018'$predicted, lty = "dashed")
lines(flex2$'2014'$x, flex2$'2014'$predicted, lty = "solid")
lines(flex2$'2015'$x, flex2$'2015'$predicted, lty = "solid")
lines(flex2$'2016'$x, flex2$'2016'$predicted, lty = "solid")
lines(flex2$'2017'$x, flex2$'2017'$predicted, lty = "solid", col = "red", lwd = 2)
lines(flex2$'2018'$x, flex2$'2018'$predicted, lty = "solid")
abline(0,1,lwd = 3) #this confirms that there is no weird growth occurring!
# visually the none linear trends are having a bit of an effect at the larger sizes... resulting in shrinkage 
# However looking at the 2018 data the curve is being fit aggressively to only a few data points. So I will follow the BIC scores and so I'll go linear.
# Interestingly again the growth in both 2017 and 2018 are the worst....was 2018 a recovery year for all species?

#check the regression outputs for the model
summary(flex_grow1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(flex_grow1, test = "Chisq")
emmeans(flex_grow1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

## Gorgonia spp.
plot(Size.t1 ~ Size.t, col = Year.t, data = Gorg_sp)

# for these models time will be used as a fixed effect. But I will run multiple regression formats to determine whether the model is non-linear, polynomial or contains random effects.
gorg_grow1 <-lm(Size.t1 ~ Size.t * Year.t, data = Gorg_sp) #fixed effect model
gorg_grow2 <- lmer(Size.t1 ~ Size.t * Year.t + (Size.t|Site), data = Gorg_sp) #includes the random effect of site (we are pooling colonies tagged in East Cabritte and Europa)
# Is again singular
# which model fits best?
AIC(gorg_grow1, gorg_grow2) #the model excluding the random effect of site is the best fit.
BIC(gorg_grow1, gorg_grow2) #the BIC scores agree.

# but is the relationship linear?
gorg_grow3 <- glm(Size.t1 ~ Size.t + I(Size.t^2) * Year.t, data = Gorg_sp)
AIC(gorg_grow1, gorg_grow3) # non-linear fits better but lets check visually.
BIC(gorg_grow1, gorg_grow3) # linear fits best for the BIC so lets check it!

gorg1 <- ggpredict(gorg_grow1, terms = c("Size.t[0:11]", "Year.t")) #predict values for a 'line of best fit'
gorg1 <- split(gorg1, gorg1$group) #split the predicted values by year so they can all be plotted seperately.
gorg2 <- ggpredict(gorg_grow3, terms = c("Size.t[0:11]", "Year.t"))
gorg2 <- split(gorg2, gorg2$group)

# Now plot these lines on the growth plot
lines(gorg1$'2013'$x, gorg1$'2013'$predicted, lty = "dashed")
lines(gorg1$'2014'$x, gorg1$'2014'$predicted, lty = "dashed")
lines(gorg1$'2015'$x, gorg1$'2015'$predicted, lty = "dashed")
lines(gorg1$'2016'$x, gorg1$'2016'$predicted, lty = "dashed")
lines(gorg1$'2017'$x, gorg1$'2017'$predicted, lty = "dashed", col = "red", lwd = 2)
lines(gorg1$'2018'$x, gorg1$'2018'$predicted, lty = "dashed")
lines(gorg2$'2013'$x, gorg2$'2013'$predicted, lty = "solid")
lines(gorg2$'2014'$x, gorg2$'2014'$predicted, lty = "solid")
lines(gorg2$'2015'$x, gorg2$'2015'$predicted, lty = "solid")
lines(gorg2$'2016'$x, gorg2$'2016'$predicted, lty = "solid")
lines(gorg2$'2017'$x, gorg2$'2017'$predicted, lty = "solid", col = "red", lwd = 2)
lines(gorg2$'2018'$x, gorg2$'2018'$predicted, lty = "solid")
abline(0,1,lwd = 3) #this confirms where the weird growth occurring!
# it appears again that the non-linear model is over fitting curves to few data points. I will therefore stick with the BIC scores and go linear.
# For this species growth was worst during 2013 interestingly

#check the regression outputs for the model
summary(gorg_grow1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(gorg_grow1, test = "Chisq")
emmeans(gorg_grow1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

## Now I can plot the annual growth characteristics for each species.

#Antillagorgia
plot(Size.t1~ Size.t, data = A_american, col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(1,11),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5) #<- size against sizet+1
abline(0,1, lty = "dashed", lwd = 2) #a reference line for the null growth threshold.
lines(Amer1$'2013'$x, Amer1$'2013'$predicted, col = "Green", lwd = 2) #add the regression lines using the model coefficients
lines(Amer1$'2014'$x, Amer1$'2014'$predicted, col = "Black", lwd = 2)
lines(Amer1$'2015'$x, Amer1$'2015'$predicted, col = "blue", lwd = 2)
lines(Amer1$'2016'$x, Amer1$'2016'$predicted, col = "Grey", lwd = 2)
lines(Amer1$'2017'$x, Amer1$'2017'$predicted, col = "red", lwd = 2)
lines(Amer1$'2018'$x, Amer1$'2018'$predicted, col = "purple", lwd = 2)
#legend(6, 4, legend = c("2013","2014", "2015","2016","2017"), fill = c("Green","Black","blue","Grey","red"), bg = NULL, bty = "n")

#Eunicea
plot(Size.t1~ Size.t, data = E_flexuosa, col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(1,11),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5) #<- size against sizet+1
abline(0,1, lty = "dashed", lwd = 2) #a reference line for the null growth threshold.
lines(flex1$'2014'$x, flex1$'2014'$predicted, col = "Black", lwd = 2) #add the regression lines using the model coefficients
lines(flex1$'2015'$x, flex1$'2015'$predicted, col = "blue", lwd = 2)
lines(flex1$'2016'$x, flex1$'2016'$predicted, col = "Grey", lwd = 2)
lines(flex1$'2017'$x, flex1$'2017'$predicted, col = "red", lwd = 2)
lines(flex1$'2018'$x, flex1$'2018'$predicted, col = "purple", lwd = 2)
#legend(8, 4, legend = c("2014", "2015","2016","2017"), fill = c("Black","blue","Grey","red"), bg = NULL, bty = "n")

#Gorgonia
plot(Size.t1~ Size.t, data = Gorg_sp[which(Gorg_sp$Year.t == '2015'),], col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(1,11),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5) #<- size against sizet+1
abline(0,1, lty = "dashed", lwd = 2) #a reference line for the null growth threshold.
lines(gorg1$'2013'$x, gorg1$'2013'$predicted, col = "Green", lwd = 2) #add the regression lines using the model coefficients
lines(gorg1$'2014'$x, gorg1$'2014'$predicted, col = "Black", lwd = 2)
lines(gorg1$'2015'$x, gorg1$'2015'$predicted, col = "blue", lwd = 2)
lines(gorg1$'2016'$x, gorg1$'2016'$predicted, col = "Grey", lwd = 2)
lines(gorg1$'2017'$x, gorg1$'2017'$predicted, col = "red", lwd = 2)
lines(gorg1$'2018'$x, gorg1$'2018'$predicted, col = "purple", lwd = 2)
legend(1,11, legend = c("2013/14","2014/15", "2015/16","2016/17","2017/18","2018/19"), fill = c("Green","Black","blue","Grey","red", "purple"), bg = NULL, bty = "n", cex = 1.2)
# one legend for the three graphs.

########## Variance in Growth #############
# How does the variability in size at t+1 change with size at t?

## A.american
# view the residuals of the chosen growth model
plot(Amer_grow1$model$Size.t,abs(resid(Amer_grow1)),xlab='size',ylab='residual') #looks pretty random - but lets see what regression analysis finds.
Amer.grow.sd = glm(abs(resid(Amer_grow1))~Amer_grow1$model$Size.t, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
Amer.grow.sd2 = lm(abs(resid(Amer_grow1))~Amer_grow1$model$Size.t) 

# which model is the best?
AIC(Amer.grow.sd,Amer.grow.sd2) # according to the AICs the gamma model fits better
BIC(Amer.grow.sd,Amer.grow.sd2) # that is agreed by the BIC.

# how does the model fit the data visually?
lines(spline(Amer_grow1$model$Size.t, exp(predict(Amer.grow.sd, list(Amer_grow1$model$Size.t))))) 

## E.flexuosa
# view the residuals of the chosen growth model
plot(flex_grow1$model$Size.t,abs(resid(flex_grow1)),xlab='size',ylab='residual') #very little variance - but lets see what regression analysis finds.
flex.grow.sd = glm(abs(resid(flex_grow1))~flex_grow1$model$Size.t, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
flex.grow.sd2 = lm(abs(resid(flex_grow1))~flex_grow1$model$Size.t) 

# which model is the best?
AIC(flex.grow.sd,flex.grow.sd2) #according to the AICs the gamma model fits better
BIC(flex.grow.sd,flex.grow.sd2) # that is agreed by the BIC.

# how does the model fit the data visually?
lines(spline(flex_grow1$model$Size.t, exp(predict(flex.grow.sd, list(flex_grow1$model$Size.t))))) 

## Gorgonia
# view the residuals of the chosen growth model
plot(gorg_grow1$model$Size.t,abs(resid(gorg_grow1)),xlab='size',ylab='residual') #very little variance - but lets see what regression analysis finds.
gorg.grow.sd = glm(abs(resid(gorg_grow1))~gorg_grow1$model$Size.t, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
gorg.grow.sd2 = lm(abs(resid(gorg_grow1))~gorg_grow1$model$Size.t) 

# which model is the best?
AIC(gorg.grow.sd,gorg.grow.sd2) #according to the AICs the gamma model fits better
BIC(gorg.grow.sd,gorg.grow.sd2) # that is agreed by the BIC.

# how does the model fit the data visually?
lines(spline(gorg_grow1$model$Size.t, exp(predict(gorg.grow.sd, list(gorg_grow1$model$Size.t)))))  

## Now plot the variance for use in publication
## Antilogorgia
plot(abs(resid(Amer_grow1)) ~ Amer_grow1$model$Size.t, col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(0,1.6),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
# add the line of best fit
abline(exp(coef(Amer.grow.sd)[1]),coef(Amer.grow.sd)[2], lwd = 2)

## Euniciea
plot(abs(resid(flex_grow1)) ~ flex_grow1$model$Size.t, col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(0,5),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
# add the line of best fit
abline(exp(coef(flex.grow.sd)[1]),coef(flex.grow.sd)[2], lwd = 2)

## Gorgonia
plot(abs(resid(gorg_grow1)) ~ gorg_grow1$model$Size.t, col = "Black", #begin by plotting data
     pch=16, cex=0.7,
     xlab="", 
     ylab="",
     ylim = c(0,5),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
# add the line of best fit
abline(exp(coef(gorg.grow.sd)[1]),coef(gorg.grow.sd)[2], lwd = 2)

########## Survival #############
# How does survival probability change with size at t over the different years?

## A. american
# plot survival to just check it
plot(Surv ~ Size.t, col = Year.t, data = A_american)
# year is always a fixed factor
ant_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = A_american)
# does the random effect of site have an effect?
ant_surv2 <- glmer(Surv ~ Size.t * Year.t + (Size.t|Site), family = "binomial", data = A_american)
# this is returning a singular fit so the model is oversaturated - ie too complex given the amount of data for each site.
# which model works best?
AIC(ant_surv1, ant_surv2) # the model excluding random effects fits best.
BIC(ant_surv1, ant_surv2) # the BIC agrees 

# Now view the survival trends on the survival plot to check the model fits
# model 1 (non-random effects)
ant.s1 <- ggpredict(ant_surv1, terms = c("Size.t[0:11]", "Year.t"))
ant.s1 <- split(ant.s1, ant.s1$group)  

lines(ant.s1$'2013'$x, ant.s1$'2013'$predicted, lty = "dashed")
lines(ant.s1$'2014'$x, ant.s1$'2014'$predicted, lty = "dashed")
lines(ant.s1$'2015'$x, ant.s1$'2015'$predicted, lty = "dashed")
lines(ant.s1$'2016'$x, ant.s1$'2016'$predicted, lty = "dashed")
lines(ant.s1$'2017'$x, ant.s1$'2017'$predicted, lty = "dashed", col = "red", lwd = 2) #Hurricane disturbance
lines(ant.s1$'2018'$x, ant.s1$'2018'$predicted, lty = "dashed") # slight recovery in the following year?
# Because of the singularity error I will follow the BIC scores.
# Survival was worst during the Hurricane year and slightly improved during the subsequent 'recovery' year.

#check the regression outputs for the model
summary(ant_surv1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(ant_surv1, test = "Chisq")
emmeans(ant_surv1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

## E. flexuosa
# plot survival to just check it
plot(Surv ~ Size.t, col = Year.t, data = E_flexuosa)
# year is always a fixed factor
flex_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = E_flexuosa)
# does the random effect of site have an effect?
flex_surv2 <- glmer(Surv ~ Size.t * Year.t + (Size.t|Site), family = "binomial", data = E_flexuosa)
# this is returning convergence warnings - may therefore not be reliable.
# which model works best?
AIC(flex_surv1, flex_surv2) # the model including random effects fits best.
BIC(flex_surv1, flex_surv2) # however the BIC disagrees suggesting there maybe so over fitting occurring (particularly since the random effects model includes a singularity error).
# So I will follow the BIC scores

# Now view the survival trends on the survival plot to check the model fit
# model 1 (non-random effects)
flex.s1 <- ggpredict(flex_surv1, terms = c("Size.t[0:11]", "Year.t"))
flex.s1 <- split(flex.s1, flex.s1$group)  

lines(flex.s1$'2014'$x, flex.s1$'2014'$predicted, lty = "dashed")
lines(flex.s1$'2015'$x, flex.s1$'2015'$predicted, lty = "dashed")
lines(flex.s1$'2016'$x, flex.s1$'2016'$predicted, lty = "dashed")
lines(flex.s1$'2017'$x, flex.s1$'2017'$predicted, lty = "dashed", col = "red", lwd = 2) # Survival of larger colonies is worst during the Hurricane year (and is low for smaller colonies - though not the worst)
lines(flex.s1$'2018'$x, flex.s1$'2018'$predicted, lty = "dashed") # Survival of small colonies recovers in recovery year but no apparent change in survival of larger colonies.

#check the regression outputs for the model
summary(flex_surv1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(flex_surv1, test = "Chisq")
emmeans(flex_surv1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

## Gorgonia
# plot survival to just check it
plot(Surv ~ Size.t, col = Year.t, data = Gorg_sp)
# year is always a fixed factor
gorg_surv1 <- glm(Surv ~ Size.t * Year.t, family = "binomial", data = Gorg_sp)
# does the random effect of site have an effect?
gorg_surv2 <- glmer(Surv ~ Size.t * Year.t + (Size.t|Site), family = "binomial", data = Gorg_sp)
# this is returning a singular fit so the model is oversaturated - ie too complex given the amount of data for each site.
# which model works best?
AIC(gorg_surv1, gorg_surv2) # the model excluding random effects fits best.
BIC(gorg_surv1, gorg_surv2) # the BIC agrees 

# Now view the survival trends on the survival plot to check the model fits
# model 1 (non-random effects)
gorg.s1 <- ggpredict(gorg_surv1, terms = c("Size.t[0:11]", "Year.t"))
gorg.s1 <- split(gorg.s1, gorg.s1$group)  

lines(gorg.s1$'2013'$x, gorg.s1$'2013'$predicted, lty = "dashed")
lines(gorg.s1$'2014'$x, gorg.s1$'2014'$predicted, lty = "dashed")
lines(gorg.s1$'2015'$x, gorg.s1$'2015'$predicted, lty = "dashed")
lines(gorg.s1$'2016'$x, gorg.s1$'2016'$predicted, lty = "dashed")
lines(gorg.s1$'2017'$x, gorg.s1$'2017'$predicted, lty = "dashed", col = "red", lwd = 2) #Hurricane disturbance
lines(gorg.s1$'2018'$x, gorg.s1$'2018'$predicted, lty = "dashed") #slight recovery in the following year.
# Survival is more variable but worst during the Hurricane year and slightly improved during the subsequent 'recovery' year.
# Equally the poor survival of large colonies compared to smaller colonies during the hurricane year implicates the broad surface area of gorgonia fans in highly volitaile waters.

#check the regression outputs for the model
summary(gorg_surv1)
# this line of code confirms whether both fixed variables have an effect on Size at t+1
anova(gorg_surv1, test = "Chisq")
emmeans(gorg_surv1, pairwise ~ Year.t)$contrasts #since year has a significant effect this line of code identifies where this lies

##########   
## A. american
plot(Surv~Size.t, type = "n", data = A_american, #now we can plot the data seperatly - start by plotting a blank plot
     # the x axis needs to match the one used when predicting values.
     xlab="", 
     ylab="",
     ylim = c(0,1),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
lines(spline(ant.s1$'2013'$predicted ~ ant.s1$'2013'$x), col = "green", lwd = 2)
lines(spline(ant.s1$'2014'$predicted ~ ant.s1$'2014'$x), col = "black", lwd = 2)
lines(spline(ant.s1$'2015'$predicted ~ ant.s1$'2015'$x), col = "blue", lwd = 2)
lines(spline(ant.s1$'2016'$predicted ~ ant.s1$'2016'$x), col = "grey", lwd = 2)
lines(spline(ant.s1$'2017'$predicted ~ ant.s1$'2017'$x), col = "red", lwd = 2)
lines(spline(ant.s1$'2018'$predicted ~ ant.s1$'2018'$x), col = "purple", lwd = 2)

## Eunicea flexuosa
plot(Surv~Size.t, type = "n", data = E_flexuosa, #the data can be plotted straight in..
     xlab="", 
     ylab="",
     ylim = c(0,1),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
lines(spline(flex.s1$'2014'$predicted ~ flex.s1$'2014'$x), col = "black", lwd = 2)
lines(spline(flex.s1$'2015'$predicted ~ flex.s1$'2015'$x), col = "blue", lwd = 2)
lines(spline(flex.s1$'2016'$predicted ~ flex.s1$'2016'$x), col = "grey", lwd = 2)
lines(spline(flex.s1$'2017'$predicted ~ flex.s1$'2017'$x), col = "red", lwd = 2)
lines(spline(flex.s1$'2018'$predicted ~ flex.s1$'2018'$x), col = "purple", lwd = 2)

## Gorgonia
plot(Surv~Size.t, type = "n", data = Gorg_sp, #now we can plot the data separately - start by plotting a blank plot
     # the x axis needs to match the one used when predicting values.
     xlab="", 
     ylab="",
     ylim = c(0,1),
     xlim = c(1,11), xaxs = "i", yaxs = "i", cex.axis = 1.5)
lines(spline(gorg.s1$'2013'$predicted ~ gorg.s1$'2013'$x), col = "green", lwd = 2)
lines(spline(gorg.s1$'2014'$predicted ~ gorg.s1$'2014'$x), col = "black", lwd = 2)
lines(spline(gorg.s1$'2015'$predicted ~ gorg.s1$'2015'$x), col = "blue", lwd = 2)
lines(spline(gorg.s1$'2016'$predicted ~ gorg.s1$'2016'$x), col = "grey", lwd = 2)
lines(spline(gorg.s1$'2017'$predicted ~ gorg.s1$'2017'$x), col = "red", lwd = 2)
lines(spline(gorg.s1$'2018'$predicted ~ gorg.s1$'2018'$x), col = "purple", lwd = 2)
legend(5, 0.5, legend = c("2013/14","2014/15", "2015/16","2016/17","2017/18","2018/19"), fill = c("green","Black","blue","Grey","red","purple"), bg = NULL, bty = "n", cex = 1.2)

##############################
#Step 3: Determine species specific size distributions (recruits & full populations)
##############################
# this part of the script involves data collected during surveys by Howard Lasker and Angela Martinez Quintana.

#####-------------------- Determine most recent recorded size distribution - all sizes (2019)
# Load quadrat survey data collected in 2019.
Adults19 <- read.csv('Adult Census_2019.csv')
Recruits19 <- read.csv('Recruit Census_2019.csv')
# these data files contains a lot of detail that isn't needed here.
# firstly remove non-spring/summer survey data (this is out of sequence).
Adults19 <- Adults19[which(substr(Adults19$Date, 4,6) %in% c("Jul","Aug")),]
Recruits19 <- Recruits19[which(substr(Recruits19$Date, 4,6) %in% c("Jul","Aug")),]

# Next, the data needs formatting so that it only contains the species I am interested in and that height is on the same scale as the tagged data set.
Adults19 <- subset(Adults19, Height.of.living.tissue..cm. != "nd" & #remove the 'none' height entries.
                   Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") & #these are the only genera that we are interested in from this datafile
                   Species.code %in% c("aa","gv", 'aa2', "efl"))
Recruits19 <- subset(Recruits19, Height..cm. != "nd" & #remove the 'none' height entries.
                     Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") & #these are the only genera that we are interested in from this datafile
                     Species.code %in% c("aa","gv", "efl", "efl "))

# ensure size variable names match across files
names(Adults19)[9] <- names(Recruits19)[10]
# And stitch together the files.
# NOTE: Recruit surveys were carried out using 0.25m2 quadrats and so need to be scaled by a factor of 4 to reflect the same spatial scale used in the adult colony surveys
pop2019 <- dplyr::bind_rows(Adults19, Recruits19, Recruits19, Recruits19, Recruits19)
# Convert height data into numeric form
pop2019$Height..cm. <- as.numeric(paste(pop2019$Height..cm.))
# and square root transform to match sizes range used across the survival and growth data
pop2019$Height..cm. <- sqrt(pop2019$Height..cm.)

# Extract the size distributions (pooled across sites) recorded for Eunicea, Antillogorgia, and Gorgonia in 2019
gorg.2019 <- pop2019[which(pop2019$Genus == "Gorgonia"),]
ant.2019 <-  pop2019[which(pop2019$Genus == "Antillogorgia"),]
flex.2019 <- pop2019[which(pop2019$Genus == "Eunicea"),]
initial.pop <- list(ant.2019$Height..cm., flex.2019$Height..cm., gorg.2019$Height..cm.) #store them

#####-------------------- Determine species-specific recruit size distributions - pooled across years
# Load recruit quadrat survey data collected between 2014 and 2018 (2019 already loaded).
Recruits14_18 <- read.csv('Recruits Censuses_2014_2018.csv')
# Again remove non-spring/summer survey data.
Recruits14_18[which(nchar(Recruits14_18$Date) < 9), "Date"] <- paste0('0', Recruits14_18[which(nchar(Recruits14_18$Date) < 9), "Date"]) # to ensure consistency in date format
Recruits14_18 <- Recruits14_18[which(substr(Recruits14_18$Date, 4,6) %in% c("Jul","Aug")),]
# Again, reformat the nessecary data.
Recruits14_18 <- subset(Recruits14_18, Height..cm. != "nd" & #remove the 'none' height entries.
                        Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") & #these are the only genera that we are interested in from this datafile
                        Species.code %in% c("aa","gv", "efl"))
# And stitch together the files.
Recruits14_18$Meter.on.transect <- as.numeric(paste(Recruits14_18$Meter.on.transect))
Recruits_Full <- dplyr::bind_rows(Recruits14_18, Recruits19)
# and convert height to numeric form 
Recruits_Full$Height..cm. <- as.numeric(paste(Recruits_Full$Height..cm.))
# and square root transform to match sizes range used across the survival and growth data
Recruits_Full$Height..cm. <- sqrt(Recruits_Full$Height..cm.)

# Split recruit data into the three different genera.
ant.recruits <- subset(Recruits_Full, Genus == "Antillogorgia")
gorg.recruits <- subset(Recruits_Full, Genus == "Gorgonia")
eunicea.recruits <- subset(Recruits_Full, Genus == "Eunicea")

# Now determine the recruit size distributions for these three octocoral species.
## Gorgonia
hist(gorg.recruits$Height..cm., breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,40),
     xaxs = "i", yaxs = "i") #could be a normal distribution - but it is continuous and none negative.
# the fitdistrplus package contains functions for determining these types of distributions 
# the package allows for the plotting of a Cullen and Frey graph to determine the distribution of the data.
descdist(gorg.recruits$Height..cm., boot = 1000) #according to this the data is very centred and close to a normal and uniform distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
g.norm <- fitdist(gorg.recruits$Height..cm., distr = "norm", method = "mle") #normal distribution
g.uni <- fitdist(gorg.recruits$Height..cm., distr = "unif", method = "mle") #uniform distribution
g.gam <- fitdist(gorg.recruits$Height..cm., distr = "gamma", method = "mle") #gamma distribution
g.bull <- fitdist(gorg.recruits$Height..cm., distr = "weibull", method = "mle") # weibull
g.log <- fitdist(gorg.recruits$Height..cm., distr = "lnorm", method = "mle") #log normal
# plot them together to see which fits best
denscomp(list(g.norm, g.uni, g.gam, g.bull, g.log)) # visually the normal and weibull distributions are the best!
summary(g.norm)$bic; summary(g.uni)$bic; summary(g.gam)$bic; summary(g.bull)$bic; summary(g.log)$bic # weibull distribution is the best.

## A. american
hist(ant.recruits$Height..cm., breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,20),
     xaxs = "i", yaxs = "i") # this one is skewed towards larger recruit sizes.
# check the distribution.
descdist(ant.recruits$Height..cm., boot = 1000) #according to this the data is not very centred and close to a normal and gamma, and lognormal distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
a.norm <- fitdist(ant.recruits$Height..cm., distr = "norm", method = "mle") #normal distribution
# I will also test the fits of more skewed distributions
a.gam <- fitdist(ant.recruits$Height..cm., distr = "gamma", method = "mle") #gamma distribution
a.bull <- fitdist(ant.recruits$Height..cm., distr = "weibull", method = "mle") # weibull
a.log <- fitdist(ant.recruits$Height..cm., distr = "lnorm", method = "mle") #log normal
# plot them all together to see which fits best
denscomp(list(a.norm,a.gam, a.bull, a.log)) #visually the normal and weibull distributions are the best!
#compare model AICs and BICs
summary(a.norm)$aic; summary(a.gam)$aic; summary(a.bull)$aic; summary(a.log)$aic # Weibull is the best mathematical fit
summary(a.norm)$bic; summary(a.gam)$bic; summary(a.bull)$bic; summary(a.log)$bic

## E. flexuosa
hist(eunicea.recruits$Height..cm., breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,30),
     xaxs = "i", yaxs = "i") # doesn't look like  a normal distribution
# check the distribution.
descdist(eunicea.recruits$Height..cm., boot = 1000) #according to this the data is very centred and close to a normal and uniform distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
e.norm <- fitdist(eunicea.recruits$Height..cm., distr = "norm", method = "mle") #normal distribution
# I will also test the fits of more skewed distributions
e.gam <- fitdist(eunicea.recruits$Height..cm., distr = "gamma", method = "mle") #gamma distribution
e.bull <- fitdist(eunicea.recruits$Height..cm., distr = "weibull", method = "mle") # weibull
e.log <- fitdist(eunicea.recruits$Height..cm., distr = "lnorm", method = "mle") #log normal
# plot them all together to see which fits best
denscomp(list(e.norm,e.gam, e.bull, e.log)) #visually the normal and weibull distributions are the best!
#compare model AICs and BICs
summary(e.norm)$aic; summary(e.gam)$aic; summary(e.bull)$aic; summary(e.log)$aic # weibull fits the best mathematically
summary(e.norm)$bic; summary(e.gam)$bic; summary(e.bull)$bic; summary(e.log)$bic

#A.american
denscomp(a.bull, #what density is being plotted with the data
         xlim = c(0,2.5),
         ylim = c(0,1.2),
         main = NULL,
         xlab = NULL,
         ylab = NULL,
         addlegend = FALSE,
         xaxs = "i", yaxs = "i",
         cex.axis = 1.5)
#E.flexuosa
denscomp(e.bull, #what density is being plotted with the data
         xlim = c(0,2.5),
         ylim = c(0,1.2),
         main = NULL,
         xlab = NULL,
         ylab = NULL,
         addlegend = FALSE,
         xaxs = "i", yaxs = "i",
         cex.axis = 1.5)
# Gorgonia
denscomp(g.bull, #what density is being plotted with the data
         xlim = c(0,2.5),
         ylim = c(0,1.2),
         main = NULL,
         #xlab = NULL,
         #ylab = NULL,
         addlegend = FALSE,
         xaxs = "i", yaxs = "i",
         cex.axis = 1.5)


##############################
#Step 4: Determine the density-dependent recruitment function.
##############################
# Ensure all recruit and adult census data from 2014 to 2019 is reloaded
# Recruitment data can be downloaded from the BCO-DMO database (https://www.bco-dmo.org/dataset/851382/data and https://www.bco-dmo.org/dataset/851382/data )
Adults1 <- read.csv('Adult Censuses_2014_2018.csv')
Adults2 <- read.csv('Adult Census_2019.csv')
Recruits1 <- read.csv('Recruits Censuses_2014_2018.csv')
Recruits2 <- read.csv('Recruit Census_2019.csv')

# Again these data files contains a lot of detail that isn't needed here.
# Remove non-spring/summer survey data.
Recruits1[which(nchar(Recruits1$Date) < 9), "Date"] <- paste0('0', Recruits1[which(nchar(Recruits1$Date) < 9), "Date"]) # to ensure consistency in date format
Adults1[which(nchar(Adults1$Date) < 9), "Date"] <- paste0('0', Adults1[which(nchar(Adults1$Date) < 9), "Date"])
# remove data not collected in June, Jul, or Aug
Recruits1 <- Recruits1[which(substr(Recruits1$Date, 4,6) %in% c("Jun","Jul","Aug")),]
Recruits2 <- Recruits2[which(substr(Recruits2$Date, 4,6) %in% c("Jun","Jul","Aug")),]
Adults1 <- Adults1[which(substr(Adults1$Date, 4,6) %in% c("Jun","Jul","Aug")),]
Adults2 <- Adults2[which(substr(Adults2$Date, 4,6) %in% c("Jun","Jul","Aug")),]

# Next, reformat the nessecary data.
Recruits1 <- subset(Recruits1,
                    Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") & #these are the only genera that we are interested in from this datafile
                    Species.code %in% c("aa","gv", "efl"))
Recruits2 <- subset(Recruits2,
                    Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") &
                    Species.code %in% c("aa","gv", "efl", "efl "))
Adults1 <- subset(Adults1,
                  Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") &
                  Species.code %in% c("aa","gv", "efl", "efl "))
Adults2 <- subset(Adults2,
                  Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia") & 
                  Species.code %in% c("aa", "aa2", "gv", "efl"))

# And stitch together the files.
Recruits1$Meter.on.transect <- as.numeric(paste(Recruits1$Meter.on.transect))
Adults2$Transect.Position <- as.character(Adults2$Transect.Position)
Recruits <- dplyr::bind_rows(Recruits1, Recruits2)
Adults <- dplyr::bind_rows(Adults1, Adults2)
# convert size data to numeric
Recruits$Height..cm. <- as.numeric(paste(Recruits$Height..cm.))
Adults$Height.of.living.tissue..cm. <- as.numeric(paste(Adults$Height.of.living.tissue..cm.))

# Now I need to use this data to populate a new data frame extracting annual counts of individuals per quadrat (across sites).
# first extract grouped summaries
Adult_counts <- Adults %>% group_by(Census.Year, Site, Transect.Position, Meter.on.transect, Genus) %>% summarise(n = n())
Recruit_counts <- Recruits %>% group_by(Census.Year, Site, Transect.Position, Meter.on.transect, Genus) %>% summarise(n = n())
# Next determine all unique transect position, site and census year combinations
Adult_transects <- unique(Adult_counts[c("Census.Year","Site","Transect.Position")])
Recruit_transects <- unique(Recruit_counts[c("Census.Year","Site","Transect.Position")])

# use these to define blank data frames to determine colonies recorded per quadrat (expanding the dataframe to include quadrats in which no individuals were recorded)
Adult_df <- data.frame(Site = rep(Adult_transects$Site, each = 10),
                       Year = rep(Adult_transects$Census.Year, each = 10),
                       Transect = rep(Adult_transects$Transect.Position, each = 10),
                       Quadrat = rep(1:10, length.out = dim(Adult_transects)[1]*10),
                       Count = rep(0, length.out = dim(Adult_transects)[1]*10))
Recruit_df <- data.frame(Site = rep(Recruit_transects$Site, each = 8),
                         Year = rep(Recruit_transects$Census.Year, each = 8),
                         Transect = rep(Recruit_transects$Transect.Position, each = 8),
                         Quadrat = rep(1:8, length.out = dim(Recruit_transects)[1]*8),
                         Count = rep(0, length.out = dim(Recruit_transects)[1]*8))                       
# Bring dfs together
Adult_list <- list(Adult_df, Adult_df, Adult_df); names(Adult_list) <- c('AA', 'GV', 'EF') # one dataframe for each species.
Recruit_list <- list(Recruit_df, Recruit_df, Recruit_df); names(Recruit_list) <- c('AA', 'GV', 'EF')

# And reformat the transect index tables
Adult_transects <- as.data.frame(Adult_transects)
Recruit_transects <- as.data.frame(Recruit_transects)

# Populate dataframe
# loop through in sequence extract species counts (if no count present for a given combination the count is entered as zero, i.e. no observations)
# Adult colonies
for(ii in 1:dim(Adult_transects)[1]){
  # progress read out
  print(ii)
  
  # extract counts associated with selected site, year and transect combination
  Adult_use <- Adult_counts[which(Adult_counts$Census.Year == Adult_transects[ii, 1] &
                                  Adult_counts$Site == Adult_transects[ii, 2] &
                                  Adult_counts$Transect.Position == Adult_transects[ii, 3]),]
  # split by species
  AA_adults <- pull(Adult_use[which(Adult_use$Genus == 'Antillogorgia'),], n)
  GV_adults <- pull(Adult_use[which(Adult_use$Genus == 'Gorgonia'),], n)
  EF_adults <- pull(Adult_use[which(Adult_use$Genus == 'Eunicea'),], n)
  
  # Insert counts
  if(length(AA_adults)>0){ for(xx in 1:length(AA_adults)){ Adult_list$AA[which(Adult_list$AA$Site == Adult_transects[ii,2] &
                                                                               Adult_list$AA$Year == Adult_transects[ii,1] &
                                                                               Adult_list$AA$Transect == Adult_transects[ii,3]), 'Count'][xx] <- AA_adults[xx] }}
  if(length(GV_adults)>0){ for(xx in 1:length(GV_adults)){ Adult_list$GV[which(Adult_list$GV$Site == Adult_transects[ii,2] &
                                                                               Adult_list$GV$Year == Adult_transects[ii,1] &
                                                                               Adult_list$GV$Transect == Adult_transects[ii,3]), 'Count'][xx] <- GV_adults[xx] }}
  if(length(EF_adults)>0){ for(xx in 1:length(EF_adults)){ Adult_list$EF[which(Adult_list$EF$Site == Adult_transects[ii,2] &
                                                                               Adult_list$EF$Year == Adult_transects[ii,1] &
                                                                               Adult_list$EF$Transect == Adult_transects[ii,3]), 'Count'][xx] <- EF_adults[xx] }}
  # clear for next run through
  rm(AA_adults, GV_adults, EF_adults, Adult_use)
}

# Recruit colonies
for(ii in 1:dim(Recruit_transects)[1]){
  # progress read out
  print(ii) 
  
  # extract counts associated with selected site, year and transect combination
  Recruit_use <- Recruit_counts[which(Recruit_counts$Census.Year == Recruit_transects[ii,1] &
                                      Recruit_counts$Site == Recruit_transects[ii,2] &
                                      Recruit_counts$Transect.Position == Recruit_transects[ii,3]),]
  # split by species
  AA_recruits <- pull(Recruit_use[which(Recruit_use$Genus == 'Antillogorgia'),], n)
  GV_recruits <- pull(Recruit_use[which(Recruit_use$Genus == 'Gorgonia'),], n)
  EF_recruits <- pull(Recruit_use[which(Recruit_use$Genus == 'Eunicea'),], n)
  
  # Insert counts
  if(length(AA_recruits)>0){ for(xx in 1:length(AA_recruits)){ Recruit_list$AA[which(Recruit_list$AA$Site == Recruit_transects[ii,2] &
                                                                                Recruit_list$AA$Year == Recruit_transects[ii,1] &
                                                                                Recruit_list$AA$Transect == Recruit_transects[ii,3]), 'Count'][xx] <- AA_recruits[xx]*4 }} # the *4 here is to ensure that the spatial resolution of recruits and adults match
  if(length(GV_recruits)>0){ for(xx in 1:length(GV_recruits)){ Recruit_list$GV[which(Recruit_list$GV$Site == Recruit_transects[ii,2] &
                                                                                Recruit_list$GV$Year == Recruit_transects[ii,1] &
                                                                                Recruit_list$GV$Transect == Recruit_transects[ii,3]), 'Count'][xx] <- GV_recruits[xx]*4 }}
  if(length(EF_recruits)>0){ for(xx in 1:length(EF_recruits)){ Recruit_list$EF[which(Recruit_list$EF$Site == Recruit_transects[ii,2] &
                                                                                Recruit_list$EF$Year == Recruit_transects[ii,1] &
                                                                                Recruit_list$EF$Transect == Recruit_transects[ii,3]), 'Count'][xx] <- EF_recruits[xx]*4 }}
  # clear for next run through
  rm(AA_recruits, GV_recruits, EF_recruits, Recruit_use)
}

# Okay so now the individual counts for adult and recruit colonies have been extracted. 
# A site mean recruit density and site mean adult density need to be estimated for each species in each year.
# How many site year combinations are being dealt with.
Adult_combos <- unique(Adult_transects[c('Census.Year', 'Site')])
Recruit_combos <- unique(Recruit_transects[c('Census.Year', 'Site')])
# define a blank storage dataframe
DD_df <- data.frame(Year = Adult_combos$Census.Year,
                    Site = Adult_combos$Site,
                    AA_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    AA_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    aa_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    aa_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    GV_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    GV_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    gv_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    gv_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    EF_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    EF_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    ef_mean = rep(NA, length.out = dim(Adult_combos)[1]),
                    ef_sd = rep(NA, length.out = dim(Adult_combos)[1]),
                    TOTAL_mean = rep(NA, length.out = dim(Adult_combos)[1]))

# Now work through each combination to determine a mean recruit and adult density
for(ii in 1:dim(Adult_combos)[1]){
  # progress read out
  print(ii)
  
  # extract site*year observations for each species
  Adult_AA <- Adult_list$AA[which(Adult_list$AA$Year == DD_df[ii,1] &
                                  Adult_list$AA$Site == DD_df[ii,2]), 'Count']*270 # the multiplication here ensures that the relationship between number of adults and recruits is being quantified at the same scale at which initial population size has been estimated
  Recruit_AA <- Recruit_list$AA[which(Recruit_list$AA$Year == DD_df[ii,1] &
                                      Recruit_list$AA$Site == DD_df[ii,2]), 'Count']*216
  Adult_GV <- Adult_list$GV[which(Adult_list$GV$Year == DD_df[ii,1] &
                                  Adult_list$GV$Site == DD_df[ii,2]), 'Count']*270
  Recruit_GV <- Recruit_list$GV[which(Recruit_list$GV$Year == DD_df[ii,1] &
                                      Recruit_list$GV$Site == DD_df[ii,2]), 'Count']*216
  Adult_EF <- Adult_list$EF[which(Adult_list$EF$Year == DD_df[ii,1] &
                                  Adult_list$EF$Site == DD_df[ii,2]), 'Count']*270
  Recruit_EF <- Recruit_list$EF[which(Recruit_list$EF$Year == DD_df[ii,1] &
                                      Recruit_list$EF$Site == DD_df[ii,2]), 'Count']*216
  
  # insert summary stats into the blank datafile
  # Adult counts
  if(length(Adult_AA)==0) { DD_df[ii, 'AA_mean'] <- 0 } else { DD_df[ii, 'AA_mean'] <- mean(Adult_AA, na.rm = TRUE) }
  if(length(Adult_AA)==0) { DD_df[ii, 'AA_sd'] <- 0 } else { DD_df[ii, 'AA_sd'] <- sd(Adult_AA, na.rm = TRUE) }
  # If no adults observed then this represents true zero.
  if(length(Adult_GV)==0) { DD_df[ii, 'GV_mean'] <- 0 } else { DD_df[ii, 'GV_mean'] <- mean(Adult_GV, na.rm = TRUE) }
  if(length(Adult_GV)==0) { DD_df[ii, 'GV_sd'] <- 0 } else { DD_df[ii, 'GV_sd'] <- sd(Adult_GV, na.rm = TRUE) }
  if(length(Adult_EF)==0) { DD_df[ii, 'EF_mean'] <- 0 } else { DD_df[ii, 'EF_mean'] <- mean(Adult_EF, na.rm = TRUE) }
  if(length(Adult_EF)==0) { DD_df[ii, 'EF_sd'] <- 0 } else { DD_df[ii, 'EF_sd'] <- sd(Adult_EF, na.rm = TRUE) }
  
  # Recruit counts
  if(length(Recruit_AA)==0) { DD_df[ii, 'aa_mean'] <- 0 } else { DD_df[ii, 'aa_mean'] <- mean(Recruit_AA, na.rm = TRUE) }
  if(length(Recruit_AA)==0) { DD_df[ii, 'aa_sd'] <- 0 } else { DD_df[ii, 'aa_sd'] <- sd(Recruit_AA, na.rm = TRUE) }
  if(length(Recruit_GV)==0) { DD_df[ii, 'gv_mean'] <- 0 } else { DD_df[ii, 'gv_mean'] <- mean(Recruit_GV, na.rm = TRUE) }
  if(length(Recruit_GV)==0) { DD_df[ii, 'gv_sd'] <- 0 } else { DD_df[ii, 'gv_sd'] <- sd(Recruit_GV, na.rm = TRUE) }
  if(length(Recruit_EF)==0) { DD_df[ii, 'ef_mean'] <- 0 } else { DD_df[ii, 'ef_mean'] <- mean(Recruit_EF, na.rm = TRUE) }
  if(length(Recruit_EF)==0) { DD_df[ii, 'ef_sd'] <- 0 } else { DD_df[ii, 'ef_sd'] <- sd(Recruit_EF, na.rm = TRUE) }
  # However, no recruit counts does not necessarily reflect a true zero - recruit surveys where carried out during less years than adult surveys.
  # This therefore needs to be corrected in the data
  survey.test <- filter(Recruit_combos, Census.Year == DD_df[ii, 'Year'], Site == DD_df[ii, 'Site']) # determine whether recruit surveys where carried out at the specified time and location
  # remove entries if not.
  if(dim(survey.test)[1] == 0) { DD_df[ii, 'aa_mean'] <- NA 
                                 DD_df[ii, 'aa_sd'] <- NA
                                 DD_df[ii, 'gv_mean'] <- NA 
                                 DD_df[ii, 'gv_sd'] <- NA
                                 DD_df[ii, 'ef_mean'] <- NA 
                                 DD_df[ii, 'ef_sd'] <- NA }
  
  # Clear for next run through
  rm(Adult_AA, Adult_GV, Adult_EF, Recruit_AA, Recruit_GV, Recruit_EF)
}

# calculate total adult counts
DD_df[,'TOTAL_mean'] <- rowSums(DD_df[,c('AA_mean','GV_mean','EF_mean')], na.rm= TRUE)
# Remove rows for which no recruit surveys were carried out
DD_df <- DD_df[which(complete.cases(DD_df)),]

# Now this dataframe can be used to explore the relationship between adult and recruit densities.
# Here we have site means severing as replicates.

## Before moving on with the analysis of the density-dependent recruitment function there are a few concepts to iron out.
# 1) In their 2015 paper, Privitera-Johnson et al. used site means as the dependent variables with the different sites representing the replicates, and estimated the relationship between juvenile and adult densities using a type 2 regression approach.
# This was to reflect the error present in both counts of adult and juvenile colonies. However type-2 regression is not considered as effective when using x variables to predict y variables. Thus it is frequently argued that OLS regression is the prefer approach in prediction exercises.
# 2) Colony densities are a form of count data. However having been converted to mean counts the data is now on a continuous scale and so is no longer suitable for a poisson regression.
# However the data still consists of non-negative values. The data is therefore now suitable for gamma analysis.
# 3) But the juvenile count data contains zero entries (not accepted in a normal gamma equation). 
# The approach used here therefore will be a zero-adjusted gamma model.

# Now determine whether patterns in juvenile densities for each species correspond more with con-specific adult densities or total adult densities
# This test will be carried out using distance correlation which evaluates for both linear and non-linear associations between variables

#### Gorgonia
dcor(DD_df[,9], DD_df[,7])
dcor(DD_df[,9], DD_df[,15])
plot(y = DD_df[,9], x = DD_df[,7], xlab = "Adults", ylab = "Recruits") # from this plot there is a clear outlier that is distorting the model and requires omitting
DD_df[14, c("GV_mean","gv_mean")] <- NA
# re-run the test
dcor(na.omit(DD_df[,9]), na.omit(DD_df[,7]))
dcor(na.omit(DD_df[,9]), log(DD_df[-14,15]))
plot(y = DD_df[,9], x = DD_df[,7], xlab = "Adults", ylab = "Recruits")
#### Antilogorgia
dcor(DD_df[,5], DD_df[,3]) 
dcor(DD_df[,5], DD_df[,15])
plot(y = DD_df[,5], x = DD_df[,15], xlab = "Adults", ylab = "Recruits")
#### Eunicea
dcor(DD_df[,13], DD_df[,11])
dcor(DD_df[,13], DD_df[,15])
plot(y = DD_df[,13], x = DD_df[,11], xlab = "Adults", ylab = "Recruits")
# Juvenile densities correlate best with total adult density in Antillogorgia but not in Gorgonia and Eunicea

# From previous efforts it appears the model breaks when eunicea is allowed exponetial recruitment, it needs some for of biological control built in.
# Therefore the model will assume inter-specific density dependance across species.

# Perform zero-adjusted gamma regression analysis.
## Gorgonia    
Gorg.juv <- na.omit(DD_df[,9]) # extract the required data
Gorg.adults <- (DD_df[-14,15])
gorg.dd.mod <- gamlss(Gorg.juv ~ Gorg.adults, family = ZAGA()) #gamlss performs the same function as gam but allows for the use of ZAGA
# this is a zero adjusted gamma which allows for data to have zero values.
gorg.dd.mod2 <- gamlss(Gorg.juv ~ poly(Gorg.adults, degree = 2, raw = TRUE), family = ZAGA()) # test model fit using a polynomial 
# which model has the better fit
LR.test(gorg.dd.mod, gorg.dd.mod2) # not a significant difference but the polynomial model carries the lower deviance. 
 # store the predicted recruitment values
gorg.recruit.store <- data.frame(Adults = Gorg.adults,
                                 Recruits = predict(gorg.dd.mod2, type = "response")) 
gorg.recruit.store <- gorg.recruit.store[order(gorg.recruit.store$Adults),] # re-order the adult numbers in ascending order ready for plotting

## Antillogorgia    
Ant.juv <- DD_df[,5]  # this adjusts the relationship to remove negative values (this adjustment will be removed to determine actual recruitment levels)
Ant.adults <- DD_df[,15] 
Ant.dd.mod <- gamlss(Ant.juv ~ Ant.adults, family = ZAGA())
Ant.dd.mod2 <- gamlss(Ant.juv ~ poly(Ant.adults, degree = 2, raw = TRUE), family = ZAGA()) 
# which model has the better fit
LR.test(Ant.dd.mod, Ant.dd.mod2) # not a significant difference but the polynomial model carries the lower deviance. 
# store the predicted recruitment values
Ant.recruit.store <- data.frame(Adults = Ant.adults,
                                Recruits = predict(Ant.dd.mod2, type = "response")) 
Ant.recruit.store <- Ant.recruit.store[order(Ant.recruit.store$Adults),] #re order the adult numbers in ascending order ready for plotting

## Euniciea    
flex.juv <- DD_df[,13]  
flex.adults <- DD_df[,15] 
flex.dd.mod <- gamlss(flex.juv ~ flex.adults, family = ZAGA()) 
flex.dd.mod2 <- gamlss(flex.juv ~ poly(flex.adults, degree = 2, raw = TRUE), family = ZAGA()) 
# which model has the better fit
LR.test(flex.dd.mod, flex.dd.mod2) # not a significant difference but the polynomial model carries the lower deviance. 
# store the predicted recruitment values
flex.recruit.store <- data.frame(Adults = flex.adults, # exp is for returning non-logged entries
                                 Recruits = predict(flex.dd.mod2, type = "response"))
flex.recruit.store <- flex.recruit.store[order(flex.recruit.store$Adults),] #re order the adult numbers in ascending order ready for plotting

# now I can plot the density-dependent trend
set_graph_pars("panel3")      

## Antillogorgia
plot(Ant.juv ~ Ant.adults, col = "black", #the data can be plotted straight in.
     pch=16, cex=0.8,
     xlab = "", # "Total adult density (60m^2),
     ylab = "", # "Recruit density (60m^2),
     ylim = c(0,250), 
     xlim = c(230, 900), 
     xaxs = "i", yaxs = "i", cex.axis = 1.5)
# y error
lines(spline(Ant.recruit.store$Recruits ~ Ant.recruit.store$Adults), col = "red", lwd = 1)

## Euniciea
plot(flex.juv ~ flex.adults, col = "black", #the data can be plotted straight in.
     pch=16, cex=0.8,
     xlab = "", # "log(Total adult density),
     ylab = "", # "log(Recruit density),
     ylim = c(0, 600), 
     xlim = c(230, 900), xaxs = "i", yaxs = "i", cex.axis = 1.5)
lines(spline(flex.recruit.store$Recruits ~ flex.recruit.store$Adults), col = "red", lwd = 1)

## Gorgonia
plot(Gorg.juv ~ Gorg.adults, col = "black", #the data can be plotted straight in.
     pch=16, cex=0.8,
     xlab = "", # "Con-specific adult density (60m^2),
     ylab = "", # "Recruit density (60m^2),
     ylim = c(0,600),
     xlim = c(230,900), xaxs = "i", yaxs = "i", cex.axis = 1.5)
lines(spline(gorg.recruit.store$Recruits ~ gorg.recruit.store$Adults), col = "red", lwd = 1)

# The use of ZAGA already includes a measure of variability that will be parameterised into the IPM models

##############################
#Step 5: Store model coefficients for construction of IPM.
##############################

# blank matrices to store the yearly parameters for each species separately
gorg.params <- matrix(NA, 15, 6)
flex.params <- matrix(NA, 15, 5)
ant.params <- matrix(NA, 15, 6)

# Fill each matrix.
# survival
for(i in 1:6){
  if (i == 1){
    ant.params[1,i] <- coefficients(ant_surv1)[1]
    ant.params[2,i] <- coefficients(ant_surv1)[2]
    ant.params[3,i] <- sigma(ant_surv1)
  } else {
    ant.params[1,i] <- coefficients(ant_surv1)[1] + coefficients(ant_surv1)[i+1] 
    ant.params[2,i] <- coefficients(ant_surv1)[2] + coefficients(ant_surv1)[i+6]
    ant.params[3,i] <- sigma(ant_surv1)
  }  
}

for(i in 1:6){
  if (i == 1){
    gorg.params[1,i] <- coefficients(gorg_surv1)[1]
    gorg.params[2,i] <- coefficients(gorg_surv1)[2]
    gorg.params[3,i] <- sigma(gorg_surv1)
  } else {
    gorg.params[1,i] <- coefficients(gorg_surv1)[1] + coefficients(gorg_surv1)[i+1] 
    gorg.params[2,i] <- coefficients(gorg_surv1)[2] + coefficients(gorg_surv1)[i+6]
    gorg.params[3,i] <- sigma(gorg_surv1)
  }  
}

for(i in 1:5){
  if (i == 1){
    flex.params[1,i] <- coefficients(flex_surv1)[1]
    flex.params[2,i] <- coefficients(flex_surv1)[2]
    flex.params[3,i] <- sigma(flex_surv1)
  } else {
    flex.params[1,i] <- coefficients(flex_surv1)[1] + coefficients(flex_surv1)[i+1] 
    flex.params[2,i] <- coefficients(flex_surv1)[2] + coefficients(flex_surv1)[i+5]
    flex.params[3,i] <- sigma(flex_surv1)
  }  
}

# growth
for(i in 1:6){
  if (i == 1){
    ant.params[4,i] <- coefficients(Amer_grow1)[1]
    ant.params[5,i] <- coefficients(Amer_grow1)[2]
  } else {
    ant.params[4,i] <- coefficients(Amer_grow1)[1] + coefficients(Amer_grow1)[i+1] 
    ant.params[5,i] <- coefficients(Amer_grow1)[2] + coefficients(Amer_grow1)[i+6]
  }  
}

for(i in 1:6){
  if (i == 1){
    gorg.params[4,i] <- coefficients(gorg_grow1)[1]
    gorg.params[5,i] <- coefficients(gorg_grow1)[2]
  } else {
    gorg.params[4,i] <- coefficients(gorg_grow1)[1] + coefficients(gorg_grow1)[i+1] 
    gorg.params[5,i] <- coefficients(gorg_grow1)[2] + coefficients(gorg_grow1)[i+6]
  }  
}

for(i in 1:5){
  if (i == 1){
    flex.params[4,i] <- coefficients(flex_grow1)[1]
    flex.params[5,i] <- coefficients(flex_grow1)[2]
  } else {
    flex.params[4,i] <- coefficients(flex_grow1)[1] + coefficients(flex_grow1)[i+1] 
    flex.params[5,i] <- coefficients(flex_grow1)[2] + coefficients(flex_grow1)[i+5]
  }  
}

# growth variability 
for(i in 1:6){
  ant.params[6,i] <- coefficients(Amer.grow.sd)[1]
  ant.params[7,i] <- coefficients(Amer.grow.sd)[2]
  gorg.params[6,i] <- coefficients(gorg.grow.sd)[1]
  gorg.params[7,i] <- coefficients(gorg.grow.sd)[2]
}  

for(i in 1:5){
  flex.params[6,i] <- coefficients(flex.grow.sd)[1]
  flex.params[7,i] <- coefficients(flex.grow.sd)[2]
} 

# Recruit size
for(i in 1:6){
  ant.params[8,i] <- coefficients(a.bull)[1]
  ant.params[9,i] <- coefficients(a.bull)[2]
  gorg.params[8,i] <- coefficients(g.bull)[1]
  gorg.params[9,i] <- coefficients(g.bull)[2]
}  

for(i in 1:5){
  flex.params[8,i] <- coefficients(e.bull)[1]
  flex.params[9,i] <- coefficients(e.bull)[2]
}  

# Recruitment
for(i in 1:6){
  ant.params[10,i] <- coefficients(Ant.dd.mod2)[1]
  ant.params[11,i] <- coefficients(Ant.dd.mod2)[2]
  ant.params[12,i] <- coefficients(Ant.dd.mod2)[3]
  ant.params[13,i] <- fitted(Ant.dd.mod2, "sigma")[1]
  ant.params[14,i] <- fitted(Ant.dd.mod2, "nu")[1]
  gorg.params[10,i] <- coefficients(gorg.dd.mod2)[1]
  gorg.params[11,i] <- coefficients(gorg.dd.mod2)[2]
  gorg.params[12,i] <- coefficients(gorg.dd.mod2)[3]
  gorg.params[13,i] <- fitted(gorg.dd.mod2, "sigma")[1]
  gorg.params[14,i] <- fitted(gorg.dd.mod2, "nu")[1]
}  

for(i in 1:5){
  flex.params[10,i] <- coefficients(flex.dd.mod2)[1]
  flex.params[11,i] <- coefficients(flex.dd.mod2)[2]
  flex.params[12,i] <- coefficients(flex.dd.mod2)[3]
  flex.params[13,i] <- fitted(flex.dd.mod2, "sigma")[1]
  flex.params[14,i] <- fitted(flex.dd.mod2, "nu")[1]
} 

# Now if tidy up the parameter matrices
rownames(ant.params) <- rownames(flex.params) <- rownames(gorg.params) <- c(
  #survival
  "surv.int","surv.slope", "surv.sd",
  #growth
  "grow.int", "grow.slope",
  #growth variability
  "grow.sd.int","grow.sd.slope",
  #recruit size
  "rcsz.par1","rcsz.par2",
  # recruitment function mean
  "mu.int","mu.slope1","mu.slope2",
  # recruitment function variability
  "sigma", 'nu',
  # Demographic ceiling parameter
  "U1")

colnames(ant.params) <- colnames(gorg.params) <- c("2013","2014","2015","2016","2017","2018")
colnames(flex.params) <- c("2014","2015","2016","2017","2018")

##############################
#Step 5: Build model functions.
##############################

# 1. Define the P kernel. Recruitment isn't built in to make a K kernel as this is a density-dependent model ---------------------------------------------------------------
mk_P_ceiling <- function(m, m.par, L, U, U1) { 
  # Building a ceiling into the model ensures that any individuals larger that U1 have the same vital rates as U1 - and is equivalent to the addition of a discrete stage. 
  # In this model U1 is equal to the largest observed size + 10% (U)
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * (outer(meshpts, pmin(meshpts, U1), P_z1z, m.par = m.par))
  return(list(meshpts = meshpts, P = P, U1 = U1)) # returning U1 here supports the implementation of a sensitivity analysis.
}

# 2. define the survival/growth kernel ------------------------------------------------------------
P_z1z = function(Size.t1, Size.t, m.par) {
  s_z(Size.t, m.par) * G_z1z(Size.t1, Size.t, m.par)
} #this is the same for all species

# 3. define the model functions for survival and growth -------------------------------------------
# Growth
G_z1z = function(Size.t1, Size.t, m.par) {
  mean <- m.par["grow.int"] + m.par["grow.slope"] * Size.t
  sd <- sqrt(pi/2) * exp((m.par["grow.sd.int"] + m.par["grow.sd.slope"] * Size.t))
  p_den_grow <- dnorm(Size.t1, mean = mean, sd = sd)
  return(p_den_grow)
} #linear growth

# Survival
s_z <- function(Size.t, m.par){
  mean <- m.par["surv.int"] + m.par["surv.slope"] * Size.t 
  sd <- m.par["surv.sd"]
  # randomly generated value from modeled distribution
  linearp <- rnorm(length(mean), mean, sd)
  # reverse logit link
  p <- 1/(1+exp(-linearp))
  # quick fix to prevent impossible probabilities
  #if(p > 1) {p = 1}
  #if(p < 0) {p = 0}
  # return output
  return(p) 
}

# 4. Define the overall recruitment function -------------------------------------------------------
R_DD <- function(Nt, m.par) {
  a <- m.par["mu.int"]
  b <- m.par["mu.slope1"]
  c <- m.par["mu.slope2"]
  linear.mu <- a + (b * Nt) + (c * (Nt^2))
  mu.mean <- exp(linear.mu) #this converts back from the log transformed regression parameters.
  # determine a randomly generated value from modeled distribution
  R <- rZAGA(n = 1, mu = mu.mean, sigma = m.par["sigma"], nu = m.par['nu'])
  # little fix to keep recruitment under biological control
  if (R < 0) {R = 0}
  return(R)
} 

# 5. Define recruit size function ------------------------------------------------------------------ 
c0_z1 = function(Size.t1, m.par) {
  shape <- m.par["rcsz.par1"]
  scale <- m.par["rcsz.par2"]                          
  p_den_rcsz <- dweibull(Size.t1, shape = shape, scale = scale) # weibull distribution.
  return(p_den_rcsz)
}

# 6.Function to iterate the IPM using the midpoint rule --------------------------------------------
Iterate <- function(nt, Nt, meshpts, P, m.par) { #this function will take the P kernel and fit the recruit size distribution to this.
  R.no <- R_DD(Nt, m.par)
  recruits <- R.no * c0_z1(meshpts,m.par) 
  return(list(main = recruits + P%*%nt, R = R.no))	
}   

##############################
#Step 6: Define final model details
##############################
#data is sqrt transformed so will not have negative values.

# define model parameters
L.list <- c(0.9*min(c(A_american$Size.t,A_american$Size.t1,sqrt(Recruits[which(Recruits$Genus == "Antillogorgia"),]$Height..cm.)), na.rm = T), 
            0.9*min(c(E_flexuosa$Size.t,E_flexuosa$Size.t1,sqrt(Recruits[which(Recruits$Genus == "Eunicea"),]$Height..cm.)), na.rm = T),
            0.9*min(c(Gorg_sp$Size.t,Gorg_sp$Size.t1,sqrt(Recruits[which(Recruits$Genus == "Gorgonia"),]$Height..cm.)), na.rm = T))
U.list <- c(1.1*max(c(A_american$Size.t,A_american$Size.t1,sqrt(Adults[which(Adults$Genus == "Antillogorgia"),]$Height.of.living.tissue..cm.)), na.rm = T), 
            1.1*max(c(E_flexuosa$Size.t,E_flexuosa$Size.t1,sqrt(Adults[which(Adults$Genus == "Eunicea"),]$Height.of.living.tissue..cm.)), na.rm = T),
            1.1*max(c(Gorg_sp$Size.t,Gorg_sp$Size.t1,sqrt(Adults[which(Adults$Genus == "Gorgonia"),]$Height.of.living.tissue..cm.)), na.rm = T))
# this ensures the max and min size ranges account for the max and mins observed across both survey types (transect and tagging), and all sites
# Store these max sizes as part of m.par
ant.params['U1',] <- max(c(A_american$Size.t,A_american$Size.t1,sqrt(Adults[which(Adults$Genus == "Antillogorgia"),]$Height.of.living.tissue..cm.)), na.rm = T)
flex.params['U1',] <- max(c(E_flexuosa$Size.t,E_flexuosa$Size.t1,sqrt(Adults[which(Adults$Genus == "Eunicea"),]$Height.of.living.tissue..cm.)), na.rm = T)
gorg.params['U1',] <- max(c(Gorg_sp$Size.t,Gorg_sp$Size.t1,sqrt(Adults[which(Adults$Genus == "Gorgonia"),]$Height.of.living.tissue..cm.)), na.rm = T)

# Now for the number of bins (m). A histogram of the relevant populations should tell us the appropriate (minimum) bin sizes.
set_graph_pars("panel1")  
#Antillogorgia
hist(initial.pop[[1]], breaks = 143)
#Flexuosa
hist(initial.pop[[2]], breaks = 117)
#Gorgonia
hist(initial.pop[[3]], breaks = 131)
# store the selected bin numbers (minimum plus a little to account for the extended size range for each integration)
# Keep bin numbers consistent across species models
m <- 200

##############################
#Step 7: Build the IPMs - moment of truth!
##############################
# recent plot panel
graphics.off()
par(mfrow = c(3,2)) #plot 6 panels at once

## Antillogorgia
Ant.2013 <- mk_P_ceiling(m = m, m.par = ant.params[,1], L = L.list[1], U = U.list[1], U1 = ant.params[15,1])
Ant.2014 <- mk_P_ceiling(m = m, m.par = ant.params[,2], L = L.list[1], U = U.list[1], U1 = ant.params[15,2])
Ant.2015 <- mk_P_ceiling(m = m, m.par = ant.params[,3], L = L.list[1], U = U.list[1], U1 = ant.params[15,3])
Ant.2016 <- mk_P_ceiling(m = m, m.par = ant.params[,4], L = L.list[1], U = U.list[1], U1 = ant.params[15,4])
Ant.2017 <- mk_P_ceiling(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1], U1 = ant.params[15,5])
Ant.2018 <- mk_P_ceiling(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1], U1 = ant.params[15,6])
# builds a kernel for each of the years separately using vital rate data corresponding to those years.

## Now plot each kernel to check its worked. 
#2013
image.plot(t(Ant.2013$P)) #remember to transpose your kernel when displaying it
abline(0,1,lwd=3, lty = "longdash")
#2014
image.plot(t(Ant.2014$P))
abline(0,1, lwd=3, lty = "longdash")
#2015
image.plot(t(Ant.2015$P))
abline(0,1, lwd=3, lty = "longdash")
#2016
image.plot(t(Ant.2016$P))
abline(0,1, lwd=3, lty = "longdash")
#2017 - Hurricane year
image.plot(t(Ant.2017$P))
abline(0,1, lwd=3, lty = "longdash")
# 2018 - Recovery year?
image.plot(t(Ant.2018$P))
abline(0,1, lwd=3, lty = "longdash")

# Just as another comparison - what do the eigenvalues look like for the P kernels of these populations (in theory 2017 and 2018 should be lowest)
# Because the models are open then lambda for the overall model is unrepresentative. Instead P provides a relative measure of the status of the existent population.
Re(eigen(Ant.2013$P)$values[1])
Re(eigen(Ant.2014$P)$values[1])
Re(eigen(Ant.2015$P)$values[1])
Re(eigen(Ant.2016$P)$values[1])
Re(eigen(Ant.2017$P)$values[1]) #and it is!!!!!
Re(eigen(Ant.2018$P)$values[1])

## Eunicea
flex.2014 <- mk_P_ceiling(m = m, m.par = flex.params[,1], L = L.list[2], U = U.list[2], U1 = flex.params[15,1])
flex.2015 <- mk_P_ceiling(m = m, m.par = flex.params[,2], L = L.list[2], U = U.list[2], U1 = flex.params[15,2])
flex.2016 <- mk_P_ceiling(m = m, m.par = flex.params[,3], L = L.list[2], U = U.list[2], U1 = flex.params[15,3])
flex.2017 <- mk_P_ceiling(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2], U1 = flex.params[15,4])
flex.2018 <- mk_P_ceiling(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2], U1 = flex.params[15,5])
# builds a kernel for each of the years separately using vital rate data corresponding to those years.

## Now plot each kernel to check its worked. 
#2014
image.plot(t(flex.2014$P))
abline(0,1, lwd=3, lty = "longdash")
#2015
image.plot(t(flex.2015$P))
abline(0,1, lwd=3, lty = "longdash")
#2016
image.plot(t(flex.2016$P))
abline(0,1, lwd=3, lty = "longdash")
#2017 - Hurricane year
image.plot(t(flex.2017$P))
abline(0,1, lwd=3, lty = "longdash")
# 2018 - Recovery year?
image.plot(t(flex.2018$P))
abline(0,1, lwd=3, lty = "longdash")
# recruitment is much bigger in flexuosa

# Just as another comparison - what do the eigenvalues look like for the P kernels of these populations (in theory 2017 and 2018 should be lowest)
# Because the models are open then lambda for the overall model is unrepresentative. Instead P provides a relative measure of the status of the existant population.
Re(eigen(flex.2014$P)$values[1])
Re(eigen(flex.2015$P)$values[1])
Re(eigen(flex.2016$P)$values[1])
Re(eigen(flex.2017$P)$values[1])
Re(eigen(flex.2018$P)$values[1]) # viability of the adult stock is no where nearly as badly hit by the hurricane.

## Gorgonia
gorg.2013 <- mk_P_ceiling(m = m, m.par = gorg.params[,1], L = L.list[3], U = U.list[3], U1 = gorg.params[15,1])
gorg.2014 <- mk_P_ceiling(m = m, m.par = gorg.params[,2], L = L.list[3], U = U.list[3], U1 = gorg.params[15,2])
gorg.2015 <- mk_P_ceiling(m = m, m.par = gorg.params[,3], L = L.list[3], U = U.list[3], U1 = gorg.params[15,3])
gorg.2016 <- mk_P_ceiling(m = m, m.par = gorg.params[,4], L = L.list[3], U = U.list[3], U1 = gorg.params[15,4])
gorg.2017 <- mk_P_ceiling(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3], U1 = gorg.params[15,5])
gorg.2018 <- mk_P_ceiling(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3], U1 = gorg.params[15,6])
# builds a kernel for each of the years separately using vital rate data corresponding to those years.

## Now plot each kernel to check its worked. 
#2013
image.plot(t(gorg.2013$P))
abline(0,1, lwd=3, lty = "longdash")
#2014
image.plot(t(gorg.2014$P))
abline(0,1, lwd=3, lty = "longdash")
#2015
image.plot(t(gorg.2015$P))
abline(0,1, lwd=3, lty = "longdash")
#2016
image.plot(t(gorg.2016$P))
abline(0,1, lwd=3, lty = "longdash")
#2017 - Hurricane year
image.plot(t(gorg.2017$P))
abline(0,1, lwd=3, lty = "longdash")
# 2018 - Recovery year?
image.plot(t(gorg.2018$P))
abline(0,1, lwd=3, lty = "longdash")
# recruitment is very low in gorgonia 

# Just as another comparison - what do the eigenvalues look like for the P kernels of these populations (in theory 2017 and 2018 should be lowest)
# Because the models are open then lambda for the overall model is unrepresentative. Instead P provides a relative measure of the status of the existant population.
Re(eigen(gorg.2013$P)$values[1])
Re(eigen(gorg.2014$P)$values[1])
Re(eigen(gorg.2015$P)$values[1])
Re(eigen(gorg.2016$P)$values[1])
Re(eigen(gorg.2017$P)$values[1])
Re(eigen(gorg.2018$P)$values[1]) # Gorgonia appears to decline prior to the hurricane impact.

######################### End of code #############################
