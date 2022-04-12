# This script is for working with the Caribbean octocoral demography data to construct density dependent IPMs.
# This script goes through the calculation of vital rates (growth survival and recruitment), model parameterisation and
# determination of adult viability (This is done simultaneously for all species).

# Last modified: Mar 2022
# Primary Author: James Cant

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

# set working directory
setwd("/Users/james/Documents/Gorgonian")

#load in the data set
soft_corals <- read.csv("IPM_data_2013_2019_FINAL.csv", stringsAsFactors = TRUE)
source("Standard graphical Pars.R")
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
plot(Size.t1~ Size.t, data = Gorg_sp, col = "Black", #begin by plotting data
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

#***********    
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
plot(Surv~Size.t, type = "n", data = Gorg_sp, #now we can plot the data seperatly - start by plotting a blank plot
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
# this part of the script involves data from Privitera-Johnson et al 2015, and unpublished data

# For Antillogorgia and Gorgonia use the overall community data set listing the number of recruits and adults recorded in the community between 2013 and 2018.
DD1 <- read.csv("Overall community data.csv")
# this data file contains a lot of detail that I wont be using. I just need the size of recruits. 
# before I can use this data it needs to be formated so that it only contains the species I am interested in and that height is on the same scale as the tagged data set.
DD1 <- subset(DD1, Height != "None" & #remove the 'none' height entries.
                Season != "Post-Irma" & # the post Irma surveys are out of the annual cycle.
                Genus %in% c("Antillogorgia", "Eunicea", "Gorgonia")) #these are the only species I am interested in from this datafile
# Now according to Lorenzo I am only working with a few species within the genus' above
DD1 <- subset(DD1, Species %in% c("aa","asp","gv", "efl", "esp"))
DD1 <- subset(DD1, !is.na(DD1$Height)) #remove any NAs
#convert the height data into the same format as that used for the growth and survival analysis.
summary(DD1$Height)
# I need to fix the entries that are not numerical values and then change the format of the data.
DD1$Height[which(DD1$Height == "<0.1")] <- 0.09
DD1$Height[which(DD1$Height == "<0.2")] <- 0.19
DD1$Height[which(DD1$Height == "<0.3")] <- 0.29
DD1$Height[which(DD1$Height == "<0.4")] <- 0.39
DD1$Height[which(DD1$Height == "<0.5")] <- 0.49
DD1$Height[which(DD1$Height == "<0.6")] <- 0.59
DD1$Height <- as.numeric(paste(DD1$Height))
# the data can now be transformed - by what every transformation was used in the tagged colony data
hist(sqrt(DD1$Height)) # this is heavily skewed but we aren't analyzing this so that is fine.
DD1$Height <- sqrt(DD1$Height)

# now break it down into the two genus groups (and remove any none recruits) (note: there is a better data frame for Eunicea recruit size data)
ant.recruits <- subset(DD1, Age.Use == "Recruit" & Genus == "Antillogorgia")
gorg.recruits <- subset(DD1, Age.Use == "Recruit" & Genus == "Gorgonia")

# now extract the same data for Eunicea - this is located within the data_recruits_size_eunicea file.
# this file again has a lot of detail I don't need
eunicea.recruits <- read.csv("data_recruits_size_eunicea.csv")
# again I need to remove entries that are not applicable
eunicea.recruits <- subset(eunicea.recruits, Height..cm. != "nd" & #remove the 'none' height entries.
                             Genus == "Eunicea" &
                             Species.code %in% c("efl", "esp")) # these are the only species that can be attributed to Eunicea flexuosa
# check the Height data
summary(eunicea.recruits$Height..cm.) 
#some entries need altering before these can be made numeric.
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.1")] <- 0.09
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.2")] <- 0.19
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.3")] <- 0.29
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.4")] <- 0.39
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.5")] <- 0.49
eunicea.recruits$Height..cm.[which(eunicea.recruits$Height..cm. == "<0.6")] <- 0.59
# now the data can be converted to numeric and formatted
eunicea.recruits$Height..cm. <- as.numeric(paste(eunicea.recruits$Height..cm.))
eunicea.recruits <- subset(eunicea.recruits, !is.na(eunicea.recruits$Height..cm.)) #remove any NAs
eunicea.recruits$Height..cm. <- sqrt(eunicea.recruits$Height..cm.)

# I will also use the Overall community data to extract the current populations for projection.
# this can be based off of the 2018 surveyed populations (it appears a full community analysis wasn't conducted in 2019)
gorg.2018 <- DD1[which(DD1$Year == 2018 & DD1$Genus == "Gorgonia" & DD1$Site %in% c("E Cabritte", "Europa")),] # the tagged colonies are only in E Cabritte and Europa.
ant.2018 <- DD1[which(DD1$Year == 2018 & DD1$Genus == "Antillogorgia" & DD1$Site %in% c("E Cabritte", "Europa")),]
flex.2018 <- DD1[which(DD1$Year == 2018 & DD1$Genus == "Eunicea" & DD1$Site %in% c("E Cabritte", "Europa")),]
initial.pop <- list(ant.2018$Height, flex.2018$Height, gorg.2018$Height) #store them

# Now I have sorted the size data I can determine the recruit size distributions for these three octocoral species.
## Gorgonia
hist(gorg.recruits$Height, breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,200),
     xaxs = "i", yaxs = "i") #could be a normal distribution - but it is continuous and none negative.
# the fitdistrplus package contains functions for determining these types of distributions 
# the package allows for the plotting of a Cullen and Frey graph to determine the distribution of the data.
descdist(gorg.recruits$Height, boot = 1000) #according to this the data is very centred and close to a normal and uniform distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
g.norm <- fitdist(gorg.recruits$Height, distr = "norm", method = "mle") #normal distribution
g.uni <- fitdist(gorg.recruits$Height, distr = "unif", method = "mle") #uniform distribution
denscomp(list(g.norm,g.uni)) # visually the normal distribution is the best!

## A. american
hist(ant.recruits$Height, breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,120),
     xaxs = "i", yaxs = "i") # this one is skewed towards larger recruit sizes.
# check the distribution.
descdist(ant.recruits$Height, boot = 1000) #according to this the data is very centred and close to a normal and uniform distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
a.norm <- fitdist(ant.recruits$Height, distr = "norm", method = "mle") #normal distribution
# I will also test the fits of more skewed distributions
a.gam <- fitdist(ant.recruits$Height, distr = "gamma", method = "mle") #gamma distribution
a.bull <- fitdist(ant.recruits$Height, distr = "weibull", method = "mle") # weibull
a.log <- fitdist(ant.recruits$Height, distr = "lnorm", method = "mle") #log normal
# plot them all together to see which fits best
denscomp(list(a.norm,a.gam, a.bull, a.log)) #visually the normal and weibull distributions are the best!
#compare model AICs and BICs
summary(a.norm)$aic; summary(a.gam)$aic; summary(a.bull)$aic; summary(a.log)$aic # Weibull is the best mathematical fit
summary(a.norm)$bic; summary(a.gam)$bic; summary(a.bull)$bic; summary(a.log)$bic

## E. flexuosa
hist(eunicea.recruits$Height, breaks = 10,
     xlab = "",
     ylab= "",
     main = NULL,
     cex.axis = 1.5, 
     xlim = c(0,2.5),
     ylim = c(0,160),
     xaxs = "i", yaxs = "i") # doesn't look like  a normal distribution - this one is heavily skewed towards smaller recruit sizes.
# check the distribution.
descdist(eunicea.recruits$Height, boot = 1000) #according to this the data is very centred and close to a normal and uniform distribution
# this plot also indicates that the data is beta distributed but this only applies to probability values.
e.norm <- fitdist(eunicea.recruits$Height, distr = "norm", method = "mle") #normal distribution
# I will also test the fits of more skewed distributions
e.gam <- fitdist(eunicea.recruits$Height, distr = "gamma", method = "mle") #gamma distribution
e.bull <- fitdist(eunicea.recruits$Height, distr = "weibull", method = "mle") # weibull
e.log <- fitdist(eunicea.recruits$Height, distr = "lnorm", method = "mle") #log normal
# plot them all together to see which fits best
denscomp(list(e.norm,e.gam, e.bull, e.log)) #visually the lognormal and gamma distributions are the best!
#compare model AICs and BICs
summary(e.norm)$aic; summary(e.gam)$aic; summary(e.bull)$aic; summary(e.log)$aic # gamma fits the best mathematically
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
denscomp(e.gam, #what density is being plotted with the data
         xlim = c(0,2.5),
         ylim = c(0,1.2),
         main = NULL,
         xlab = NULL,
         ylab = NULL,
         addlegend = FALSE,
         xaxs = "i", yaxs = "i",
         cex.axis = 1.5)

# Gorgonia
denscomp(g.norm, #what density is being plotted with the data
         xlim = c(0,2.5),
         ylim = c(0,1.2),
         main = NULL,
         xlab = NULL,
         ylab = NULL,
         addlegend = FALSE,
         xaxs = "i", yaxs = "i",
         cex.axis = 1.5)

##############################
#Step 4: Determine the density-dependent recruitment function.
##############################
# Load the data file Recruit_Adult density data
DD2 <- read.csv("Recruit_Adult density Data.csv")

# Correct the dataframe formatting
DD2$Year <- as.factor(DD2$Year)
DD2$G.Juveniles <- as.numeric(DD2$G.Juveniles)
DD2$G.Adults <- as.numeric(DD2$G.Adults)
DD2$A.Juveniles <- as.numeric(DD2$A.Juveniles)
DD2$A.Adults <- as.numeric(DD2$A.Adults)
DD2$E.Juveniles <- as.numeric(DD2$E.Juveniles)
DD2$E.Adults <- as.numeric(DD2$E.Adults)

# Now I need to adjust the scale of the counts. Counts were made in m2 quadrats so I need to scale this up to reflect the 120m2 area from which the community data is being pooled.
Adjusted.data <- data.frame(Site = DD2$Site,
                            Year = DD2$Year,
                            G.Juveniles = DD2$G.Juveniles*120,
                            G.Adults = DD2$G.Adults*120,
                            A.Juveniles = DD2$A.Juveniles*120,
                            A.Adults = DD2$A.Adults*120,
                            E.Juveniles = DD2$E.Juveniles*120,
                            E.Adults = DD2$E.Adults*120)
Adjusted.data$total.Adults <- (Adjusted.data$G.Adults + Adjusted.data$A.Adults + Adjusted.data$E.Adults)

## Before moving on with the analysis of the density-dependent recruitment function there are a few concepts to iron out.
# 1) In their 2015 paper, Privitera-Johnson et al. used site means as the dependent variables with the different sites representing the replicates, and estimated the relationship between juvenile and adult densities using a type 2 regression approach.
# This was to reflect the error present in both counts of adult and juvenile colonies. However type-2 regression is not considered as effective when using x variables to predict y variables. Thus it is frequently argued that OLS regression is the prefer approach in prediction exercises.
# 2) Recruit densities are a form of count data. i.e. discrete non-negative values. Therefore the relationship between adult density and recruit density will be determined using poisson regression.
# 3) The juvenile count data is zero inflated. There are zeros arising because the (1) conditions at the time prevented juvenile settlement (effective zero), and there were no juvenilles to settle (false zero).
# The approach used here therefore will be a zero-inflated poisson model.

# Now determine whether patterns in juvenile densities for each species correspond more with con-specific adult densities or total adult densities
# This test will be carried out using distance correlation which evaluates for both linear and non-linear associations between variables

# Antilogorgia
dcor(Adjusted.data$A.Juveniles, Adjusted.data$A.Adults)
dcor(Adjusted.data$A.Juveniles, Adjusted.data$total.Adults)
# Eunicea
dcor(Adjusted.data$E.Juveniles, Adjusted.data$E.Adults)
dcor(Adjusted.data$E.Juveniles, Adjusted.data$total.Adults)
# Gorgonia
dcor(Adjusted.data$G.Juveniles, Adjusted.data$G.Adults)
dcor(Adjusted.data$G.Juveniles, Adjusted.data$total.Adults)
# Across all three species juvenile counts correlated better with total adult density.

# Perform zero-inflated regression analysis.
## Antilogorgia
# determine relationship
Ant.rec.mod <- gamlss(A.Juveniles ~ total.Adults, data = Adjusted.data, family = "ZIP") # zero inflated poisson approach
# estimate predicted pattern
Ant.rec <- ggpredict(Ant.rec.mod, terms = c("total.Adults[all]"))

## Eunicea
# determine relationship
Flex.rec.mod <- gamlss(E.Juveniles ~ total.Adults, data = Adjusted.data, family = "ZIP",
                       control = gamlss.control(n.cyc = 50)) # zero inflated poisson approach
# estimate predicted pattern
Flex.rec <- ggpredict(Flex.rec.mod, terms = c("total.Adults[all]"))

## Gorgonia
# determine relationship
Gorg.rec.mod <- gamlss(G.Juveniles ~ total.Adults, data = Adjusted.data, family = "ZIP") # zero inflated poisson approach
# estimate predicted pattern
Gorg.rec <- ggpredict(Gorg.rec.mod, terms = c("total.Adults[all]"))

# All species return an increasing number of juveniles with increasing adult density - this fits evidence of a positive influence on larvae settlement generated by the presence of higher adult densities.
# These patterns will not need plotting - they are visually uninformative

### But what about the variance in recruit density with adult density
## Antilogorgia
# view the residuals of the recruitment model
plot(Adjusted.data$total.Adults, abs(resid(Ant.rec.mod)), xlab='Adult density',ylab='residual') 
Ant.DD.sd <- glm(abs(resid(Ant.rec.mod))~Adjusted.data$total.Adults, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
Ant.DD.sd2 = lm(abs(resid(Ant.rec.mod))~Adjusted.data$total.Adults) 
# which model is the best?
AIC(Ant.DD.sd,Ant.DD.sd2) #according to the AICs the gamma model fits better
BIC(Ant.DD.sd,Ant.DD.sd2) # that is agreed by the BIC.
# how does the model fit the data visually?
lines(spline(Adjusted.data$total.Adults, predict(Ant.DD.sd, type = "response"))) 

## Eunicea
# view the residuals of the recruitment model
plot(Adjusted.data$total.Adults, abs(resid(Flex.rec.mod)),xlab='Adult density',ylab='residual') 
flex.DD.sd <- glm(abs(resid(Flex.rec.mod))~Adjusted.data$total.Adults, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
flex.DD.sd2 = lm(abs(resid(Flex.rec.mod))~Adjusted.data$total.Adults) 
# which model is the best?
AIC(flex.DD.sd,flex.DD.sd2) #according to the AICs the gamma model fits better
BIC(flex.DD.sd,flex.DD.sd2) # that is agreed by the BIC.
# how does the model fit the data visually?
lines(spline(Adjusted.data$total.Adults, predict(flex.DD.sd, type = "response"))) 

## Gorgonia
# view the residuals of the recruitment model
plot(Adjusted.data$total.Adults,abs(resid(Gorg.rec.mod)),xlab='Adult density',ylab='residual') 
gorg.DD.sd <- glm(abs(resid(Gorg.rec.mod))~Adjusted.data$total.Adults, family = Gamma(link = "log")) #gamma will prevent variance from dropping below zero - also allows for a none-linear trend.
#but is the relationship linear?
gorg.DD.sd2 = lm(abs(resid(Gorg.rec.mod))~Adjusted.data$total.Adults) 
# which model is the best?
AIC(gorg.DD.sd,gorg.DD.sd2) #according to the AICs the gamma model fits better
BIC(gorg.DD.sd,gorg.DD.sd2) # that is agreed by the BIC.
# how does the model fit the data visually?
lines(spline(Adjusted.data$total.Adults, predict(gorg.DD.sd, type = "response"))) 

##############################
#Step 5: Store model coefficients for construction of IPM.
##############################

# blank matrices to store the yearly parameters for each species separately
gorg.params <- matrix(NA, 14, 6)
flex.params <- matrix(NA, 14, 5)
ant.params <- matrix(NA, 14, 6)

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
  gorg.params[8,i] <- coefficients(g.norm)[1]
  gorg.params[9,i] <- coefficients(g.norm)[2]
}  

for(i in 1:5){
  flex.params[8,i] <- coefficients(e.gam)[1]
  flex.params[9,i] <- coefficients(e.gam)[2]
}  

# Recruitment
for(i in 1:6){
  ant.params[10,i] <- coefficients(Ant.rec.mod)[1]
  ant.params[11,i] <- coefficients(Ant.rec.mod)[2]
  ant.params[12,i] <- fitted(Ant.rec.mod, "sigma")[1]
  ant.params[13,i] <- coefficients(Ant.DD.sd)[1]
  ant.params[14,i] <- coefficients(Ant.DD.sd)[2]
  gorg.params[10,i] <- coefficients(Gorg.rec.mod)[1]
  gorg.params[11,i] <- coefficients(Gorg.rec.mod)[2]
  gorg.params[12,i] <- fitted(Gorg.rec.mod, "sigma")[1]
  gorg.params[13,i] <- coefficients(gorg.DD.sd)[1]
  gorg.params[14,i] <- coefficients(gorg.DD.sd)[2]
}  

for(i in 1:5){
  flex.params[10,i] <- coefficients(Flex.rec.mod)[1]
  flex.params[11,i] <- coefficients(Flex.rec.mod)[2]
  flex.params[12,i] <- fitted(Flex.rec.mod, "sigma")[1]
  flex.params[13,i] <- coefficients(flex.DD.sd)[1]
  flex.params[14,i] <- coefficients(flex.DD.sd)[2]
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
  "mu.int","mu.slope","sigma",
  # recruitment function boundaries
  "rec.sd.int","rec.sd.slope")

colnames(ant.params) <- colnames(gorg.params) <- c("2013","2014","2015","2016","2017","2018")
colnames(flex.params) <- c("2014","2015","2016","2017","2018")

##############################
#Step 5: Build model functions.
##############################

# 1. Define the P kernel. Recruitment isn't built in to make a K kernel as this is a density-dependent model ---------------------------------------------------------------
mk_P <- function(m, m.par, L, U) {
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
  return(list(meshpts = meshpts, P = P))
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
  b <- m.par["mu.slope"]
  linear.mu <- a + (b * Nt) 
  mu.mean <- exp(linear.mu) #this converts back from the log transformed regression parameters.
  mu.sd <- sqrt(pi/2) * exp((m.par["rec.sd.int"] + m.par["rec.sd.slope"] * Nt)) # determine variance in recruitment mean
  # randomly generated value from modeled distribution
  mu <- rnorm(1, mu.mean, mu.sd)
  # The function now needs to determine if recruitment is zero or not.
  R <- rZIP(n = 1, mu = mu, sigma = m.par["sigma"])
  # little fix to keep recruitment under biological control
  if (R < 0) {R = 0}
  return(R)
} 

# 5. Define recruit size function ------------------------------------------------------------------ 
c0_z1.ant = function(Size.t1, m.par) {
  shape <- m.par["rcsz.par1"]
  scale <- m.par["rcsz.par2"]                          
  p_den_rcsz <- dweibull(Size.t1, shape = shape, scale = scale) # weibull distribution.
  return(p_den_rcsz)
}

# recruit size
c0_z1.flex = function(Size.t1, m.par) {
  shape <- m.par["rcsz.par1"]
  rate <- m.par["rcsz.par2"]                          
  p_den_rcsz <- dgamma(Size.t1, shape = shape, rate = rate) # gamma distribution.
  return(p_den_rcsz)
}

# recruit size
c0_z1.gorg = function(Size.t1, m.par) {
  mean <- m.par["rcsz.par1"]
  sd <- m.par["rcsz.par2"]                          
  p_den_rcsz <- dnorm(Size.t1, mean = mean, sd = sd) # normal distribution.
  return(p_den_rcsz)
}

# 6.Function to iterate the IPM using the midpoint rule
Iterate.a <- function(nt, Nt, meshpts, P, m.par) { #this function will take the P kernel and fit the recruit size distribution to this.
  h <- meshpts[2]-meshpts[1]
  R.no <- R_DD(Nt, m.par)
  recruits <- R.no * c0_z1.ant(meshpts,m.par) 
  return(recruits + P%*%nt)	
}    # Antillogorgia

Iterate.e <- function(nt, Nt, meshpts, P, m.par) { #this function will take the P kernel and fit the recruit size distribution to this.
  h <- meshpts[2]-meshpts[1]
  R.no <- R_DD(Nt, m.par)
  recruits <- R.no * c0_z1.flex(meshpts,m.par) 
  return(recruits + P%*%nt)	
}    # Euniciea

Iterate.g <- function(nt, Nt, meshpts, P, m.par) { #this function will take the P kernel and fit the recruit size distribution to this.
  h <- meshpts[2]-meshpts[1]
  R.no <- R_DD(Nt, m.par)
  recruits <- R.no * c0_z1.gorg(meshpts,m.par) 
  return(recruits + P%*%nt)	
}    # gorgonia

##############################
#Step 6: Define final model details
##############################
#data is sqrt transformed so will not have negative values.

# define model parameters
L.list <- c(0.9*min(c(A_american$Size.t,A_american$Size.t1,DD1[which(DD1$Genus == "Antillogorgia"),]$Height), na.rm = T), 
            0.9*min(c(E_flexuosa$Size.t,E_flexuosa$Size.t1,DD1[which(DD1$Genus == "Eunicea"),]$Height), na.rm = T),
            0.9*min(c(Gorg_sp$Size.t,Gorg_sp$Size.t1,DD1[which(DD1$Genus == "Gorgonia"),]$Height), na.rm = T))
U.list <- c(1.1*max(c(A_american$Size.t,A_american$Size.t1,DD1[which(DD1$Genus == "Antillogorgia"),]$Height), na.rm = T), 
            1.1*max(c(E_flexuosa$Size.t,E_flexuosa$Size.t1,DD1[which(DD1$Genus == "Eunicea"),]$Height), na.rm = T),
            1.1*max(c(Gorg_sp$Size.t,Gorg_sp$Size.t1,DD1[which(DD1$Genus == "Gorgonia"),]$Height), na.rm = T))
# this ensures the max and min size ranges account for the max and mins observed across both survey types (transect and tagging), and all sites

# Now for the number of bins (m). A histogram of the relevant populations should tell us the appropriate (minimum) bin sizes.
set_graph_pars("panel1")  
#Antillogorgia
hist(initial.pop[[1]], breaks = 143)
#Flexuosa
hist(initial.pop[[2]], breaks = 117)
#Gorgonia
hist(initial.pop[[3]], breaks = 117)
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
Ant.2013 <- mk_P(m = m, m.par = ant.params[,1], L = L.list[1], U = U.list[1])
Ant.2014 <- mk_P(m = m, m.par = ant.params[,2], L = L.list[1], U = U.list[1])
Ant.2015 <- mk_P(m = m, m.par = ant.params[,3], L = L.list[1], U = U.list[1])
Ant.2016 <- mk_P(m = m, m.par = ant.params[,4], L = L.list[1], U = U.list[1])
Ant.2017 <- mk_P(m = m, m.par = ant.params[,5], L = L.list[1], U = U.list[1])
Ant.2018 <- mk_P(m = m, m.par = ant.params[,6], L = L.list[1], U = U.list[1])
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
flex.2014 <- mk_P(m = m, m.par = flex.params[,1], L = L.list[2], U = U.list[2])
flex.2015 <- mk_P(m = m, m.par = flex.params[,2], L = L.list[2], U = U.list[2])
flex.2016 <- mk_P(m = m, m.par = flex.params[,3], L = L.list[2], U = U.list[2])
flex.2017 <- mk_P(m = m, m.par = flex.params[,4], L = L.list[2], U = U.list[2])
flex.2018 <- mk_P(m = m, m.par = flex.params[,5], L = L.list[2], U = U.list[2])
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
gorg.2013 <- mk_P(m = m, m.par = gorg.params[,1], L = L.list[3], U = U.list[3])
gorg.2014 <- mk_P(m = m, m.par = gorg.params[,2], L = L.list[3], U = U.list[3])
gorg.2015 <- mk_P(m = m, m.par = gorg.params[,3], L = L.list[3], U = U.list[3])
gorg.2016 <- mk_P(m = m, m.par = gorg.params[,4], L = L.list[3], U = U.list[3])
gorg.2017 <- mk_P(m = m, m.par = gorg.params[,5], L = L.list[3], U = U.list[3])
gorg.2018 <- mk_P(m = m, m.par = gorg.params[,6], L = L.list[3], U = U.list[3])
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

########################## End of code #############################