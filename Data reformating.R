# This script is for transposing the 'Raw' octocoral dataframe to convert each longitudinal record into a stacked datafile 
# consisting of individual rows documenting the separate annual transitions of tagged colonies.
# This script will also determine size and year at death (if applicable) for each individual

# IMPORTANT: this script will return two output files:
   # (1) a temporal record of individual survival and size (survival analysis)
   # (2) a file documenting the initial and final size of tagged individuals followed over annual intervals (IPM analysis)

# Primary Author: James Cant 
# Email: james.cant91@gmail.com
# ---------------------------------------------------------------------------------

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(plyr)

# set working directory
setwd("DATA FILE DIRECTORY PATH")

# Are you checking outliers?
#Outliers = TRUE
Outliers = FALSE

# load in the required data (na.strings fills empty cells) - the line of code used here depends on whether this script is being used to process the data with or without outliers checked.
if(Outliers == FALSE){
  Raw_data <- read.csv("DATA FILE", na.strings=c("", "NA"), numerals = "no.loss") 
} else {
  Raw_data <- read.csv("DATA FILE", na.strings=c("", "NA"), numerals = "no.loss") 
}

####################################################
# STEP 1: Tidy Data & standardize column notation
####################################################
# In this script the only manual effort is defining the columns and row indexing corresponding 
# with the species and site records within the .csv file

# Split data file into individual site/species files
# EU = Europa, EC = East Cabritte
# AA = Antilogorgia, EFL = Eunicea, GV = Gorgonia
EU_AA_R <- Raw_data[1:101, 1:11]
EU_EFL_R <- Raw_data[1:151, 13:22]
EU_GV_R <- Raw_data[1:123, 24:41]
EC_AA_R <- Raw_data[1:146, 43:52]
EC_EFL_R <- Raw_data[1:156, 54:63]
EC_GV_R <- Raw_data[1:80, 65:78]

# Reassign column names
colnames(EU_AA_R) <- c("Site","Species","ID","t2013","t2014","t2015","t2016","t2017","t2017.5","t2018","t2019")
colnames(EU_EFL_R) <- c("Site","Species","ID","t2014","t2015","t2016","t2017","t2017.5","t2018","t2019")
colnames(EU_GV_R) <- c("Site","Species","ID","t2013","t2013w","t2014","t2014w","t2015","t2015w","t2016","t2016w","t2017","t2017w","t2017.5","t2018","t2018w","t2019","t2019w")
colnames(EC_AA_R) <- c("Site","Species","ID","t2014","t2015","t2016","t2017","t2017.5","t2018","t2019")
colnames(EC_EFL_R) <- c("Site","Species","ID","t2014","t2015","t2016","t2017","t2017.5","t2018","t2019")
colnames(EC_GV_R) <- c("Site","Species","ID","t2015","t2015w","t2016","t2016w","t2017","t2017w","t2017.5","t2018","t2018w","t2019","t2019w")

# Add columns corresponding with missing years into each file
EU_EFL_R$t2013 <- NA
EC_EFL_R$t2013 <- NA
EC_AA_R$t2013 <- NA
EC_GV_R$t2013 <- EC_GV_R$t2013w <- EC_GV_R$t2014 <- EC_GV_R$t2014w <- NA

# And reorder columns to resolve annual sequence
# define order sequence
order1 <- c("Site","Species","ID","t2013","t2014","t2015","t2016","t2017","t2017.5","t2018","t2019")
order2 <- c("Site","Species","ID","t2013","t2013w","t2014","t2014w","t2015","t2015w","t2016","t2016w","t2017","t2017w","t2017.5","t2018","t2018w","t2019","t2019w")
# re-order
EU_AA_R <- EU_AA_R[order1]
EU_EFL_R <- EU_EFL_R[order1]
EU_GV_R <- EU_GV_R[order2]
EC_AA_R <- EC_AA_R[order1]
EC_EFL_R <- EC_EFL_R[order1]
EC_GV_R <- EC_GV_R[order2]

####################################################
# STEP 2: Transpose each separate species/site file and stack back together (ready for temporal survival analysis)
####################################################
# This section will also determine each individuals size at death, and their average lifetime recorded growth rate

# bring together the datafiles so that the transposing can be applied using a loop
R_data_list <- list(EU_AA_R, EU_EFL_R, EU_GV_R, EC_AA_R, EC_EFL_R, EC_GV_R)
# define output storage
R_data_trans <- list()

# run loop
for(ii in 1:length(R_data_list)){
  # reshape the data from its current wide format into a long (stacked) format
  if(ii %in% c(1,2,4,5)) { 
    data.temp <- reshape(R_data_list[[ii]], 
                         direction = "long",
                         idvar=c("Site","Species","ID"),
                         varying = c("t2013", "t2014","t2015","t2016","t2017","t2017.5","t2018","t2019"), 
                         v.names = "Size", timevar = "Year")
  } else {
    # the approach is slightly different for Gorgonia ventalina as we are only interested in the height data
    colNamesToDrop <- c("t2013w","t2014w","t2015w","t2016w","t2017w","t2018w","t2019w") # single out width columns for removal
    colNamesToKeep <- ! names(R_data_list[[ii]]) %in% colNamesToDrop
    data.temp <- reshape(R_data_list[[ii]][, colNamesToKeep], 
                         direction = "long",
                         idvar=c("Site","Species","ID"),
                         varying = c("t2013", "t2014","t2015","t2016","t2017","t2017.5","t2018","t2019"), 
                         v.names = "Size", timevar = "Year")
  }
  
  # bring together measurements from the same individuals
  data.temp <- data.temp[order(data.temp$ID, data.temp$Year),]
  
  # Edit the year column to improve clarity
  data.temp$Year <- rep(c("2013","2014","2015","2016","2017","2017.5","2018","2019"), length.out = length(data.temp$Year))
  
  # Add mortality column (it is necessary to define mortality in a separate column in order to return the size data to its numeric form)
  data.temp$Dead <- NA
  data.temp[!is.na(data.temp$Size),]$Dead <- "L"
  data.temp[which(data.temp$Size == "DEAD"),]$Dead <- "D"
  #data.temp[which(data.temp$Size == "NFX"),]$Dead <- NA # It is preferred to keep colonies listed as alive even they are missed during surveys but rediscovered later.
  
  # return size data to its numeric form
  data.temp$Size <- suppressWarnings(as.numeric(data.temp$Size)) # This will return a warning as NFX and DEAD entries are converted to NA. 
                                                                 # This is fine as these have already been store elsewhere
  
  # Calculate size at death
  data.temp$Size.death <- NA
  # loop through each entry, determine if it corresponds with a death, and if it does store the previous size record
  for(xx in 1:dim(data.temp)[1]){
    if(data.temp$Dead[xx] == "D" & !is.na(data.temp$Dead[xx])) { data.temp$Size.death[xx] <- data.temp$Size[xx-1]}
  }
  
  # determine growth trajectories
  # annual rates
  data.temp$growth.rate <- NA
  for (jj in 1:(dim(data.temp)[1]-1)) {
    if(data.temp$ID[jj] == data.temp$ID[jj+1]) { data.temp$growth.rate[jj] <- data.temp$Size[jj+1] - data.temp$Size[jj] }
  }
  # average observed across the full record of the tagged individual
  data.temp$avg.growth <- NA
  for (ss in 1:dim(data.temp)[1]) {
    if(!is.na(data.temp$Size.death[ss])) { data.temp$avg.growth[ss] <- mean(data.temp[which(data.temp$ID == data.temp$ID[ss]),]$growth.rate, na.rm = TRUE) }
  }
  
  # Finally, add other details associated with the data
  # Month timeline
  data.temp$Months <- rep(c(0,12,24,36,48,51,60,72), length.out = dim(data.temp)[1])
  # Pre- or Post-hurricane survey?
  data.temp$Hurricane <- ifelse(data.temp$Year %in% c("2013","2014","2015","2016","2017"), "Pre", "Post")
  
  # Reorder columns to tidy up!
  order <- c("Site","Species","ID","Year","Months","Hurricane","Size","Dead","Size.death","growth.rate","avg.growth")
  data.temp <- data.temp[order]
  # And store the newly transposed data
  R_data_trans[[ii]] <- data.temp
}

# Stitch together the stack files for each of the different species/site combinations into one single data file
New_data <- ldply(R_data_trans, data.frame)
# Ensure formatting across each variable
# Factor variables
New_data$Site <- as.factor(New_data$Site)
New_data$Species <- as.factor(New_data$Species)
New_data$ID <- as.factor(New_data$ID)
New_data$Year <- as.factor(New_data$Year)
New_data$Hurricane <- as.factor(New_data$Hurricane)
New_data$Dead <- as.factor(New_data$Dead)

# and save the file

write.csv(New_data, "NEW DATA FILE", row.names = FALSE)

####################################################
# STEP 3: Transpose each separate species/site file and stack back together (ready for IPM analysis)
####################################################
# For the IPM the data needs stacking but with two size columns per row - one defining size at t and the other size at t+1
# Redefine new storage
R_data_IPM <- list()

# Re-run loop
for(ii in 1:length(R_data_list)){
  # first determine the data being processed
  data.temp <- R_data_list[[ii]]
  
  # For the IPMs we are only interested in full annual intervals - therefore the data collected in Nov 2017 - needs removing.
  # To do this, any records of DEAD need to be shifted into the corresponding 2018 cell first.
  data.temp[which(data.temp$t2017.5 == "DEAD"),]$t2018 <- "DEAD"
  # and remove the intermediate survey data 
  data.temp$t2017.5 <- NULL
 
  # Now reshape the data from its current wide format into a long (stacked) format
  if(ii %in% c(1,2,4,5)) { 
    data.temp <- reshape(data.temp, 
                         direction = "long",
                         idvar=c("Site","Species","ID"),
                         varying = c("t2013", "t2014","t2015","t2016","t2017","t2018","t2019"), 
                         v.names = "Size.t", timevar = "Year.t")
  } else {
    # the approach is slightly different for Gorgonia ventalina as we are only interested in the height data
    colNamesToDrop <- c("t2013w","t2014w","t2015w","t2016w","t2017w","t2018w","t2019w") # single out width columns for removal
    colNamesToKeep <- ! names(data.temp) %in% colNamesToDrop
    data.temp <- reshape(data.temp[, colNamesToKeep], 
                         direction = "long",
                         idvar=c("Site","Species","ID"),
                         varying = c("t2013", "t2014","t2015","t2016","t2017","t2018","t2019"), 
                         v.names = "Size.t", timevar = "Year.t")
  }
  
  # bring together measurements from the same individuals
  data.temp <- data.temp[order(data.temp$ID, data.temp$Year.t),]
  
  # Edit the year column to improve clarity
  data.temp$Year.t <- rep(c("2013","2014","2015","2016","2017","2018","2019"), length.out = length(data.temp$Year.t))
  
  # Add in a column defining the size individuals reached the following year (year t+1)
  data.temp$ID2 <- c(data.temp$ID[2:dim(data.temp)[1]],NA) # generate temporary ID column
  data.temp$Size.t1 <- c(data.temp$Size.t[2:dim(data.temp)[1]],NA)
  # use mismatched ID entries to remove data records inappropriately transferred from different individuals
  data.temp[which(data.temp$ID != data.temp$ID2),]$Size.t1 <- NA
  # remove temporary ID column
  data.temp$ID2 <- NULL
  
  # remove Year.t 2019 (there is no year.t in 2019 as individuals where not followed beyond this year - at least not in this dataset)
  data.temp <- data.temp[which(data.temp$Year.t != 2019),]
  
  # Add survival column (it is necessary to define survival in a separate column in order to return the size data to its numeric form)
  # NB: Survival != mortality
  data.temp$Surv <- NA
  data.temp[which(!is.na(data.temp$Size.t) & !is.na(data.temp$Size.t1)),]$Surv <- 1 # if there is data entered for both year t and year t+1 assign survival
  data.temp[which(data.temp$Size.t1 == "DEAD"),]$Surv <- 0 # However change this survival to mortality if the entry recorded is 'DEAD' 
  data.temp[which(data.temp$Size.t1 == "NFX"),]$Surv <- NA # and remove the entry if the individual was missed in year t+1 (NFX)
  data.temp[which(data.temp$Size.t == "NFX"),]$Surv <- NA # and remove the entry if the individual was missed in year t (NFX)
  
  # return size data to its numeric form
  data.temp$Size.t <- suppressWarnings(as.numeric(data.temp$Size.t)) # This will return a warning as NFX and DEAD entries are converted to NA. This is fine as these have already been store elsewhere
  data.temp$Size.t1 <- suppressWarnings(as.numeric(data.temp$Size.t1))
  
  # Finally, define whether the time interval is Pre-, During, or Post-hurricane?
  data.temp$Hurricane <- NA
  data.temp[which(data.temp$Year.t %in% c("2013","2014","2015","2016")),]$Hurricane <- "Pre" 
  data.temp[which(data.temp$Year.t %in% c("2017","2017.5")),]$Hurricane <- "Yes" 
  data.temp[which(data.temp$Year.t  == "2018"),]$Hurricane <- "Post" 
  
  # Reorder columns to tidy up!
  order <- c("Site","Species","ID","Year.t","Hurricane","Size.t","Size.t1","Surv")
  data.temp <- data.temp[order]
  # And store the newly transposed data
  R_data_IPM[[ii]] <- data.temp
}

# Stitch together the stack files for each of the different species/site combinations into one single data file
New_IPM_data <- ldply(R_data_IPM, data.frame)
# Ensure formatting across each variable
# Factor variables
New_IPM_data$Site <- as.factor(New_IPM_data$Site)
New_IPM_data$Species <- as.factor(New_IPM_data$Species)
New_IPM_data$ID <- as.factor(New_IPM_data$ID)
New_IPM_data$Year.t <- as.factor(New_IPM_data$Year.t)
New_IPM_data$Hurricane <- as.factor(New_IPM_data$Hurricane)

# and save the file
if(Outliers == FALSE){
  write.csv(New_IPM_data, "IPM_DATA", row.names = FALSE)
} else {
  write.csv(New_IPM_data, "IPM_DATA_FINAL", row.names = FALSE)
}

############### ------------------------------------------------------ END OF CODE -------------------------------------
