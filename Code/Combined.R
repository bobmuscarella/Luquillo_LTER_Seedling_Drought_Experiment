### Combined code for Seedling Drought Experiment

# Read data from Github
path <- "https://raw.github.com/bobmuscarella/Luquillo_LTER_Seedling_Drought_Experiment/master/Data/"

growth_url <- paste0(path, "LUQ_DroughtExp_Seedling_growth.csv")
survive_url <- paste0(path, "LUQ_DroughtExp_Seedling_survival.csv")
traits_url <- paste0(path, "LUQ_DroughtExp_Seedling_traits.csv")
photo_url <- paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv")

grow <- read.csv(growth_url)
surv <- read.csv(survive_url)
trait <- read.csv(traits_url)
photo <- read.csv(photo_url)


#########################
### Survival analysis ###
#########################

library(survival)
library(survminer)
library(coxme)
library(ggplot2)

### Run species-by-species survival models in a loop
surv_results <- list()
for (i in seq_along(levels(surv$Species))){
  # Subset the data for the 'i-th' species
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]
  # Run the models (we had experimented with other random effects but AIC selected this)
  surv_results[[i]] <- coxme(Surv(Days, Status) ~ Moisture + Densiometer + May_LA 
                             + (1|Plot), data=focdat)
  names(surv_results)[i] <- levels(surv$Species)[i]
}

# Build a summary table for survival model results
surv_results




################## Graphing ###########################
PREMON_test.cox <- predict(surv_results$PREMON)

# The following causes the error (x and y lengths differ)
plot(surv$Moisture[surv$Species %in% "PREMON"], PREMON_test.cox)

# This is because observations with NA for any co-variate are excluded when fitting the model
# Some rows have NA for the initial size (May_LA) variable
# and that gives you the 'unequal' lengths error

# You can remove those records and then look at it
PREMON_moist <- surv$Moisture[surv$Species %in% "PREMON" & !is.na(surv$May_LA)]


UREBAC_moist <- seedlingdata$Moisture[seedlingdata$species %in% "UREBAC" & 
                                        !is.na(seedlingdata$May_LA)]

MANBID_moist <- seedlingdata$Moisture[seedlingdata$species %in% "MANBID" & 
                                        !is.na(seedlingdata$May_LA)]

TETBAL_moist <- seedlingdata$Moisture[seedlingdata$species %in% "TETBAL" & 
                                        !is.na(seedlingdata$May_LA)]
INLAU_moist <- seedlingdata$Moisture[seedlingdata$species %in% "INLAU" & 
                                       !is.na(seedlingdata$May_LA)]
GUAGUI_moist <- seedlingdata$Moisture[seedlingdata$species %in% "GUAGUI" & 
                                        !is.na(seedlingdata$May_LA)]
##BOB QUESTION: INSTEAD OF THE LINEAR PREDICTOR I THINK IT MAKES MORE SENSE TO 
## GET THE SURVIVAL PROBABILITY, WHICH ACCORDING TO THE HELP PAGE WE GET BY 
## type='expected' AND THEN ITS exp(-expected), WHICH I TRIED OUT BELOW AND IT WORKED
## THE Y-AXIS UNITS MAKE SENSE (SURVIVAL PROB SHOULD BE 0-1) AND THE OVERALL SHAPE MAKES SENSE

pred <- predict(results$PREMON$RE_none, type='expected')
Pred<- exp(-pred)
plot(moist, Pred)
#plot(moist, pred)

##BOB QUESTION: BUT WHEN I TRY TO DO THE SAME THING WHEN PREDICTING TO 
## NEW DATA IT PRODUCES A GRAPH THAT DOESN'T MAKE SENSE (Y AXIS DOESN'T MAKE SENSE AND SHAPE SEEMS WEIRD),
## NOT SURE IF THE type='expected' IS IN THE CORRECT SPOT OR SOMETHING ELSE? 
## MAYBE THIS HAS TO DO WITH THE NEW DATA NOT BEING 'SUPPORTED' WITH THE PREDICT FUNCTION?
## NOT SURE HOW ELSE WE CAN SHOW THE FITTED LINE
# You can also predict to new data and that way see the fitted line

#PREMON
premon_try <- (predict(results$PREMON$RE_none, 
                       newdata=data.frame(Moisture=seq(0,50, length.out = 51), type='expected',
                                          densiometer=rep(
                                            mean(seedlingdata$densiometer
                                                 [seedlingdata$species=="PREMON"]), 51),
                                          May_LA=rep(
                                            mean(seedlingdata$May_LA
                                                 [seedlingdata$species=="PREMON"], na.rm=T), 51))))

##TAKE THE NEG EXPONENT
Premon_predict<- exp(-premon_predict)
Urebac_predict<- exp(-urebac_predict)


##THE POINTS SEEM CORRECT BUT THE LINE ISN'T CORRECT
plot(seq(0,50, length.out = 51), Test, type='l', ylim=c(0,1))
test2 <- exp(-predict(results$PREMON$RE_none, type='expected'))
points(moist, test2, col=4)

plot(seq(0,50, length.out = 51), premon_predict, ylim=c(0,1))
test2 <- exp(-predict(results$PREMON$RE_none, type='expected'))
points(moist, test2, col=4)


Premon_predict <- exp(-predict(results$PREMON$RE_none, type='expected'))
Urebac_predict <- exp(-predict(results$UREBAC$RE_none, type='expected'))
Manbid_predict<- exp(-predict(results$MANBID$RE_none, type='expected'))
Tetbal_predict <- exp(-predict(results$TETBAL$RE_none, type='expected'))
Inlau_predict<- exp(-predict(results$INLAU$RE_none, type='expected'))
Guagui_predict<- exp(-predict(results$GUAGUI$RE_none, type='expected'))


plot(PREMON_moist, Premon_predict, col=1, pch = 19, ylab="Survival probability", xlab="Soil moisture",ylim=c(0,1) )
points(UREBAC_moist, Urebac_predict, col=2, pch = 19 )
points(MANBID_moist, Manbid_predict, col=3, pch = 19 )
points(TETBAL_moist, Tetbal_predict, col=4, pch = 19 )
points(INLAU_moist, Inlau_predict, col=5, pch = 19 )
points(GUAGUI_moist, Guagui_predict, col=6, pch = 19 )
legend(35, 0.4, legend=c("PREMON", "UREBAC", "MANBID", "TETBAL", "INLAU", "GUAGUI"),
       col=c("1", "2", "3", "4", "5", "6"), pch = 19)


points(moist, Premon_predict, col=1)
points(moist, Manbid_predict, col=3)
points(moist, Tetbal_predict, col=4)
points(moist, Inlau_predict, col=5)





# Note there that the default 'type' of prediction for the cox model is "lp", which
# apparently means 'linear predictor'.  There are other options; see these help files:
?predict.coxph # This allows for 'newdata' and different 'types'
?predict.coxme # This says 'newdata' is not yet supported but it allows for different 'types' of predictions







##############################################
##############################################
### Growth analysis ##########################
##############################################
##############################################


# Here is the growth script and data, which hasn’t been changed since you so expertly crafted it. The script needs to be modified so the coefficients (and 95% CI) for soil moisture, densitometer & starting size are saved for each species from the best model to be used in plotting with the PCA scores.  


# setwd("C:\\Users\\Matlaga\\Documents\\R")

# For Bob
# setwd("/Users/au529793/Downloads/fw")

library(lme4)

d <- read.table("Seedling growth june 12.txt",header=TRUE,sep="\t")

# Create the new data.frame to have growth during each interval

df <- data.frame()
for(i in unique(d$Order)){
  focdat <- d[d$Order==i,]
  species <- as.factor(focdat$Species[-1])
  growth <- (focdat$Leaf_area[-1] - focdat$Leaf_area[-9]) / (focdat$days[-1] - focdat$days[-9])
  startsize <- focdat$Leaf_area[-9]
  # Note here we are using the grand mean soil moisture...
  moisture <- focdat$Moisture[-1]
  densiometer <- focdat$densiometer[-1]
  plot <- as.factor(focdat$Shelter[-1])
  indv <- as.factor(focdat$Order[-1])
  interval <- 1:8
  tmpdf <- data.frame(species, growth, moisture, densiometer, startsize, plot, indv, interval)
  df <- rbind(df, tmpdf)
  df <- df[!is.na(df$growth),]
}


cols <- RColorBrewer::brewer.pal(8, "Dark2")
plot(df$moisture, df$growth, col=scales::alpha(cols[df$species],0.5), pch=16)
abline(h=0, lty=2)

fits <- list()
for(sp in 1:length(levels(df$species))){
  print(levels(df$species)[sp])
  tmpdf <- df[df$species %in% levels(df$species)[sp],]
  head(tmpdf)
  fits[[sp]] <- lmer(growth ~ moisture + densiometer + startsize + (1|plot) + (1|indv),
                     data=tmpdf)
  names(fits)[sp] <- levels(df$species)[sp]
}

# Extract the coefficient for soil moisture
moist_coeffs <- unlist(lapply(fits, function(x) fixef(x)[2]))

# Get 95% confidence intervals on the soil moisture effect
moist_coeffs_cis <- do.call(rbind, lapply(fits, function(x) confint(x)['moisture',]))

# Plot the point estimates with 95% CIs
par(mar=c(6,5,2,2))
plot(moist_coeffs, axes=F, xlab=NA, 
     ylab='Estimated coefficient (growth x moisture)', 
     pch=21, bg=cols, cex=2, ylim=range(cis))
segments(1:8, cis[,1], 1:8, cis[,2], col=cols, lwd=2)
abline(h=0, lty=2)
axis(1, labels=levels(df$species), at=1:8, las=2)
axis(2)
box()







