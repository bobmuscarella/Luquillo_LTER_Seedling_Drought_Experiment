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



### Survival analysis
library(survival)
library(survminer)
library(coxme)
library(ggplot2)

# the `seedlingdata` object has some weird extra rows at the end... I remove those here:
seedlingdata <- survival[1:960,]

# Make sure the 'species' column is a factor and drop the blank factor level that originated from the extra rows above
seedlingdata$species <- droplevels(as.factor(seedlingdata$species))

### Run all the models in a loop

# Make a list to hold the results
results <- list()

# Run the loop, all the models for one species at a time
for (i in seq_along(levels(seedlingdata$species))){
  
  # Subset the data for the 'i-th' species
  focdat <- seedlingdata[seedlingdata$species %in% levels(seedlingdata$species)[i],]
  
  # Run the models (including a model with only plot as rand eff)
  fit1 <- coxph(Surv(days, status) ~ Moisture + densiometer + May_LA, data=focdat, model=T)
  fit2 <- coxme(Surv(days, status) ~ Moisture + densiometer + May_LA 
                + (1|Plot), data=focdat)
  fit3 <- coxme(Surv(days, status) ~ Moisture + densiometer + May_LA 
                + (1|Plot/ID), data=focdat)
  
  # Make the results into a smaller list of three that goes into one space of the 'results' list
  results[[i]] <- list(RE_none=fit1, RE_plot=fit2, RE_plotID=fit3)
  
  # Name each element of the results list as the species
  names(results)[i] <- levels(seedlingdata$species)[i]
}

# ANOVA for each species, all at the same time
lapply(results, function(x) anova(x[[1]], x[[2]], x[[3]]))

# Summarize all the models
lapply(results, function(x) lapply(x, summary))


################## Graphing ###########################
PREMON_test.cox <- predict(results$PREMON$RE_plot)

# The following causes the error (x and y lengths differ)
plot(seedlingdata$Moisture[seedlingdata$species %in% "PREMON"], PREMON_test.cox)

# This is because observations with NA for any co-variate are excluded when fitting the model
# Some rows have NA for the initial size (May_LA) variable
# and that gives you the 'unequal' lengths error

# You can remove those records and then look at it
PREMON_moist <- seedlingdata$Moisture[seedlingdata$species %in% "PREMON" & 
                                        !is.na(seedlingdata$May_LA)]

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







