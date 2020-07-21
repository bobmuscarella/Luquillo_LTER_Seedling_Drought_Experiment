###########################################################################
###########################################################################
###########################################################################
### Combined code for Seedling Drought Experiment #########################
###########################################################################
###########################################################################
###########################################################################

#####################
### General setup ###
#####################

# Read data from Github
path <- "https://raw.github.com/bobmuscarella/Luquillo_LTER_Seedling_Drought_Experiment/master/Data/"

grow <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_growth.csv"))
surv <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_survival.csv"))
trait <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_traits.csv"))
photo <- read.csv(paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv"))

# Load packages
require(survival)
require(survminer)
require(coxme)
require(ggplot2)
require(RColorBrewer)
library(lme4)

# Set a color palette for species
cols <- brewer.pal(8, "Dark2")


#########################
### Survival analysis ###
#########################


### Run species-by-species survival models in a loop
surv_fits <- surv_fits_plot <- list()
options(na.action=na.exclude) # retain NA in predictions

for (i in seq_along(levels(surv$Species))){
  # Subset the data for the 'i-th' species
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]
  # Run the models with plot random effect (we tested other random effects but AIC selected this)
  surv_fits[[i]] <- coxme(Surv(Days, Status) ~ Moisture + Densiometer + May_LA
                             + (1|Plot), data=focdat)
  # Run models without random effect in order to be able to plot predicted response
  # (not currently possible with coxme models)
  surv_fits_plot[[i]] <- coxph(Surv(Days, Status) ~ Moisture + Densiometer + May_LA, 
                             data=focdat, model=T)
  # Name the items in the resulting list
  names(surv_fits)[i] <- levels(surv$Species)[i]
}

### Function to extract summary table from Cox ME models
extract_coxme_table <- function (mod){
  beta <- mod$coefficients
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2)
  p <- signif(1 - pchisq((beta/se)^2, 1), 2)
  table = data.frame(cbind(beta, se, z, p))
  return(table)
}

### Make a summary table for survival model results
t1 <- do.call(rbind, lapply(surv_fits, extract_coxme_table))
t1 <- round(t1, 3)
t1$Species <- unlist(lapply(strsplit(rownames(t1), "\\."), function(x)x[[1]]))
t1$Variable <- unlist(lapply(strsplit(rownames(t1), "\\."), function(x)x[[2]]))
t1 <- cbind(t1, round(do.call(rbind, lapply(surv_fits, confint)),3))
t1 <- t1[,c("Species","Variable","beta","2.5 %", "97.5 %","se","z","p")]

write.csv(t1, "Desktop/Table1.csv", row.names = F)






##############################################
##############################################
### Growth analysis ##########################
##############################################
##############################################

# The original growth data has plant size at each census.
# We first create the new data.frame to have growth during each interval.

growdf <- data.frame()
for(i in unique(grow$Order)){
  focdat <- grow[grow$Order==i,]
  species <- as.factor(focdat$Species[-1])
  growth <- (focdat$Leaf_area[-1] - focdat$Leaf_area[-9]) / (focdat$Days[-1] - focdat$Days[-9])
  startsize <- focdat$Leaf_area[-9]
  moisture <- focdat$Moisture[-1]
  densiometer <- focdat$Densiometer[-1]
  plot <- as.factor(focdat$Shelter[-1])
  indv <- as.factor(focdat$Order[-1])
  interval <- 1:8
  tmpdf <- data.frame(species, growth, moisture, densiometer, startsize, plot, indv, interval)
  growdf <- rbind(growdf, tmpdf)
  growdf <- growdf[!is.na(growdf$growth),]
}

# Fit the growth models
grow_fits <- list()

# grow_fits0 <- list()
# grow_fits1 <- list()
# grow_fits2 <- list()
# grow_fits3 <- list()
# grow_fits4 <- list()

for(sp in 1:length(levels(growdf$species))){
  
  print(levels(growdf$species)[sp])

  tmpdf <- growdf[growdf$species %in% levels(growdf$species)[sp],]
  
  # grow_fits0[[sp]] <- lmer(growth ~ moisture + densiometer + startsize + (1|plot) + (1|indv), 
  #                          data=tmpdf)
  # 
  # grow_fits1[[sp]] <- lmer(growth ~ moisture + I(moisture^2) + densiometer + startsize 
  #                          + (1|plot) + (1|indv), data=tmpdf)
  # 
  # grow_fits2[[sp]] <- lmer(growth ~ moisture + densiometer + startsize + (1|plot),
  #                         data=tmpdf)
  # 
  # grow_fits3[[sp]] <- lm(growth ~ moisture + densiometer + startsize, 
  #                        data=tmpdf)
  # 
  grow_fits[[sp]] <- lm(growth ~ moisture + densiometer + startsize, 
                         data=tmpdf)
  
  names(grow_fits)[sp] <- levels(growdf$species)[sp]
}

# moist_coeffs <- unlist(lapply(grow_fits, function(x) fixef(x)[2]))

# round(do.call(cbind, lapply(grow_fits, function(x) fixef(x))),3)

# Extract coefficients and compute confidence intervals
Species <- rep(names(grow_fits), each=4)
Variable <- rep(c("Intercept","Moisture","Densiometer","Start_size"), 8)
# Estimate <- round(do.call(c, lapply(grow_fits, fixef)),3)
Estimate <- round(do.call(c, lapply(grow_fits, coef)),3)
# cis <- do.call(rbind, lapply(grow_fits, function(x) round(confint(x)[-c(1:3),],3)))
cis <- do.call(rbind, lapply(grow_fits, function(x) round(confint(x),3)))

t2 <- as.data.frame(cbind(Species, Variable, Estimate, cis))
rownames(t2) <- NULL
t2



### Plot the point estimates with 95% CIs
# par(mar=c(6,5,2,2))
# plot(moist_coeffs, axes=F, xlab=NA, 
#      ylab='Estimated coefficient (growth x moisture)', 
#      pch=21, bg=cols, cex=2, ylim=range(cis))
# segments(1:8, cis[,1], 1:8, cis[,2], col=cols, lwd=2)
# abline(h=0, lty=2)
# axis(1, labels=levels(growdf$species), at=1:8, las=2)
# axis(2)
# box()


###############################################
### Plotting Effects on Survival and Growth ###
###############################################

pdf("/Users/au529793/Desktop/Figure1.pdf", width = 6, height = 8)

# Predict survival as a function of soil moisture by species
pred_moist <- list()
pred_light <- list()
pred_size <- list()

for (i in seq_along(levels(surv$Species))){
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]
  newdata_moist <- data.frame(Moisture = 0:50, 
                              Densiometer=mean(focdat$Densiometer),
                              May_LA=mean(focdat$May_LA, na.rm=T), Days=228, Status=1)
  pred_moist[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_moist, type='expected'))
  newdata_light <- data.frame(Moisture=mean(focdat$Moisture),
                              Densiometer=0:25,
                              May_LA=mean(focdat$May_LA, na.rm=T), Days=228, Status=1)
  pred_light[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_light, type='expected'))
  newdata_size <- data.frame(Moisture=mean(focdat$Moisture),
                             Densiometer=mean(focdat$Densiometer),
                             May_LA=seq(min(surv$May_LA, na.rm=T), 
                                        max(surv$May_LA, na.rm=T), 
                                        length.out = 51), Days=228, Status=1)
  pred_size[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_size, type='expected'))
}

par(mfcol=c(3,2), mar=c(4,5,1,0.5), oma=c(1,1,3,1))

plot(0:50, ylim=c(0,100), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab="Predicted Survival (%)")
mtext("Survival", 3, 1)
for(i in 1:length(pred_moist)){
  lty <- ifelse(t1$p[t1$Variable=="Moisture"][i]<0.05, 1, 2)
  lwd <- ifelse(t1$p[t1$Variable=="Moisture"][i]<0.05, 3, 2)
  points(0:50, 100*pred_moist[[i]], col=cols[i], type='l', lwd=lwd, lty=lty)
}

plot(0:25, ylim=c(0,100), pch=NA,
     xlim=c(-0.1,25), xlab="Canopy openness (%)", 
     ylab="Predicted Survival (%)")
for(i in 1:length(pred_light)){
  lty <- ifelse(t1$p[t1$Variable=="Densiometer"][i]<0.05, 1, 2)
  lwd <- ifelse(t1$p[t1$Variable=="Densiometer"][i]<0.05, 3, 2)
  points(0:25, 100*pred_light[[i]], col=cols[i], type='l', lwd=lwd, lty=lty)
}

plot(1,1, xlim=range(surv$May_LA, na.rm=T), 
     ylim=c(0,100), pch=NA,
     xlab=bquote("Initial leaf area (cm"^2*")"), 
     ylab="Predicted Survival (%)")
for(i in 1:length(pred_size)){
  x <- seq(min(surv$May_LA, na.rm=T), 
           max(surv$May_LA, na.rm=T), 
           length.out = 51)
  lty <- ifelse(t1$p[t1$Variable=="May_LA"][i]<0.05, 1, 2)
  lwd <- ifelse(t1$p[t1$Variable=="May_LA"][i]<0.05, 3, 2)
  points(x, 100*pred_size[[i]], col=cols[i], type='l', lwd=lwd, lty=lty)
}
legend('bottomright', legend=levels(surv$Species),
       cex=0.7, bty='n', lty=1, col=cols, lwd=2)

plot(0:50, ylim=c(-0.5,1), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab=bquote("Pred. Growth (cm"^2~"day"^-1*")"))
for(i in seq_along(levels(growdf$species))){
  tmpdf <- growdf[growdf$species == levels(growdf$species)[i],]
  newdata <- data.frame(moisture=0:50,
                        densiometer=mean(tmpdf$densiometer),
                        startsize=mean(tmpdf$startsize))
  lty <- ifelse(summary(grow_fits[[i]])$coeff["moisture", 4] < 0.05, 1, 2)
  lwd <- ifelse(summary(grow_fits[[i]])$coeff["moisture", 4] < 0.05, 3, 2)
  pred <- predict(grow_fits[[i]], newdata=newdata, re.form=NA)
  points(0:50, pred, col=cols[i], type='l', lwd=lwd, lty=lty)
}
legend('topleft', legend=levels(surv$Species),
       cex=0.7, bty='n', lty=1, col=cols, lwd=2)
mtext("Growth", 3, 1)


plot(0:25, ylim=c(-0.25,1), pch=NA,
     xlab="Canopy openness (%)", 
     ylab=bquote("Pred. Growth (cm"^2~"day"^-1*")"))
for(i in seq_along(levels(growdf$species))){
  tmpdf <- growdf[growdf$species == levels(growdf$species)[i],]
  newdata <- data.frame(moisture=mean(tmpdf$moisture),
                        densiometer=0:25,
                        startsize=mean(tmpdf$startsize))
  lty <- ifelse(summary(grow_fits[[i]])$coeff["densiometer", 4] < 0.05, 1, 2)
  lwd <- ifelse(summary(grow_fits[[i]])$coeff["densiometer", 4] < 0.05, 3, 2)
  pred <- predict(grow_fits[[i]], newdata=newdata, re.form=NA)
  points(0:25, pred, col=cols[i], type='l', lwd=lwd, lty=lty)
}

plot(0:50, ylim=c(-1,1.5), pch=NA,
     xlab=bquote("Initial leaf area (cm"^2*")"), 
     ylab=bquote("Pred. Growth (cm"^2~"day"^-1*")"))
for(i in seq_along(levels(growdf$species))){
  tmpdf <- growdf[growdf$species == levels(growdf$species)[i],]
  newdata <- data.frame(moisture=mean(tmpdf$moisture),
                        densiometer=mean(tmpdf$densiometer),
                        startsize=0:50)
  lty <- ifelse(summary(grow_fits[[i]])$coeff["startsize", 4] < 0.05, 1, 2)
  lwd <- ifelse(summary(grow_fits[[i]])$coeff["startsize", 4] < 0.05, 3, 2)
  pred <- predict(grow_fits[[i]], newdata=newdata, re.form=NA)
  points(0:50, pred, col=cols[i], type='l', lwd=lwd, lty=lty)
}

dev.off()




