###########################################################################
###########################################################################
### Combined code for Seedling Drought Experiment #########################
###########################################################################
###########################################################################

#####################
### General setup ###
#####################

### LOAD PACKAGES
require(survival)
require(survminer)
require(coxme)
require(lme4)
library(MuMIn)
require(factoextra)
require(vegan)
require(plotfunctions)
# require(ggplot2)
require(RColorBrewer)
require(scales)
require(dplyr)


species_labs <- c("Cecropia schreberiana",
            "Guarea guidonia",
            "Inga laurina",
            "Manilkara bidentata",
            "Prestoea montana",
            "Schefflera morototoni",
            "Tetragastris balsamifera",
            "Urera baccifera")
sp_labs <- c("C. schreberiana",
            "G. guidonia",
            "I. laurina",
            "M. bidentata",
            "P. montana",
            "S. morototoni",
            "T. balsamifera",
            "U. baccifera")


### READ DATA FROM GITHUB
path <- "https://raw.github.com/bobmuscarella/Luquillo_LTER_Seedling_Drought_Experiment/master/Data/"

### Full soil moisture measurements
moist <- read.csv(paste0(path, "LUQ_DroughtExp_Soil_moisture_complete.csv"))

# Growth data
grow <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_growth.csv"))
grow$Date <- as.Date(as.character(grow$Date), format="%m/%d/%y")
grow$ID <- paste(grow$Plot, grow$Position, sep = '.')
grow$Species <- as.factor(grow$Species)

### Survival data
surv <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_survival.csv"))
surv$ID <- paste(surv$Plot, surv$Position, sep = '.')
surv$Species <- as.factor(surv$Species)

### Trait data
trait <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_traits.csv"))
trait$ID <- paste(trait$Plot, trait$Position, sep = '.')
trait$Species <- as.factor(trait$Species)

### Photosynthesis data
# photo <- read.csv(paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv"))
# photo$ID <- paste(photo$Plot, photo$Position, sep = '.')

### Nutrient data
# nutrient <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_leafnutrients.csv"))
# # Add empty factor levels to keep species order consistent with other datasets 
# # (Because not all species were included in the nutrient analysis)
# nutrient$Species <- factor(nutrient$Species, levels = levels(grow$Species))
# nutrient$ID <- paste(nutrient$Plot, nutrient$Position, sep = '.')

# Pre-dawn water potential data
# pdwp


### SET A COLOR PALETTE TO IDENTIFY SPECIES
cols <- brewer.pal(8, "Dark2")


### LOAD CUSTOM HELPER FUNCTION (to extract summary table from Cox ME models)
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


####################################################
### Figure of soil moisture through time by plot ###
####################################################

### Get average moisture (from 4 measurements) per plot per time period
moist_summary <- moist %>% group_by(Plot) %>% summarize_each(funs=mean)
moist_summary$Treatment <- moist$Treatment[match(moist_summary$Plot, moist$Plot)]

moist_summary_min <- moist %>% group_by(Plot) %>% summarize_each(funs=min)
moist_summary_min$Treatment <- moist$Treatment[match(moist_summary_min$Plot, moist$Plot)]

### Remove one survey (July 29, 2019) where plots had been mixed up in the data
moist_summary <- moist_summary[,-7]
moist_summary_min <- moist_summary_min[,-7]

plot(apply(moist_summary[,-c(1:2)], 1, mean, na.rm=T), 
     apply(moist_summary_min[,-c(1:2)], 1, min, na.rm=T))
cor.test(apply(moist_summary[,-c(1:2)], 1, mean, na.rm=T), 
         apply(moist_summary_min[,-c(1:2)], 1, min, na.rm=T))

labs <- as.Date(c("2019-05-31",
                  "2019-06-18",
                  "2019-07-1",
                  "2019-07-15",
                  # "2019-07-29",
                  "2019-08-15",
                  "2019-08-26",
                  "2019-09-09",
                  "2019-09-23",
                  "2019-10-23",
                  "2019-11-22",
                  "2019-12-05",
                  "2020-01-04"))

# pdf("Figures/FigureS0_Soil_Moisture_through_time.pdf")

plot(labs, moist_summary[1,3:ncol(moist_summary)], 
     pch=NA, ylim=c(0,55), las=2, 
     xlab="Month", ylab='Soil Moisture %')

for(i in 1:nrow(moist_summary)){
  lines(labs, moist_summary[i,3:ncol(moist_summary)], 
        col=ifelse(moist_summary$Treatment[i]=='D', 2, 4))
}

legend('topleft', legend = c("Control","Drought"), 
       lty=1, col=c(4,2), bty='n', lwd=3)

# dev.off()


min_moist <- apply(moist_summary_min[,-c(1:2)], 1, min, na.rm=T)
names(min_moist) <- moist_summary_min$Plot
grow$min_moist <- min_moist[match(grow$Plot, names(min_moist))]

# Get the 10% and 90%-iles of soil moisture from the data to use in defining sensitivity
moist_percentiles <- quantile(moist_summary[,-c(1:2)], na.rm=T, probs=c(0.1, 0.9))

######################################################################################
### Figure of individual growth trajectories over time with soil moisture and mortality ###
######################################################################################

# pdf("Figures/Figure1.pdf", height=6, width=8)

# Color based on soil moisture in each plot
grow$soil_col <- brewer.pal(10, "Spectral")[cut(grow$Moisture, 10)]

# Reorder the data for plotting convenience
grow <- grow[order(grow$Species, grow$Order, grow$Census),]

par(mfrow=c(2,4), mar=c(4,4,2,1))
for (sp in 1:length(unique(grow$Species))){
  focsp <- levels(grow$Species)[sp]
  tmpdat <- grow[grow$Species %in% focsp,]
  
  # Make a blank plot to hold the lines drawn below
  plot(tmpdat$Date, tmpdat$Leaf_area, pch=NA, 
       xlab=NA, ylab="Leaf Area (cm2)",
       main=sp_labs[sp], las=2, font.main=3)
  mtext(LETTERS[sp], 3, -1.5, at=as.Date("2019-6-10"))
  
  for(i in 1:length(unique(tmpdat$Order))){
    focdat <- tmpdat[tmpdat$Order == unique(tmpdat$Order)[i],]
    focdat <- focdat[!is.na(focdat$Leaf_area),]
    lines(focdat$Date, focdat$Leaf_area, col=focdat$soil_col, lwd=0.5)
  }
  
  for(i in 1:length(unique(tmpdat$Order))){
    focdat <- tmpdat[tmpdat$Order == unique(tmpdat$Order)[i],]
    focdat <- focdat[!is.na(focdat$Leaf_area),]
    if(nrow(focdat)<9){
      points(focdat$Date[nrow(focdat)], focdat$Leaf_area[nrow(focdat)], 
             pch=4, cex=0.35, lwd=0.75)
    }
  }
  
  # Add a legend to the panel with CECSCH (it fits best here)
  if (focsp=="CECSCH"){
    plotfunctions::gradientLegend(
      round(range(grow$Moisture),0),
      n.seg=1,
      color = brewer.pal(10, "Spectral"),
      
      pos=c(18200, 22, 
            18220, 32), coords=TRUE,
      
      # pos=c(as.Date("2019-11-1"), 22, 
      #       as.Date("2019-11-20"), 32), coords=TRUE,
      dec=0, border.col = NA, fit.margin = F
    )
    text(as.Date("2019-11-15"), 33.5, "Soil", cex=0.8, font=1, pos=3)
    text(as.Date("2019-11-15"), 32, "Moisture %", cex=0.8, font=1, pos=3)
  }
}

# dev.off()


################################
### Plot raw survival curves ###
################################

# pdf("Figures/FigureS6.pdf", height=6, width=8)

par(mfrow=c(1,1), mar=c(5,5,2,1))

surv_curves_fit <- survfit(Surv(Days, Status) ~ Species, data = surv)  
plot(surv_curves_fit, col=cols, lwd=2, axes=F, xlim=c(0,250), 
     xlab="Days", ylab="Percent Surviving")
legend('bottomleft', legend=sp_labs,
       cex=0.7, bty='n', lty=1, col=cols, lwd=2, text.font=3)
axis(2, at=seq(0,1,by=0.2), labels=seq(0,100,by=20))
axis(1)

# dev.off()


#########################
### Survival analysis ###
#########################

### Run species-by-species survival models in a loop
surv_fits <- surv_fits_plot <- surv_fits_int <-list()
options(na.action=na.exclude) # retain NA in predictions

for (i in seq_along(levels(surv$Species))){
  # Subset the data for the ith species
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]

  # Run the models with plot random effect (we tested other random effects but AIC selected this)
  surv_fits[[i]] <- coxme(Surv(Days, Status) ~ Moisture + Densiometer + Start_LA
                             + (1|Plot), data=focdat)
  surv_fits_int[[i]] <- coxme(Surv(Days, Status) ~ Moisture + Densiometer + Moisture*Densiometer
                          + Start_LA + (1|Plot), data=focdat)
    
  # Run models without random effect in order to be able to plot predicted response
  # (not currently possible with coxme models)
  surv_fits_plot[[i]] <- coxph(Surv(Days, Status) ~ Moisture + Densiometer #+ Moisture*Densiometer
                               + Start_LA, data=focdat, model=T)
  
  # Name the items in the resulting list
  names(surv_fits)[i] <- levels(surv$Species)[i]
  names(surv_fits_int)[i] <- levels(surv$Species)[i]
  names(surv_fits_plot)[i] <- levels(surv$Species)[i]
}


# Check AIC to see if interaction is useful:  
x <- data.frame(cbind(do.call(rbind, lapply(surv_fits, AIC)),
                      do.call(rbind, lapply(surv_fits_int, AIC))))
names(x) <- c("No interaction", "With interation")
x$deltaAIC <- x[,1] - x[,2]

for(i in 1:8){
  print(anova(surv_fits[[i]], surv_fits_int[[i]]))
}


### Make a summary table for survival model results
t1 <- do.call(rbind, lapply(surv_fits, extract_coxme_table))
t1 <- round(t1, 3)
t1$Species <- unlist(lapply(strsplit(rownames(t1), "\\."), function(x)x[[1]]))
t1$Variable <- unlist(lapply(strsplit(rownames(t1), "\\."), function(x)x[[2]]))
t1 <- cbind(t1, round(do.call(rbind, lapply(surv_fits, confint)),3))
t1 <- t1[,c("Species","Variable","beta","2.5 %", "97.5 %","se","z","p")]

# write.csv(t1, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/Table1_survival_summary.csv", row.names = F)

##############################################
### Growth analysis ##########################
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
  id <- focdat$ID[-1]
  plot <- as.factor(focdat$Plot[-1])
  indv <- as.factor(focdat$Order[-1])
  interval <- 1:8
  tmpdf <- data.frame(id, species, growth, moisture, densiometer, startsize, plot, indv, interval)
  growdf <- rbind(growdf, tmpdf)
  growdf <- growdf[!is.na(growdf$growth),]
}


growdf$moisture.z <- scale(growdf$moisture)

### Fit the growth models
grow_fits_int <- grow_fits <- grow_fits2 <- list()

for(sp in 1:length(levels(growdf$species))){
  print(levels(growdf$species)[sp])
  tmpdf <- growdf[growdf$species %in% levels(growdf$species)[sp],]

  # With plot and individual random effects
  grow_fits[[sp]] <- lmer(growth ~ moisture + densiometer + startsize 
                           + (1|plot) + (1|indv), data=tmpdf)
  
  # With plot random effect only
  grow_fits2[[sp]] <- lmer(growth ~ moisture + densiometer + startsize 
                           + (1|plot), data=tmpdf)
  
  # With interaction term
  grow_fits_int[[sp]] <- lmer(growth ~ moisture + densiometer 
                              + moisture*densiometer + startsize 
                              + (1|plot) + (1|indv), data=tmpdf)
  
  
  names(grow_fits_int)[sp] <- levels(growdf$species)[sp]
  names(grow_fits)[sp] <- levels(growdf$species)[sp]
  names(grow_fits2)[sp] <- levels(growdf$species)[sp]
}

aic_comp <- data.frame(no_int = unlist(lapply(grow_fits, AIC)),
                       no_int2 = unlist(lapply(grow_fits2, AIC)),
                       with_int = unlist(lapply(grow_fits_int, AIC)))

# No support for including an interaction term
aic_comp$no_int_preferred <- aic_comp$no_int < aic_comp$with_int
aic_comp$no_indv_randeff_preferred <- aic_comp$no_int2 < aic_comp$no_int

aic_comp

# Extract coefficients and compute confidence intervals
Species <- rep(names(grow_fits), each=4)
Variable <- rep(c("Intercept","Moisture","Densiometer","Start_size"), 8)
# Species <- rep(names(grow_fits), each=5)
# Variable <- rep(c("Intercept","Moisture","Densiometer","SM*Light","Start_size"), 8)

Estimate <- round(do.call(c, lapply(grow_fits, fixef)),3)

growcis <- do.call(rbind, lapply(grow_fits, function(x) round(confint(x)[-(1:3),],3)))

t2 <- as.data.frame(cbind(Species, Variable, Estimate, growcis))
t2$`2.5 %` <- as.numeric(as.character(t2$`2.5 %`))
t2$`97.5 %` <- as.numeric(as.character(t2$`97.5 %`))
t2$sig <- sign(t2$`2.5 %`)==sign(t2$`97.5 %`)
rownames(t2) <- NULL

r2 <- round(do.call(rbind, lapply(grow_fits, r.squaredGLMM)), 2)
r2

# write.csv(t2, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/Table2_growth_summary_20211103.csv", row.names = F)

########################################################################
### VISUALIZE THE INTERACTION TERMS FOR GROWTH AND SURVIVAL MODELS

# pdf("/Users/au529793/Desktop/growth_interactions.pdf", height=8, width=6)
# 
# par(mfrow=c(4,2), mar=c(4,5,1,1))
# 
# for(i in 1:length(grow_fits)){
# tmp <- grow_fits[[i]]
# 
# ndlolight <- data.frame(moisture=seq(0,100,length.out=100), 
#                       densiometer=rep(3, 100),
#                       startsize=rep(50,100))
# 
# ndhilight <- data.frame(moisture=seq(0,100,length.out=100), 
#                       densiometer=rep(15, 100),
#                       startsize=rep(50,100),
#                       plot=rep(1,100))
# 
# lolight_pred <- predict(tmp, newdata=ndlolight,  re.form=NA)
# hilight_pred <- predict(tmp, newdata=ndhilight,  re.form=NA)
# 
# plot(seq(0,100,length.out=100), lolight_pred, type='l', 
#      ylim=range(c(hilight_pred, lolight_pred)), lwd=3,
#      xlab='Soil moisture (%)', ylab='Pred. growth', main=names(grow_fits)[i])
# lines(seq(0,100,length.out=100), hilight_pred, col=2, lwd=3)
# legend('topleft', legend=c('Low light (3%)', 'High light (15%)'), 
#        lty=1, col=1:2, bty='n', lwd=3)
# }
# 
# dev.off()
# 
# 
# ndlosm <- data.frame(moisture=rep(8, 100), 
#                         densiometer=seq(3,15,length.out=100),
#                         startsize=rep(50,100))
# 
# ndhism <- data.frame(moisture=rep(40, 100), 
#                         densiometer=seq(3,15,length.out=100),
#                         startsize=rep(50,100),
#                         plot=rep(1,100))
# 
# ndlosm_pred <- predict(grow_fits$INGLAU, newdata=ndlosm, re.form=NA)
# ndhism_pred <- predict(grow_fits$INGLAU, newdata=ndhism, re.form=NA)
# 
# plot(seq(3,15,length.out=100), ndlosm_pred, type='l', 
#      ylim=range(c(ndlosm_pred, ndhism_pred)), main='INGLAU', lwd=3,
#      xlab='Canopy openness (%)', ylab='Pred. growth')
# lines(seq(3,15,length.out=100), ndhism_pred, col=2, lwd=3)
# legend('topleft', legend=c('Low moisture (8%)', 'High moisture (40%)'), 
#        lty=1, col=1:2, bty='n', lwd=3)


# Predict survival as a function of soil moisture by species
pred_moist <- list()
pred_moist_lolight <- list()
pred_moist_hilight <- list()

for (i in seq_along(levels(surv$Species))){
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]
  newdata_moist <- data.frame(Moisture = 0:50,
                              Densiometer=mean(focdat$Densiometer),
                              Start_LA=mean(focdat$Start_LA, na.rm=T), Days=228, Status=1)
  pred_moist[[i]] <- predict(surv_fits_plot[[i]], newdata=newdata_moist, type='survival')
  
  newdata_moist_lolight <- data.frame(Moisture = 0:50,
                                      Densiometer=3,
                                      Start_LA=mean(focdat$Start_LA, na.rm=T), Days=228, Status=1)
  pred_moist_lolight[[i]] <- predict(surv_fits_plot[[i]], newdata=newdata_moist_lolight, 
                                     type='survival')
  
  newdata_moist_hilight <- data.frame(Moisture = 0:50,
                                      Densiometer=15,
                                      Start_LA=mean(focdat$Start_LA, na.rm=T), Days=228, Status=1)
  pred_moist_hilight[[i]] <- predict(surv_fits_plot[[i]], newdata=newdata_moist_hilight, 
                                     type='survival')
  
  names(pred_moist_hilight)[i] <- as.character(unique(focdat$Species))
  
}

#### SURVIVAL INTERACTIONS
# # pdf("/Users/au529793/Desktop/survival_interactions.pdf", height=8, width=6)
# 
# par(mfrow=c(4,2), mar=c(4,5,1,1))
# 
# for(i in 1:length(pred_moist_hilight)){
#   plot(0:50, ylim=c(0,100), pch=NA,
#        xlab="Soil Moisture (%)", 
#        ylab="Pred. Survival (%)", 
#        main=names(pred_moist_hilight)[i])
#   lines(0:50, 100*pred_moist_hilight[[i]], col=2, lwd=3)
#   lines(0:50, 100*pred_moist_lolight[[i]], col=1, lwd=3)
#   legend('topleft', legend=c('Low light (3%)', 'High light (15%)'), 
#          lty=1, col=1:2, bty='n', lwd=3)
# }
# 
# dev.off()




###########################################
### Figure Showing Effects on Survival and Growth ###
###########################################



# pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure2_new.pdf", width=7, height=4)

par(mfcol=c(1,2), mar=c(5,5,3,0.5))

plot(0:50, ylim=c(0,100), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab="Pred. Survival (%)")
mtext("A", 3, -1.5, at=3, cex=1.5)
for(i in 1:length(pred_moist)){
  lty <- ifelse(sum(t1[t1$Variable=='Moisture', c("2.5 %","97.5 %")][i,]>=0)!=1, 1, 2)
  points(0:50, 100*pred_moist[[i]], col=cols[i], type='l', lwd=3, lty=lty)
}

abline(v=moist_percentiles, lty=1, col='grey')

plot(0:50, ylim=c(-0.75,1.5), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab=bquote("Pred. Leaf Area Growth (cm"^2~"day"^-1*")"))
mtext("B", 3, -1.5, at=50, cex=1.5)
polygon(c(-10,70,70,-10), c(0,0,-2,-2), col='lightgrey', lty=0)
box()
for(i in seq_along(levels(growdf$species))){
  focsp <- levels(growdf$species)[i]
  tmpdf <- growdf[growdf$species == focsp,]
  newdata <- data.frame(moisture=0:50,
                        densiometer=mean(tmpdf$densiometer),
                        startsize=mean(tmpdf$startsize))
  lty <- ifelse(sum(t2[t2$Species==focsp & t2$Variable=='Moisture', c("2.5 %","97.5 %")]>=0)!=1, 1, 2)
  lwd <- ifelse(sum(t2[t2$Species==focsp & t2$Variable=='Moisture', c("2.5 %","97.5 %")]>=0)!=1, 3, 2)
  pred <- predict(grow_fits[[i]], newdata=newdata, re.form=NA, type='response')
  points(0:50, pred, col=cols[i], type='l', lwd=lwd, lty=lty)
}

segments(moist_percentiles[1], -10, moist_percentiles[1], 0.8, col='darkgrey')
abline(v=moist_percentiles[2], col='darkgrey')

arrows(moist_percentiles[1], -0.5, moist_percentiles[2], -0.5, 
       len=0.1, code=3, lwd=1.5)
text(mean(moist_percentiles), -0.43, 'Moisure percentiles
(10% and 90%)', cex=0.5)

legend('topleft', legend=sp_labs, text.font=3,
       cex=0.7, bty='n', lty=1, col=cols, lwd=2)

dev.off()



#############################################
### Plot the point estimates with 95% CIs ###
#############################################

# pdf("Figures/FigureS1.pdf", height=4, width=7)
# 
# par(mfrow=c(1,2), mar=c(6,4,1,1), oma=c(0.5,0.5,2,0.5))
# 
# pt.pch <- ifelse(apply(t1[t1$Variable=="Moisture",c('2.5 %','97.5 %')]<0, 1, 
#                        function(x) x[1]==x[2]),
#                  16, 21)
# plot(-as.numeric(as.character(t1$beta[t1$Variable=='Moisture'])), 
#      axes=F, xlab=NA,
#      ylab='Estimated coefficient',
#      pch=pt.pch, col=cols, cex=2, lwd=3,
#      ylim=c(-max(t1$`97.5 %`[t1$Variable=="Moisture"]),
#             -min(t1$`2.5 %`[t1$Variable=="Moisture"])),
#      xlim=c(0,9),
#      main="Moisture effect on Survival")
# segments(1:8, -as.numeric(as.character(t1$`97.5 %`[t1$Variable=='Moisture'])),
#          1:8, -as.numeric(as.character(t1$`2.5 %`[t1$Variable=='Moisture'])), 
#          col=cols, lwd=2)
# abline(h=0, lty=2)
# axis(1, labels=levels(growdf$species), at=1:8, las=2)
# axis(2)
# box()
# 
# 
# pt.pch <- ifelse(apply(t2[t2$Variable=="Moisture",c('2.5 %','97.5 %')]<0, 1, 
#                        function(x) x[1]==x[2]), 16, 21)
# plot(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
#      axes=F, xlab=NA,
#      ylab='Estimated coefficient',
#      pch=pt.pch,
#      col=cols, cex=2, 
#      ylim=range(growcis[rownames(growcis)=="moisture",]),
#      xlim=c(0,9),
#      main="Moisture effect on Growth", lwd=3)
# segments(1:8, as.numeric(as.character(t2$`2.5 %`[t2$Variable=='Moisture'])),
#          1:8, as.numeric(as.character(t2$`97.5 %`[t2$Variable=='Moisture'])),
#          col=cols, lwd=2)
# abline(h=0, lty=2)
# axis(1, labels=levels(growdf$species), at=1:8, las=2)
# axis(2)
# box()
# 
# 
# dev.off()







### USE FITTED MODELS TO PREDICT GROWTH / SURVIVAL AT DIFFERENT SOIL MOISTURE VALUES
gmat <- matrix(ncol=8, nrow=length(seq(0, 100, 0.1)))
rownames(gmat) <- seq(0, 100, 0.1)
for(i in 1:length(unique(growdf$species))){
  focsp <- levels(growdf$species)[i]
  tmpdf <- growdf[growdf$species == focsp,]
  newdata <- data.frame(moisture=seq(0, 100, 0.1),
                        densiometer=mean(tmpdf$densiometer),
                        startsize=mean(tmpdf$startsize))
  gmat[,i] <- predict(grow_fits[[i]], newdata=newdata, re.form=NA, type='response')
}

smat <- matrix(ncol=8, nrow=length(seq(0, 100, 0.1)))
rownames(smat) <- seq(0, 100, 0.1)
for(i in 1:length(unique(surv$Species))){
  focsp <- levels(surv$Species)[i]
  tmpdf <- surv[surv$Species == focsp,]
  newdata_moist <- data.frame(Moisture=seq(0, 100, 0.1),
                              Densiometer=mean(tmpdf$Densiometer),
                              Start_LA=mean(tmpdf$Start_LA, na.rm=T), Days=228, Status=1)
  smat[,i] <- predict(surv_fits_plot[[i]], newdata=newdata_moist, type='survival')
}



### COMPUTE VALUE OF SOIL MOISTURE WHERE SURVIVAL (ACROSS 0-100 SOIL MOISTURE) IS 50% MAX SURVIVAL
surv_tol <- vector()
for(sp in 1:8){
  s50 <- max(smat[,sp]) - (max(smat[,sp]) - min(smat[,sp]))/2
  surv_tol[sp] <- which.min(abs(smat[,sp] - s50))
}
names(surv_tol) <- levels(surv$Species)
surv_tol

# ### COMPUTE SURVIVAL TOLERANCE AS SURVIVAL AT LOW (ARBITRARY 15% SOIL MOISTURE)
# surv_15 <- round(smat['15',], 3)

### COMPUTE SURVIVAL TOLERANCE AS SURVIVAL AT LOW (lower 10%-ile SOIL MOISTURE)
round(moist_percentiles[1], 1)
surv_tol <- round(smat['13.5',], 3)

### COMPUTE GROWTH TOLERANCE AS GROWTH AT LOW (lower 10%-ile SOIL MOISTURE)
grow_tol <- round(gmat['13.5',], 3)

### COMPUTE SURVIVAL SENSITIVITY AS SURVIVAL AT HIGH minus LOW moisture (90%-ile minus 10%-ile SOIL MOISTURE)
round(moist_percentiles, 1)
surv_sens <- round(smat['39.1',] - smat['13.5',], 3)

### COMPUTE GROWTH SENSITIVITY AS GROWTH AT HIGH minus LOW moisture (90%-ile minus 10%-ile SOIL MOISTURE)
grow_sens <- round(gmat['39.1',] - gmat['13.5',], 3)

par(mfrow=c(2,2))
plot(grow_tol, grow_sens, bg=cols, pch=21, cex=1.5)
plot(surv_tol, surv_sens, bg=cols, pch=21, cex=1.5)
plot(grow_tol, surv_tol, bg=cols, pch=21, cex=1.5)
plot(grow_sens, surv_sens, bg=cols, pch=21, cex=1.5)




### GROW TOLERANCE AS PREDICTED GROWTH AT ZERO SOIL MOISTURE (INTERCEPT)
# grow_tol <- round(unlist(lapply(grow_fits, function(g) fixef(g)[1])), 3)
# grow_15 <- round(gmat[16,], 3)
# names(grow_15) <- levels(surv$Species)
# grow_15


# grow_beta <- as.numeric(t2$Estimate[t2$Variable=="Moisture"])
# names(grow_beta) <- levels(surv$Species)
# grow_beta

# surv_beta <- as.numeric(t1$beta[t1$Variable=="Moisture"])
# names(surv_beta) <- levels(surv$Species)
# surv_beta





par(mfrow=c(2,3), mar=c(4,4,1,1))

plot(grow_sens, grow_tol, pch=16, col=cols, cex=1.5)
cor.test(grow_sens, grow_tol)

plot(surv_sens, surv_tol, pch=16, col=cols, cex=1.5)
cor.test(surv_sens, surv_tol)

plot(grow_sens, surv_sens, pch=16, col=cols, cex=1.5)
cor.test(grow_sens, surv_sens)

plot(grow_tol, surv_sens, pch=16, col=cols, cex=1.5)
cor.test(grow_tol, surv_sens)

plot(grow_sens, surv_tol, pch=16, col=cols, cex=1.5)
cor.test(grow_sens, surv_tol)
     
plot(grow_tol, surv_tol, pch=16, col=cols, cex=1.5)
cor.test(grow_tol, surv_tol)

cor.test(grow_tol, surv_tol)
cor.test(grow_sens, surv_sens)




######################################################
######################################################
############### Trait data analysis ##################
######################################################
######################################################


### Get average survival by species
avg_surv <- round(100*(1-tapply(surv$Status, surv$Species, mean)),1)

growdf <- growdf[order(growdf$species, growdf$id, growdf$interval),]

### Get average survival by species
g <- vector()
for(i in 1:length(unique(growdf$id))){
  focdat <- growdf[growdf$id %in% unique(growdf$id)[i],]
  g[i] <- sum(focdat$growth)
  names(g)[i] <- as.character(focdat$species[1])
}
avg_grow <- tapply(g, names(g), mean)



head(trait)


# leaf area (cm2)
# leaf mass per area (LMA; kg m-2)
# leaf thickness (mm)
# leaf dry matter content (%)
# roots:shoot mass ratio (%)
# specific root length (SRL, cm g-1)
# root tissue density (RTD, g cm-3)
# total root system length (cm)
# average root system diameter (mm)
# maximum rooting depth (cm)
# number of root tips

trait$StemDryMassFraction[is.na(trait$StemDryMassFraction)] <- 0
trait$LeafDryMassFraction[is.na(trait$LeafDryMassFraction)] <- 0
trait$RootDryMassFraction[is.na(trait$RootDryMassFraction)] <- 0

trait$RootDryMass[is.na(trait$RootDryMass)] <- 0
trait$LeafDryMass[is.na(trait$LeafDryMass)] <- 0
trait$StemDryMass[is.na(trait$StemDryMass)] <- 0

trait$RSratio <- trait$RootDryMassFraction / 
  (trait$StemDryMassFraction + trait$LeafDryMassFraction)

trait$total_drymass <- trait$RootDryMass + trait$LeafDryMass + trait$StemDryMass

foctraits <- c("Plot",
               "Species",
               "ID",
               "Moisture",
               "LeafArea",
               "LMA", 
               "LeafThicknessMean",
               "LDMC",
               "RSratio",
               "SRL",
               "RTD",
               "RootLength",
               "RootAvgDiam",
               "RootDepth",
               "RootTips"#, "total_drymass"
               )

pctraits <- trait[,colnames(trait) %in% foctraits]
pctraits <- pctraits[!is.na(rowSums(pctraits[,!colnames(pctraits) %in% c("Plot",
                                                                         "Species",
                                                                         "ID",
                                                                         "Moisture")])),]
pctraits <- pctraits[,match(foctraits, names(pctraits))]

# par(mfrow=c(3,4))
# for(i in 1:ncol(pctraits)){
#   hist(pctraits[,i], xlab=names(pctraits)[i], main=NA)
#   hist(log10(pctraits[,i]), col=2, xlab=names(pctraits)[i], main=NA)
# }

pctraits$logLeafArea <- log10(pctraits$LeafArea)
pctraits$logLMA <- log10(pctraits$LMA)
pctraits$logLeafThicknessMean <- log10(pctraits$LeafThicknessMean)
pctraits$logRSratio <- log10(pctraits$RSratio)
pctraits$logSRL <- log10(pctraits$SRL)
pctraits$logRTD <- log10(pctraits$RTD)
pctraits$logRootLength <- log10(pctraits$RootLength)
pctraits$logRootAvgDiam <- log10(pctraits$RootAvgDiam)
pctraits$logRootDepth <- log10(pctraits$RootDepth)
pctraits$logRootTips <- log10(pctraits$RootTips)
#pctraits$logTotalDryMass <- log10(pctraits$total_drymass)

logfoctraits <- c("Plot",
                  "Species",
                  "ID",
                  "Moisture",
                  "logLeafArea",
                  "logLMA", 
                  "logLeafThicknessMean",
                  "LDMC",
                  "logRSratio",
                  "logSRL",
                  "logRTD",
                  "logRootLength",
                  "logRootAvgDiam",
                  "logRootDepth",
                  "logRootTips"
                  #,"logTotalDryMass"
                  )

pctraits <- pctraits[,names(pctraits) %in% logfoctraits]

### PCA with leaf and root traits
ord <- prcomp(~ logLeafArea +
                logLMA +
                logLeafThicknessMean +
                LDMC +
                logRSratio +
                logSRL +
                logRTD +
                logRootLength +
                logRootAvgDiam +
                logRootDepth +
                logRootTips #+ logTotalDryMass
              , center = TRUE, scale = TRUE, data=pctraits)

summary(ord) # percent of variance explained by each PC AXIS


### Varimax
scaled_data <- scale(pctraits[,-c(1:4)]) 

scaled_data <- scaled_data[,!colnames(scaled_data) %in% "logTotalDryMass"]

ord_new <- psych::principal(scaled_data, rotate="Varimax", nfactors=2, scores=TRUE)


d <- vegan::vegdist(ord_new$scores, method='gower')
adonis2(d ~ pctraits$Species)


##################################################################
### FIGURE OF ORDINATION
##################################################################
pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure3_v3.pdf",
    width=8, height=4)

par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(ord_new$scores[,1:2], pch=16, col="grey", cex=0.75, asp=1)
mtext("A", 3, -1.25, adj=0.05, font=2)
abline(v=0, h=0, lty=3)

l.x <- ord_new$weights[,1]*6.5
l.y <- ord_new$weights[,2]*6.5

arrows(rep(0,11), rep(0,11), l.x, l.y, len=0.1, lwd=1.25)

# Label position
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

arrownames <- c("LDMC",
                "Leaf area",
                "LMA",
                "Leaf thickness",
                "Root:Shoot",
                "SRL",
                "Root tissue density",
                "Root length",
                "Root diameter",
                "Root depth",
                "Root tips"#, "TotalMass"
                )

l.pos["logRootLength"] <- "2"
l.pos["logRootDepth"] <- "4"
l.pos["LDMC"] <- "3"
l.pos["logLMA"] <- "4"
l.pos["logRTD"] <- "1"

text(l.x, l.y, labels=arrownames, col="blue", pos=l.pos,, cex=0.75)

plot(ord_new$scores[,1:2], pch=16, col=scales::alpha(cols[pctraits$Species], 0.75), cex=0.75, asp=1)
abline(v=0, h=0, lty=3)

points(tapply(ord_new$scores[,1], pctraits$Species, mean),
       tapply(ord_new$scores[,2], pctraits$Species, mean), 
       pch=23, bg=cols, cex=1.2, lwd=1.25)

legend('topleft', legend=sp_labs, text.font=3, bty='n', pch=21, pt.bg=cols, cex=0.6, pt.cex=1.25)
mtext("B", 3, -1.25, adj=0.95, font=2)

dev.off()
##################################################################
##################################################################






##################################################################
### FIGURE OF SPECIES TRAITS VS. AVERAGE DEMOGRAPHIC RATES
##################################################################

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureX_new.pdf",
    width=10, height=11)

par(mfcol=c(6,4), mar=c(4,4,0.5,0.5))

for(i in 5:length(foctraits)){
  plot(tapply(pctraits[,logfoctraits[i]], pctraits$Species, mean, na.rm=T), avg_surv, 
       pch=16, col=cols, cex=1.5, main=logfoctraits[i])
  plot(tapply(pctraits[,logfoctraits[i]], pctraits$Species, mean, na.rm=T), avg_grow, 
       pch=17, col=cols, cex=1.5)
}

plot(tapply(ord_new$scores[,1], pctraits$Species, mean), avg_surv, 
     pch=16, col=cols, cex=1.5, main=logfoctraits[i])
plot(tapply(ord_new$scores[,2], pctraits$Species, mean), avg_grow, 
     pch=17, col=cols, cex=1.5)
dev.off()


sptrait_surv_cors <- sptrait_grow_cors <- list()

sptrait_surv_cors[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), avg_surv)
sptrait_surv_cors[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), avg_surv)

sptrait_grow_cors[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), avg_grow)
sptrait_grow_cors[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), avg_grow)

for(i in 5:length(foctraits)){
  sptrait_surv_cors[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                              pctraits$Species, 
                                              mean, na.rm=T), avg_surv)
  sptrait_grow_cors[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                              pctraits$Species, 
                                              mean, na.rm=T), avg_grow)
}

survcors <- round(do.call(rbind, lapply(sptrait_surv_cors, function(x) c(x[[4]], x[[3]]))), 2)
growcors <- round(do.call(rbind, lapply(sptrait_grow_cors, function(x) c(x[[4]], x[[3]]))), 2)
rownames(survcors) <- rownames(growcors) <- c('rc1', 'rc2', logfoctraits[-c(1:4)])
colnames(survcors) <- colnames(growcors) <- c('cor','p')


#################################################################################
### CORRELATIONS BETWEEN SPECIES-LEVEL SENSITIVTY / TOLERANCE AND TRAITS
#################################################################################

sptrait_surv_sens <- sptrait_grow_sens <- sptrait_surv_tol <- sptrait_grow_tol <- list()

sptrait_surv_sens[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), surv_sens)
sptrait_surv_tol[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), surv_tol)

sptrait_surv_sens[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), surv_sens)
sptrait_surv_tol[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), surv_tol)

sptrait_grow_sens[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), grow_sens)
sptrait_grow_tol[[1]] <- cor.test(tapply(ord_new$scores[,1], pctraits$Species, mean), grow_tol)

sptrait_grow_sens[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), grow_sens)
sptrait_grow_tol[[2]] <- cor.test(tapply(ord_new$scores[,2], pctraits$Species, mean), grow_tol)

for(i in 5:length(foctraits)){
  sptrait_surv_tol[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                              pctraits$Species, 
                                              mean, na.rm=T), surv_tol)
  sptrait_surv_sens[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                              pctraits$Species, 
                                              mean, na.rm=T), surv_sens)
  
  sptrait_grow_tol[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                              pctraits$Species, 
                                              mean, na.rm=T), grow_tol)
  sptrait_grow_sens[[i-2]] <- cor.test(tapply(pctraits[,logfoctraits[i]], 
                                             pctraits$Species, 
                                             mean, na.rm=T), grow_sens)
}


tol_sens <- cbind(round(do.call(rbind, 
                                lapply(sptrait_surv_sens, function(x) c(x[[4]], x[[3]]))), 2),
                  round(do.call(rbind, 
                                lapply(sptrait_surv_tol, function(x) c(x[[4]], x[[3]]))), 2),
                  round(do.call(rbind, 
                                lapply(sptrait_grow_sens, function(x) c(x[[4]], x[[3]]))), 2),
                  round(do.call(rbind, 
                                lapply(sptrait_grow_tol, function(x) c(x[[4]], x[[3]]))), 2))

rownames(tol_sens) <- c('rc1', 'rc2', logfoctraits[-c(1:4)])
colnames(tol_sens) <- paste(rep(c('surv_sens', 'surv_tol', 'grow_sens', 'grow_tol'), each=2),
                            rep(c('cor', 'p'),4), sep='-')

tol_sens
tol_sens[1:2,]

colnames(growcors) <- paste("grow_avg-", colnames(growcors), sep="")
colnames(survcors) <- paste("surv_avg-", colnames(survcors), sep="")


write.csv(cbind(growcors, survcors, tol_sens), 
          "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/TableX_demographic-traits-correlations_20230831.csv")
      



############################################################
##### FIGURE 4 ON DEMOGRAPHIC METRICS AND RC1 AND RC2 ######
############################################################

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure4_20230831.pdf", height=6, width=12)

par(mfcol=c(2,4), mar=c(4,5,1,1))

plot(tapply(ord_new$scores[,1], pctraits$Species, mean),
     surv_sens, pch=16, col=cols, cex=2, xlab="RC1", ylab="Survival sensitivity")
mtext("A", 3, -2, adj=0.95, font=2)
mod <- lm(surv_sens ~ tapply(ord_new$scores[,1], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

legend('topleft', legend=sp_labs, bty='n', pch=16, 
       col=cols, cex=1, text.font=3, pt.cex=2)

plot(tapply(ord_new$scores[,2], pctraits$Species, mean),
     surv_sens, pch=16, col=cols, cex=2, xlab="RC2", ylab="Survival sensitivity")
mtext("B", 3, -2, adj=0.05, font=2)
mod <- lm(surv_sens ~ tapply(ord_new$scores[,2], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,1], pctraits$Species, mean),
     surv_tol, pch=16, col=cols, cex=2, xlab="RC1", ylab="Survival tolerance")
mtext("C", 3, -2, adj=0.05, font=2)
mod <- lm(surv_tol ~ tapply(ord_new$scores[,1], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,2], pctraits$Species, mean),
     surv_tol, pch=16, col=cols, cex=2, xlab="RC2", ylab="Survival tolerance")
mtext("D", 3, -2, adj=0.05, font=2)
mod <- lm(surv_tol ~ tapply(ord_new$scores[,2], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,1], pctraits$Species, mean),
     grow_sens, pch=16, col=cols, cex=2, xlab="RC1", ylab="Growth sensitivity")
mtext("E", 3, -2, adj=0.05, font=2)
mod <- lm(grow_sens ~ tapply(ord_new$scores[,1], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,2], pctraits$Species, mean),
     grow_sens, pch=16, col=cols, cex=2, xlab="RC2", ylab="Growth sensitivity")
mtext("F", 3, -2, adj=0.05, font=2)
mod <- lm(grow_sens ~ tapply(ord_new$scores[,2], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,1], pctraits$Species, mean),
     grow_tol, pch=16, col=cols, cex=2, xlab="RC1", ylab="Growth tolerance")
mtext("G", 3, -2, adj=0.95, font=2)
mod <- lm(grow_tol ~ tapply(ord_new$scores[,1], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

plot(tapply(ord_new$scores[,2], pctraits$Species, mean),
     grow_tol, pch=16, col=cols, cex=2, xlab="RC2", ylab="Growth tolerance")
mtext("H", 3, -2, adj=0.95, font=2)
mod <- lm(grow_tol ~ tapply(ord_new$scores[,2], pctraits$Species, mean))
abline(mod, lty=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 3),
       col=ifelse(summary(mod)$coef[2,4] < 0.05, 1, 'grey'),
       lwd=ifelse(summary(mod)$coef[2,4] < 0.05, 2, 1))

dev.off()

############################################################
############################################################


i = 7
message(names(trait)[i])
plot(trait$Moisture, log(trait[,i]), col=cols[trait$Species])
cor.test(trait$Moisture, log(trait[,i]))

dim(trait)
dim(pctraits)

pctraits <- cbind(pctraits, ord_new$scores)

plot(moist, pctraits$RC1, col=cols[pctraits$Species])




plot(total_drymass[match(pctraits$ID, trait$ID)], 
     pctraits$RC1, col=cols[pctraits$Species], log='x')


moist <- trait$Moisture[!is.na(rowSums(trait[,colnames(trait) %in% foctraits][,-1]))]

itvtraits <- pctraits
itvtraits$moist <- moist
itvtraits <- itvtraits[!itvtraits$Species %in% "CECSCH",]
itvtraits$Species <- droplevels(itvtraits$Species)

itv_res <- data.frame()
itv_mods <- list()
for(t in 2:(ncol(itvtraits)-1)){
  itv_mods[[t-1]] <- list()
  names(itv_mods)[t-1] <- names(itvtraits)[t]
  for(sp in 1:length(levels(itvtraits$Species))){
    focdat <- itvtraits[itvtraits$Species == levels(itvtraits$Species)[sp],]
    itv_mods[[t-1]][[sp]] <- lm(focdat[,t] ~ focdat$moist)
    
    itv_res <- rbind(itv_res, data.frame(trait=names(itvtraits)[t], 
                              species=levels(itvtraits$Species)[sp], 
                              int=summary(itv_mods[[t-1]][[sp]])$coef[1,1],
                              slope=summary(itv_mods[[t-1]][[sp]])$coef[2,1],
                              p=summary(itv_mods[[t-1]][[sp]])$coef[2,4]))
  }
  names(itv_mods[[t-1]]) <- levels(itvtraits$Species)
}

itv_res[,3:5] <- round(itv_res[,3:5], 3)
itv_res[itv_res$p<0.05,]

plot(moist, pctraits$RC1, col=cols[pctraits$Species])
abline(lm(pctraits$RC1[pctraits$Species=="TETBAL"] ~ moist[pctraits$Species=="TETBAL"]))



############################################################
#### FIGURE 5 ON ITV OF RC1 AND RC2 ALONG THE SOIL MOISTURE GRADIENT
############################################################

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure5_new.pdf", height=5, width=9)

par(mfrow=c(1,2))
plot(pctraits$Moisture, pctraits$RC1, pch=16, 
     col=scales::alpha(cols[pctraits$Species], 0.75),
     xlab="Soil Moisture (%)", ylab="RC1", xlim=c(0,50))
for(i in 1:length(itv_mods$RC1)){
  p <- summary(itv_mods$RC1[[i]])$coef[2,4]
  abline(itv_mods$RC1[[i]], col=cols[1+i], 
         lty=ifelse(p<0.05, 1, 2), lwd=ifelse(p<0.05, 3, 1))
}
mtext("A", 3, -1.5, adj=0.05, font=2)

plot(moist, pctraits$RC2, pch=16, 
     col=scales::alpha(cols[pctraits$Species], 0.75),
     xlab="Soil Moisture (%)", ylab="RC2", xlim=c(0,50))
for(i in 1:length(itv_mods$RC2)){
  p <- summary(itv_mods$RC2[[i]])$coef[2,4]
  abline(itv_mods$RC2[[i]], col=cols[1+i], 
         lty=ifelse(p<0.05, 1, 2), lwd=ifelse(p<0.05, 3, 1))
}
mtext("B", 3, -1.5, adj=0.05, font=2)

legend('bottomleft', legend=sp_labs[-1], bty='n', pch=16, 
       col=cols[-1], cex=0.6, text.font=3, pt.cex=1)


dev.off()

############################################################
############################################################




### PLOT OF FINAL TOTAL DRY MASS 
trait$RootDryMass[is.na(trait$RootDryMass)] <- 0
trait$LeafDryMass[is.na(trait$LeafDryMass)] <- 0
trait$StemDryMass[is.na(trait$StemDryMass)] <- 0

total_drymass <- trait$RootDryMass + trait$LeafDryMass + trait$StemDryMass

plot(trait$Moisture, total_drymass, log='xy', 
     col=scales::alpha(cols[trait$Species], 0.75), pch=16)
summary(lm(total_drymass[trait$Species=="MANBID"] ~ trait$Moisture[trait$Species=="MANBID"]))

b <- boxplot(total_drymass ~ trait$Species, col=cols, log='y',
        ylab="Total dry mass (g)", xlab="Species", axes=F)
axis(2)
axis(1, labels=F, at=1:8)
legend('topleft', legend=sp_labs, bty='n', pch=15, 
       col=cols, cex=0.75, text.font=3, pt.cex=1)


foctraits <- c("Species",
               "LeafArea",
               "LMA", 
               "LeafThicknessMean",
               "LDMC",
               "RSratio",
               "SRL",
               "RTD",
               "RootLength",
               "RootAvgDiam",
               "RootDepth",
               "RootTips")

pctraits <- traits[!is.na(rowSums(trait[,-1])),]
pctraits <- pctraits[,match(foctraits, names(pctraits))]




###############################################
### Box plots of individual (logged) traits ###
###############################################
pdf("~/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS2_v2.pdf", height=9, width=8)

par(mfrow=c(6,4), mar=c(1,4,0.25,0.25), oma=c(5,0,0.5,0))

# splabs <- c("Cecropia schreberiana",
#             "Guarea guidonia",
#             "Inga laurina",
#             "Manilkara bidentata",
#             "Prestoea montana",
#             "schefflera morototoni",
#             "Tetragastris balsamifera",
#             "Urera baccifera")


for(i in 1:24){
  column <- c(5:20,24:31)[i]
  b <- boxplot(trait[,column] ~ trait$Species, col=cols, 
               log='y', axes=F, lwd=0.5, 
               ylab=NA, xlab=NA)
  axis(1, labels=F, at=1:8)
  axis(2)
  mtext(colnames(trait)[column], 2, 2, cex=0.5)
  mtext(LETTERS[i], 3, -0.75, at=1, cex=0.75)
  if(column > 27){
    # axis(1, labels=levels(trait$Species), at=1:8, las=2, cex.axis=0.5)
    axis(1, labels=splabs, at=1:8, las=2, font=3, cex.axis=0.5)
  }
}

dev.off()

############################################################################
### Linear regressions for individual traits vs. soil moisture (LOGGED!) ###
############################################################################
pdf("Figures/FigureS4.pdf", height=9, width=8)

par(mfrow=c(6,4), mar=c(1,4,0.25,0.25), oma=c(3,0,0.5,0))

for(i in 1:24){
  column <- c(5:20,24:31)[i]
  plot(trait$Moisture, trait[,column], col=cols[trait$Species],
       axes=F, lwd=0.5, ylab=NA, xlab=NA, log='y', cex=0.5)
  axis(1, labels=F)
  axis(2)
  mtext(colnames(trait)[column], 2, 2, cex=0.5)
  mtext(LETTERS[i], 3, -1, at=12, cex=0.75)
  if(column > 27){
    axis(1)
    mtext("Soil Moisture (%)", 1, 2.25, cex=0.75)
  }
  for(sp in 1:8){
    spdat <- trait[trait$Species==levels(trait$Species)[sp],]
    if(sum(!is.na(spdat[,column]))>5){
      fit <- lm(log10(spdat[,column]) ~ spdat$Moisture)
      if(summary(fit)$coefficients[2,4] <= 0.05){
        abline(fit, col=cols[sp], lwd=2)
      }
    }
  }
}

dev.off()


#####################################
### PCA with leaf and root traits ###
#####################################

### Remove one UREBAC outlier
trait <- trait[!(trait$Plot==12 & trait$Position==11),]

### PCA with leaf and root traits
ord <- prcomp(~ LeafArea +
                LMA +
                LDMC +
                LeafThicknessMean +
                RootDryMassFraction +
                RootLength +
                RootDepth +
                RootAvgDiam +
                RootTips +
                SRL +
                RTD,
              center = TRUE, scale. = TRUE, data=trait)

biplot(ord)
summary(ord)

cor(log(trait$LMA), ord$x[,1], use='p')



### visualize the variables used in the PCA
a1 <- fviz_pca_var(ord, col.var="contrib",
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE)

###  visualize results for individuals
b1 <- fviz_pca_ind(ord,
                   label = "none", # hide individual labels
                   habillage = trait$Species, # color by groups
                   palette = cols,
                   # Concentration ellipses
                   addEllipses = T)

### Put the plots together
PCAplot1 <- ggarrange(a1, b1, ncol = 2, nrow = 1)

########################
### Plot PCA Results ###
########################
# pdf("Figures/Figure3.pdf", height=5, width=10)
# PCAplot1
# dev.off()
########################

### Contribution of variables
tres.var <- get_pca_var(ord)
contrib1 <- tres.var$contrib
contrib1

### Extract PCA coordinates for individuals and plotting it against soil moisture
trait$PCA1 <- get_pca_ind(ord)$coord[,1]
trait$PCA2 <- get_pca_ind(ord)$coord[,2]


###############################################
### Plot PCA axis 1 and 2 vs. Soil Moisture ###
###############################################

pdf("Figures/FigureS3.pdf", height=5, width=10)

par(mfrow=c(1,2))
plot(trait$Moisture, trait$PCA1, pch=21,
     col=cols[trait$Species], xlim=c(0,50),
     xlab="Soil Moisture (%)",
     ylab="PC Dim 1")
mtext("A", 3, -1.2, at=3, cex=1.25)
for(sp in 1:8){
  spdat <- trait[trait$Species==levels(trait$Species)[sp],]
  if(sum(!is.na(spdat$PCA1))>5){
    fit <- lm(spdat$PCA1 ~ spdat$Moisture)
    if(summary(fit)$coefficients[2,4] <= 0.05){
      abline(fit, col=cols[sp], lwd=2)
    }
  }
}

legend('bottomleft', legend=levels(trait$Species),
       cex=0.7, bty='n', col=cols, pch=21, pt.lwd=1.5, pt.cex=1.25)

plot(trait$Moisture, trait$PCA2, pch=21, 
     col=cols[trait$Species], xlim=c(0,50),
     xlab="Soil Moisture (%)",
     ylab="PC Dim 2")
mtext("B", 3, -1.2, at=3, cex=1.25)
for(sp in 1:8){
  spdat <- trait[trait$Species==levels(trait$Species)[sp],]
  if(sum(!is.na(spdat$PCA2))>5){
    fit <- lm(spdat$PCA2 ~ spdat$Moisture)
    if(summary(fit)$coefficients[2,4] <= 0.05){
      abline(fit, col=cols[sp], lwd=2)
    }
  }
}

dev.off()

###############################################



#############################################################################
### Plot PCA axis 1 and 2 vs. Soil Moisture effect on Survival and Growth ###
#############################################################################

pdf("Figures/Figure4.pdf")

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(-t1$beta[t1$Variable=='Moisture'], 
     tapply(trait$PCA1, trait$Species, mean, na.rm=T),
     bg=cols, pch=21,
     xlim=c(-0.01,0.65),
     ylim=c(-6,6),
     cex=2,
     xlab="Moisture effect on Survival",
     ylab="PC Dim 1")
legend('topright', legend=levels(surv$Species),
       cex=0.75, bty='n', pch=21, pt.bg = cols, pt.cex=1)
mtext("A", 3, -1.2, at=-0.01)
segments(-t1$`2.5 %`[t1$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, mean, na.rm=T),
         -t1$`97.5 %`[t1$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, mean, na.rm=T), 
         col=cols, lwd=2)
segments(-t1$beta[t1$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, quantile, 0.025, na.rm=T),
         -t1$beta[t1$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, quantile, 0.975, na.rm=T),
         col=cols, lwd=2)
abline(h=0, v=0, lty=3)

plot(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
     tapply(trait$PCA1, trait$Species, mean, na.rm=T),
     bg=cols, pch=21,
     xlim=c(-0.005,.05),
     ylim=c(-6,5),
     cex=2,
     xlab="Moisture effect on LA Growth",
     ylab="PC Dim 1")
mtext("B", 3, -1.2, at=-0.0045)
segments(t2$`2.5 %`[t2$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, mean, na.rm=T),
         t2$`97.5 %`[t2$Variable=='Moisture'], 
         tapply(trait$PCA1, trait$Species, mean, na.rm=T),
         col=cols, lwd=2)
segments(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
         tapply(trait$PCA1, trait$Species, quantile, 0.025, na.rm=T),
         as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
         tapply(trait$PCA1, trait$Species, quantile, 0.975, na.rm=T), 
         col=cols, lwd=2)
abline(h=0, v=0, lty=3)

plot(-t1$beta[t1$Variable=='Moisture'], 
     tapply(trait$PCA2, trait$Species, mean, na.rm=T),
     bg=cols, pch=21,
     xlim=c(-0.01,0.65),
     ylim=c(-4,5),
     cex=2,
     xlab="Moisture effect on Survival",
     ylab="PC Dim 2")
mtext("C", 3, -1.2, at=-0.01)
segments(-t1$`2.5 %`[t1$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, mean, na.rm=T),
         -t1$`97.5 %`[t1$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, mean, na.rm=T), 
         col=cols, lwd=2)
segments(-t1$beta[t1$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, quantile, 0.025, na.rm=T),
         -t1$beta[t1$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, quantile, 0.975, na.rm=T),
         col=cols, lwd=2)
abline(h=0, v=0, lty=3)

plot(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
     tapply(trait$PCA2, trait$Species, mean, na.rm=T),
     bg=cols, pch=21,
     xlim=c(-0.005,.05),
     ylim=c(-3,5),
     cex=2,
     xlab="Moisture effect on LA Growth",
     ylab="PC Dim 2")
mtext("D", 3, -1.2, at=-0.0045)
segments(t2$`2.5 %`[t2$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, mean, na.rm=T),
         t2$`97.5 %`[t2$Variable=='Moisture'], 
         tapply(trait$PCA2, trait$Species, mean, na.rm=T),
         col=cols, lwd=2)
segments(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
         tapply(trait$PCA2, trait$Species, quantile, 0.025, na.rm=T),
         as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
         tapply(trait$PCA2, trait$Species, quantile, 0.975, na.rm=T), 
         col=cols, lwd=2)
abline(h=0, v=0, lty=3)

dev.off()





############################################
### Look at change in WUE (Isotope data) ###
############################################

pdf("Figures/FigureS5_v2.pdf")

par(mar=c(5,5,1,1))
plot(nutrient$Moisture, nutrient$iWUE_Lambers, pch=21,
     bg=cols[nutrient$Species], xlim=c(0,45),
     xlab="Soil Moisture (%)",
     ylab=bquote("i"*italic("WUE")~"("*plain(mu)~mol~mol^-1*")"))

splabs <- c("G. guidonia",
            "I. laurina",
            "M . bidentata",
            "P. montana",
            "S. morototoni",
            "T. balsamifera")

legend('bottomleft', legend=splabs, text.font=3, #sort(unique(nutrient$Species)),
       cex=0.75, bty='n', pch=21, 
       pt.bg = cols[which(levels(nutrient$Species) %in% sort(unique(nutrient$Species)))], 
       pt.cex=1.25, pt.lwd=1.25)

t3 <- data.frame()
wue_mods <- list()
for(sp in 1:8){
  focsp <- levels(nutrient$Species)[sp]
  focnut <- nutrient[nutrient$Species %in% focsp,]
  if(nrow(focnut)>2){
    wue_mods[[sp]] <- lm(focnut$iWUE_Lambers ~ focnut$Moisture)
    lty <- ifelse(summary(wue_mods[[sp]])$coef[2,4] < 0.05, 1, 2)
    lwd <- ifelse(summary(wue_mods[[sp]])$coef[2,4] < 0.05, 2, 1)
    segments(min(focnut$Moisture), 
             predict(wue_mods[[sp]])[which(focnut$Moisture==min(focnut$Moisture))],
             max(focnut$Moisture),
             predict(wue_mods[[sp]])[which(focnut$Moisture==max(focnut$Moisture))], 
             col=cols[sp], lty=lty, lwd=lwd)
    t3 <- rbind(t3, as.data.frame(cbind(focsp, round(summary(wue_mods[[sp]])$coef, 4))))
  }
}


dev.off()



t3 <- t3[grepl("Moisture", rownames(t3)),]
rownames(t3) <- NULL

write.csv(t3, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/TableS1_WUE_fits.csv", row.names = F)





########################################
### Look at PRE DAWN WATER POTENTIAL ###
########################################

# Pre-dawn leaf water potential measurements from some seedlings
wp <- read.csv(paste0(path, "LUQ_DroughtExp_WaterPotential.csv"))

# Remove a dying plant
wp <- wp[!wp$Notes %in% 'Seedling dying',]

# Get mean value per plant of PLWP (of several measurements made per plant)
wp$WP_mean <- apply(wp[,-c(1:2, 7:8)], 1, mean, na.rm=T)

# Add other variables to water potential data 
wp$ID <- paste(wp$Plot, wp$Position, sep=".")
wp$Species <- as.factor(grow$Species[match(wp$ID, grow$ID)])
wp$Moisture <- grow$Moisture[match(wp$ID, grow$ID)]


### PLOT IT!!!

pdf("Figures/FigureSx_water_potential_v1.pdf")

par(mfrow=c(1,1), mar=c(5,5,1,1))

plot(wp$Moisture, wp$WP_mean, 
     pch=16, col=scales::alpha(cols[wp$Species], 0.6), 
     cex=1.5,  xlim=c(0,45),
     ylab=bquote(Psi['PLWP']~(MPa)),
     xlab='Soil moisture (%)', cex.lab=1.25)
points(tapply(wp$Moisture, wp$Plot, mean, na.rm=T), 
       tapply(wp$WP_mean, wp$Plot, mean, na.rm=T),
       pch=25, lwd=3, cex=1.5)

legend('topleft', legend=c("C. schreberiana",
                           "G. guidonia",
                           "I. laurina",
                           "M. bidentata",
                           "P. montana",
                           "S. morototoni",
                           "T. balsamifera"),
       pch=16, col=cols[1:8], bty='n', pt.cex = 1.5, text.font=3, cex=0.8)

wp <- droplevels(wp[!is.na(wp$Species),])

dev.off()



pdf("Figures/FigureSx_water_potential_v2.pdf")

par(mfrow=c(1,1), mar=c(5,5,1,1))

plot(wp$Moisture, wp$WP_mean, 
     pch=16, col=scales::alpha(cols[wp$Species], 0.6), 
     cex=1.5,  xlim=c(0,45),
     ylab=bquote(Psi['PDWP']~(MPa)),
     xlab='Soil moisture (%)', cex.lab=1.25)
points(tapply(wp$Moisture, wp$Plot, mean, na.rm=T), 
       tapply(wp$WP_mean, wp$Plot, mean, na.rm=T),
       pch=25, lwd=3, cex=1.5)

legend('topleft', legend=c("C. schreberiana",
                           "G. guidonia",
                           "I. laurina",
                           "M. bidentata",
                           "P. montana",
                           "S. morototoni",
                           "T. balsamifera"),
       pch=16, col=cols[1:8], bty='n', pt.cex = 1.5, text.font=3, cex=0.8)

wp <- droplevels(wp[!is.na(wp$Species),])

dev.off()









