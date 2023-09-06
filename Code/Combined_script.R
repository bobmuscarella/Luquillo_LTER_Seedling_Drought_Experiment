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
require(factoextra)
require(vegan)
require(plotfunctions)
require(ggplot2)
require(RColorBrewer)
require(scales)
require(dplyr)


### READ DATA FROM GITHUB
path <- "https://raw.github.com/bobmuscarella/Luquillo_LTER_Seedling_Drought_Experiment/master/Data/"

# Full soil moisture measurements
moist <- read.csv(paste0(path, "LUQ_DroughtExp_Soil_moisture_complete.csv"))

# Growth data
grow <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_growth.csv"))
grow$Date <- as.Date(as.character(grow$Date), format="%m/%d/%y")
grow$ID <- paste(grow$Plot, grow$Position, sep = '.')
grow$Species <- as.factor(grow$Species)

# Survival data
surv <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_survival.csv"))
surv$ID <- paste(surv$Plot, surv$Position, sep = '.')
surv$Species <- as.factor(surv$Species)

# Trait data
trait <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_traits.csv"))
trait$ID <- paste(trait$Plot, trait$Position, sep = '.')
trait$Species <- as.factor(trait$Species)

# Photosynthesis data
# photo <- read.csv(paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv"))
# photo$ID <- paste(photo$Plot, photo$Position, sep = '.')

# Nutrient data
nutrient <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_leafnutrients.csv"))
# Add empty factor levels to keep species order consistent with other datasets 
# (Because not all species were included in the nutrient analysis)
nutrient$Species <- factor(nutrient$Species, levels = levels(grow$Species))
nutrient$ID <- paste(nutrient$Plot, nutrient$Position, sep = '.')

# Pre-dawn water potential data
# pdwp


### SET A COLOR PALETTE TO IDENTIFY SPECIES
cols <- brewer.pal(8, "Dark2")


### LOAD HELPER FUNCTION (to extract summary table from Cox ME models)
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
       main=focsp, las=2)
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
    text(as.Date("2019-11-15"), 33.5, "Soil", cex=0.8, font=2, pos=3)
    text(as.Date("2019-11-15"), 32, "Moisture %", cex=0.8, font=2, pos=3)
  }
}

# dev.off()


################################
### Plot raw survival curves ###
################################

# pdf("Figures/FigureS6.pdf", height=6, width=8)

surv_curves_fit <- survfit(Surv(Days, Status) ~ Species, data = surv)  
plot(surv_curves_fit, col=cols, lwd=2, axes=F, xlim=c(0,250), 
     xlab="Days", ylab="Percent Surviving")
legend('bottomleft', legend=sort(unique(surv$Species)),
       cex=0.7, bty='n', lty=1, col=cols, lwd=2)
axis(2, at=seq(0,1,by=0.2), labels=seq(0,100,by=20))
axis(1)

# dev.off()


#########################
### Survival analysis ###
#########################

### Run species-by-species survival models in a loop
surv_fits <- surv_fits_plot <- list()
surv_fits_int <- list()
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
  # surv_fits_plot[[i]] <- coxph(Surv(Days, Status) ~ Moisture + Densiometer + Moisture:Densiometer
                               # + Start_LA, data=focdat, model=T)
  surv_fits_plot[[i]] <- coxph(Surv(Days, Status) ~ Moisture + Densiometer + Moisture*Densiometer
                               + Start_LA, data=focdat, model=T)
  
  # Name the items in the resulting list
  names(surv_fits)[i] <- levels(surv$Species)[i]
  names(surv_fits_int)[i] <- levels(surv$Species)[i]
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

# Fit the growth models
grow_fits_int <- grow_fits <- grow_fits2 <- list()
for(sp in 1:length(levels(growdf$species))){
  print(levels(growdf$species)[sp])
  tmpdf <- growdf[growdf$species %in% levels(growdf$species)[sp],]
  grow_fits_int[[sp]] <- lmer(growth ~ moisture + densiometer 
                              + moisture*densiometer + startsize 
                              + (1|plot) + (1|indv), data=tmpdf)

  grow_fits[[sp]] <- lmer(growth ~ moisture + densiometer + startsize 
                           + (1|plot) + (1|indv), data=tmpdf)
  
  grow_fits2[[sp]] <- lmer(growth ~ moisture + densiometer + startsize 
                          + (1|plot), data=tmpdf)
  
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

library(MuMIn)
r.squaredGLMM(grow_fits[[1]])

# write.csv(t2, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/Table2_growth_summary_20211103.csv", row.names = F)


# pdf("/Users/au529793/Desktop/growth_interactions.pdf", height=8, width=6)

par(mfrow=c(4,2), mar=c(4,5,1,1))

for(i in 1:length(grow_fits)){
tmp <- grow_fits[[i]]

ndlolight <- data.frame(moisture=seq(0,100,length.out=100), 
                      densiometer=rep(3, 100),
                      startsize=rep(50,100))

ndhilight <- data.frame(moisture=seq(0,100,length.out=100), 
                      densiometer=rep(15, 100),
                      startsize=rep(50,100),
                      plot=rep(1,100))

lolight_pred <- predict(tmp, newdata=ndlolight,  re.form=NA)
hilight_pred <- predict(tmp, newdata=ndhilight,  re.form=NA)

plot(seq(0,100,length.out=100), lolight_pred, type='l', 
     ylim=range(c(hilight_pred, lolight_pred)), lwd=3,
     xlab='Soil moisture (%)', ylab='Pred. growth', main=names(grow_fits)[i])
lines(seq(0,100,length.out=100), hilight_pred, col=2, lwd=3)
legend('topleft', legend=c('Low light (3%)', 'High light (15%)'), 
       lty=1, col=1:2, bty='n', lwd=3)
}

dev.off()
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




###########################################
### Figure Showing Effects on Survival and Growth ###
###########################################

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


# pdf("Figures/Figure2.pdf", width=7, height=4)

par(mfcol=c(1,2), mar=c(5,5,3,0.5))

plot(0:50, ylim=c(0,100), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab="Pred. Survival (%)")
mtext("A", 3, -1.5, at=3, cex=1.5)
for(i in 1:length(pred_moist)){
  lty <- ifelse(sum(t1[t1$Variable=='Moisture', c("2.5 %","97.5 %")][i,]>=0)!=1, 1, 2)
  lwd <- ifelse(sum(t1[t1$Variable=='Moisture', c("2.5 %","97.5 %")][i,]>=0)!=1, 3, 2)
  points(0:50, 100*pred_moist[[i]], col=cols[i], type='l', lwd=lwd, lty=lty)
}

plot(0:50, ylim=c(-0.75,1.5), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab=bquote("Pred. Leaf Area Growth (cm"^2~"day"^-1*")"))
mtext("B", 3, -1.5, at=50, cex=1.5)
polygon(c(-10,70,70,-10), c(0,0,-2,-2), col='grey', lty=0)
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

legend('topleft', legend=levels(surv$Species),
       cex=0.7, bty='n', lty=1, col=cols, lwd=2)

dev.off()






# pdf("/Users/au529793/Desktop/survival_interactions.pdf", height=8, width=6)

par(mfrow=c(4,2), mar=c(4,5,1,1))

for(i in 1:length(pred_moist_hilight)){
  plot(0:50, ylim=c(0,100), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab="Pred. Survival (%)", 
     main=names(pred_moist_hilight)[i])
  lines(0:50, 100*pred_moist_hilight[[i]], col=2, lwd=3)
  lines(0:50, 100*pred_moist_lolight[[i]], col=1, lwd=3)
  legend('topleft', legend=c('Low light (3%)', 'High light (15%)'), 
         lty=1, col=1:2, bty='n', lwd=3)
}

dev.off()











#############################################
### Plot the point estimates with 95% CIs ###
#############################################

pdf("Figures/FigureS1.pdf", height=4, width=7)

par(mfrow=c(1,2), mar=c(6,4,1,1), oma=c(0.5,0.5,2,0.5))

pt.pch <- ifelse(apply(t1[t1$Variable=="Moisture",c('2.5 %','97.5 %')]<0, 1, 
                       function(x) x[1]==x[2]),
                 16, 21)
plot(-as.numeric(as.character(t1$beta[t1$Variable=='Moisture'])), 
     axes=F, xlab=NA,
     ylab='Estimated coefficient',
     pch=pt.pch, col=cols, cex=2, lwd=3,
     ylim=c(-max(t1$`97.5 %`[t1$Variable=="Moisture"]),
            -min(t1$`2.5 %`[t1$Variable=="Moisture"])),
     xlim=c(0,9),
     main="Moisture effect on Survival")
segments(1:8, -as.numeric(as.character(t1$`97.5 %`[t1$Variable=='Moisture'])),
         1:8, -as.numeric(as.character(t1$`2.5 %`[t1$Variable=='Moisture'])), 
         col=cols, lwd=2)
abline(h=0, lty=2)
axis(1, labels=levels(growdf$species), at=1:8, las=2)
axis(2)
box()


pt.pch <- ifelse(apply(t2[t2$Variable=="Moisture",c('2.5 %','97.5 %')]<0, 1, 
                       function(x) x[1]==x[2]), 16, 21)
plot(as.numeric(as.character(t2$Estimate[t2$Variable=='Moisture'])), 
     axes=F, xlab=NA,
     ylab='Estimated coefficient',
     pch=pt.pch,
     col=cols, cex=2, 
     ylim=range(growcis[rownames(growcis)=="moisture",]),
     xlim=c(0,9),
     main="Moisture effect on Growth", lwd=3)
segments(1:8, as.numeric(as.character(t2$`2.5 %`[t2$Variable=='Moisture'])),
         1:8, as.numeric(as.character(t2$`97.5 %`[t2$Variable=='Moisture'])),
         col=cols, lwd=2)
abline(h=0, lty=2)
axis(1, labels=levels(growdf$species), at=1:8, las=2)
axis(2)
box()


dev.off()






grow_fits[[1]]@resp












######################################################
######################################################
############### Trait data analysis ##################
######################################################
######################################################

###############################################
### Box plots of individual (logged) traits ###
###############################################
pdf("~/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS2_v2.pdf", height=9, width=8)

par(mfrow=c(6,4), mar=c(1,4,0.25,0.25), oma=c(5,0,0.5,0))

splabs <- c("Cecropia schreberiana",
            "Guarea guidonia",
            "Inga laurina",
            "Manilkara bidentata",
            "Prestoea montana",
            "schefflera morototoni",
            "Tetragastris balsamifera",
            "Urera baccifera")


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
              center = TRUE, scale = TRUE, data=trait)

summary(ord)

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
            "M. bidentata",
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









