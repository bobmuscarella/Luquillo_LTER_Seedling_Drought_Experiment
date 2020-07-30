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
grow$Date <- as.Date(as.character(grow$Date), format="%m/%d/%y")
grow$ID <- paste(grow$Plot, grow$Position, sep = '.')

surv <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_survival.csv"))
surv$ID <- paste(surv$Plot, surv$Position, sep = '.')

trait <- read.csv(paste0(path, "LUQ_DroughtExp_Seedling_traits.csv"))
trait$ID <- paste(trait$Plot, trait$Position, sep = '.')

photo <- read.csv(paste0(path, "LUQ_DroughtExp_Survey_photosynthesis.csv"))
photo$ID <- paste(photo$Plot, photo$Position, sep = '.')

# Load packages
require(survival)
require(survminer)
require(coxme)
require(ggplot2)
require(RColorBrewer)
library(lme4)
library(scales)
library(ggplot2)
library(factoextra)

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
  
  ############################################################
  ### EXCLUDE INDIVIDUALS WHO DIED IN THE FIRST 2 INTERVALS...
  # focdat <- focdat[focdat$Days>14,]
  ############################################################
  
  # Run the models with plot random effect (we tested other random effects but AIC selected this)
  surv_fits[[i]] <- coxme(Surv(Days, Status) ~ Moisture + Densiometer + Start_LA
                             + (1|Plot), data=focdat)
  # Run models without random effect in order to be able to plot predicted response
  # (not currently possible with coxme models)
  surv_fits_plot[[i]] <- coxph(Surv(Days, Status) ~ Moisture + Densiometer + Start_LA, 
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

write.csv(t1, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/Table1_survival_summary.csv", row.names = F)


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
  tmpdf <- data.frame(id, species, growth, moisture, densiometer, maysize, startsize, plot, indv, interval)
  growdf <- rbind(growdf, tmpdf)
  growdf <- growdf[!is.na(growdf$growth),]
}

# Fit the growth models
grow_fits <- list()

for(sp in 1:length(levels(growdf$species))){
  print(levels(growdf$species)[sp])
  tmpdf <- growdf[growdf$species %in% levels(growdf$species)[sp],]
  
  ############################################################
  ### Remove the first two intervals to see how that affects the starting size effect
  # tmpdf <- tmpdf[tmpdf$interval>1,]
  ############################################################
  
  grow_fits[[sp]] <- lmer(growth ~ moisture + densiometer + startsize 
                          + (1|plot), data=tmpdf)
  names(grow_fits)[sp] <- levels(growdf$species)[sp]
}

# Extract coefficients and compute confidence intervals
Species <- rep(names(grow_fits), each=4)
Variable <- rep(c("Intercept","Moisture","Densiometer","Start_size"), 8)
# Estimate <- round(do.call(c, lapply(grow_fits, coef)),3)
Estimate <- round(do.call(c, lapply(grow_fits, fixef)),3)
# growcis <- do.call(rbind, lapply(grow_fits, function(x) round(confint(x),3)))
growcis <- do.call(rbind, lapply(grow_fits, function(x) round(confint(x)[-(1:2),],3)))

t2 <- as.data.frame(cbind(Species, Variable, Estimate, growcis))
t2$`2.5 %` <- as.numeric(as.character(t2$`2.5 %`))
t2$`97.5 %` <- as.numeric(as.character(t2$`97.5 %`))
rownames(t2) <- NULL

write.csv(t2, "/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Tables/Table2_growth_summary.csv", row.names = F)


######################################################################################
### Plot individual growth trajectories over time with soil moisture and mortality ###
######################################################################################

grow$soil_col <- brewer.pal(10, "Spectral")[cut(grow$Moisture, 10)]

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure1.pdf", height=6, width=8)

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
    legend("topright", legend=c("High","Mid","Low"), 
           col=brewer.pal(10, "Spectral")[c(9,5,2)], 
           bty="n", lty=1, lwd=2, title="Soil Moisture")
  }
}

dev.off()

###########################################
### Plot Effects on Survival and Growth ###
###########################################

# Predict survival as a function of soil moisture by species
pred_moist <- list()
pred_light <- list()
pred_size <- list()

for (i in seq_along(levels(surv$Species))){
  focdat <- surv[surv$Species %in% levels(surv$Species)[i],]
  newdata_moist <- data.frame(Moisture = 0:50, 
                              Densiometer=mean(focdat$Densiometer),
                              Start_LA=mean(focdat$Start_LA, na.rm=T), Days=228, Status=1)
  pred_moist[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_moist, type='expected'))
  newdata_light <- data.frame(Moisture=mean(focdat$Moisture),
                              Densiometer=0:25,
                              Start_LA=mean(focdat$Start_LA, na.rm=T), Days=228, Status=1)
  pred_light[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_light, type='expected'))
  newdata_size <- data.frame(Moisture=mean(focdat$Moisture),
                             Densiometer=mean(focdat$Densiometer),
                             Start_LA=seq(min(surv$Start_LA, na.rm=T), 
                                        max(surv$Start_LA, na.rm=T), 
                                        length.out = 51), Days=228, Status=1)
  pred_size[[i]] <- exp(-predict(surv_fits_plot[[i]], newdata=newdata_size, type='expected'))
}


pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure2.pdf", width=7, height=4)

par(mfcol=c(1,2), mar=c(5,5,3,0.5))

plot(0:50, ylim=c(0,100), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab="Pred. Survival (%)")
mtext("Survival", 3, 1)
for(i in 1:length(pred_moist)){
  lty <- ifelse(sum(t1[t1$Variable=='Moisture', c("2.5 %","97.5 %")][i,]>=0)!=1, 1, 2)
  lwd <- ifelse(sum(t1[t1$Variable=='Moisture', c("2.5 %","97.5 %")][i,]>=0)!=1, 3, 2)
  points(0:50, 100*pred_moist[[i]], col=cols[i], type='l', lwd=lwd, lty=lty)
}

plot(0:50, ylim=c(-0.75,1.5), pch=NA,
     xlab="Soil Moisture (%)", 
     ylab=bquote("Pred. Growth (cm"^2~"day"^-1*")"))
mtext("Growth", 3, 1)
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


#############################################
### Plot the point estimates with 95% CIs ###
#############################################

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS1.pdf", height=4, width=7)


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
                       function(x) x[1]==x[2]),
                 16, 21)
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


######################################################
######################################################
######################################################
############### Trait data analysis ##################
######################################################
######################################################
######################################################


### Boxplots of individual traits
pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS2.pdf", height=9, width=8)

par(mfrow=c(6,4), mar=c(1,4,0.25,0.25), oma=c(3,0,0,0))

for(i in c(5:20,24:31)){
  b <- boxplot(trait[,i] ~ trait$Species, col=cols, log='y', axes=F, lwd=0.5, 
               ylab=NA, xlab=NA)
  axis(1, labels=F, at=1:8)
  axis(2)
  mtext(colnames(trait)[i], 2, 2, cex=0.5)
  if(i>27){
    axis(1, labels=levels(trait$Species), at=1:8, las=2, cex.axis=0.5)
  }
}

dev.off()


### Linear regressions for individual traits vs. soil moisture (LOGGED!)
pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS3.pdf", height=9, width=8)

par(mfrow=c(6,4), mar=c(1,4,0.25,0.25), oma=c(3,0,0.5,0))

for(i in c(5:20,24:31)){
  plot(trait$Moisture, trait[,i], col=cols[trait$Species],
       axes=F, lwd=0.5, ylab=NA, xlab=NA, log='y', cex=0.5)
  axis(1, labels=F)
  axis(2)
  mtext(colnames(trait)[i], 2, 2, cex=0.5)
  if(i>27){
    axis(1)
    mtext("Soil Moisture (%)", 1, 2.25, cex=0.75)
  }
  for(sp in 1:8){
    spdat <- trait[trait$Species==levels(trait$Species)[sp],]
    if(sum(!is.na(spdat[,i]))>5){
      fit <- lm(log10(spdat[,i]) ~ spdat$Moisture)
      if(summary(fit)$coefficients[2,4] <= 0.05){
        abline(fit, col=cols[sp], lwd=2)
      }
    }
  }
}
dev.off()




##########################################
### PCA with leaf and root traits
##########################################

### remove NA values and 2 urebac outliers (56,94)
# PCA_data <- trait[-c(56,94,194,195,324,339,401,497),]

trait <- trait[trait$Plot!=12 & trait$Position!=11,]

### PCA with only leaf and root traits
library(vegan)

ord <- prcomp(~ LeafArea  +
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

pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure3.pdf", height=5, width=10)
PCAplot1
dev.off()


### Contribution of variables
tres.var <- get_pca_var(ord)
contrib1 <- tres.var$contrib
contrib1

### Extract PCA coordinates for individuals and plotting it against soil moisture
trait$PCA1 <- get_pca_ind(ord)$coord[,1]
trait$PCA2 <- get_pca_ind(ord)$coord[,2]


pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/FigureS3.pdf", height=5, width=10)

par(mfrow=c(1,2))
plot(trait$Moisture, trait$PCA1, pch=16,
     col=cols[trait$Species], xlim=c(0,50), cex=0.5,
     xlab="Soil Moisture (%)",
     ylab="PCA 1")
for(sp in 1:8){
  spdat <- trait[trait$Species==levels(trait$Species)[sp],]
  if(sum(!is.na(spdat[,i]))>5){
    fit <- lm(spdat$PCA1 ~ spdat$Moisture)
    if(summary(fit)$coefficients[2,4] <= 0.05){
      abline(fit, col=cols[sp], lwd=2)
    }
  }
}

legend('bottomleft', legend=levels(trait$Species),
       cex=0.8, bty='n', col=cols, pch=16, pt.lwd=1.5)

plot(trait$Moisture, trait$PCA2, pch=16, 
     col=cols[trait$Species], xlim=c(0,50), cex=0.5,
     xlab="Soil Moisture (%)",
     ylab="PCA 2")
for(sp in 1:8){
  spdat <- trait[trait$Species==levels(trait$Species)[sp],]
  if(sum(!is.na(spdat[,i]))>5){
    fit <- lm(spdat$PCA2 ~ spdat$Moisture)
    if(summary(fit)$coefficients[2,4] <= 0.05){
      abline(fit, col=cols[sp], lwd=2)
    }
  }
}

dev.off()









pdf("/Users/au529793/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure4.pdf")

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(-t1$beta[t1$Variable=='Moisture'], 
     tapply(trait$PCA1, trait$Species, mean, na.rm=T),
     bg=cols, pch=21,
     xlim=c(-0.01,0.65),
     ylim=c(-6,6),
     cex=2,
     xlab="Moisture effect on Survival",
     ylab="PCA1")
legend('topright', legend=levels(surv$Species),
       cex=0.5, bty='n', pch=21, pt.bg = cols, pt.cex=1)
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
     xlab="Moisture effect on Growth",
     ylab="PCA1")
legend('topright', legend=levels(surv$Species),
       cex=0.5, bty='n', pch=21, pt.bg = cols, pt.cex=1)
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
     ylab="PCA2")
legend('topright', legend=levels(surv$Species),
       cex=0.5, bty='n', pch=21, pt.bg = cols, pt.cex=1)
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
     xlab="Moisture effect on Growth",
     ylab="PCA2")
legend('topright', legend=levels(surv$Species),
       cex=0.5, bty='n', pch=21, pt.bg = cols, pt.cex=1)
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







# pdf("/Users/au529793/Desktop/Figure8.v2.pdf")
# 
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# boxplot(trait$PCA1 ~ trait$Species, 
#         at=-t1$beta[t1$Variable=='Moisture'], 
#         col=cols, 
#         xlim=c(0,0.4),
#         ylim=c(-8,7),
#         boxwex=0.025, axes=F, 
#         xlab="Moisture effect on Survival",
#         ylab="PCA1")
# axis(1); axis(2); box()
# legend('topright', legend=levels(surv$Species),
#        cex=0.5, bty='n', pch=21, pt.bg = cols, pt.cex = 1.05)
# 
# boxplot(trait$PCA2 ~ trait$Species, 
#         at=-unlist(lapply(surv_fits, function(x) coefficients(x)[1])), col=cols, 
#         xlim=c(0,0.4),
#         ylim=c(-8,7),
#         boxwex=0.025, axes=F, 
#         xlab="Moisture effect on Survival",
#         ylab="PCA2")
# axis(1); axis(2); box()
# 
# boxplot(trait$PCA1 ~ trait$Species, 
#         at=unlist(lapply(grow_fits, function(x) fixef(x)[2])), col=cols, 
#         xlim=c(-0.01,0.05),
#         ylim=c(-8,7),
#         boxwex=0.005, axes=F, 
#         xlab="Moisture effect on Growth",
#         ylab="PCA1")
# axis(1); axis(2); box()
# 
# boxplot(trait$PCA2 ~ trait$Species, 
#         at=unlist(lapply(grow_fits, function(x) fixef(x)[2])), col=cols, 
#         xlim=c(-0.01,0.05),
#         ylim=c(-8,7),
#         boxwex=0.005, axes=F, 
#         xlab="Moisture effect on Growth",
#         ylab="PCA2")
# axis(1); axis(2); box()
# 
# dev.off()



