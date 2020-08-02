### Master trait dataset

### read in the dataset
library(readxl)
LUQDrought_RootTraitData_2020_06_04 <- read_excel("C:/Users/roell/Dropbox/LUQ Drought EXP EcoPhys - Matlaga/Datasets/Trait Data_ROEL/ROOT TRAITS/LUQDrought_RootTraitData_2020-06-04.xlsx", 
                                                  col_types = c("numeric", "numeric", "text", 
                                                                "numeric", "numeric", "text", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric"))
Traitdata <- LUQDrought_RootTraitData_2020_06_04 

### Boxplots of leaf an stem traits

### Color coding
library(RColorBrewer)
library(scales)
Traitdata$col <- brewer.pal(8, "Dark2")[as.factor(Traitdata$Species)]

### Manual color coding
Color1 <- c("#1B9E77","#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
Color2 <- c("#D95F02","#7570B3", "#E7298A", "#E6AB02", "#A6761D")
Color3 <- c("#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
Color4 <- c("#D95F02","#7570B3", "#E7298A", "#E6AB02", "#A6761D", "#666666")

### Compute mass fractions
Traitdata$`Stem dry mass (g)`[Traitdata$Species=="Premon"] <- 0

Traitdata$drymass <-rowSums(cbind(Traitdata$`RootDryMass(g)`, Traitdata$`Stem dry mass (g)`, Traitdata$`Leaf dry mass (g)`, na.rm=TRUE))

Traitdata$leaf_fraction <- (Traitdata$`Leaf dry mass (g)`/Traitdata$drymass)
Traitdata$stem_fraction <- (Traitdata$`Stem dry mass (g)`/Traitdata$drymass)
Traitdata$root_fraction <- (Traitdata$`RootDryMass(g)`/Traitdata$drymass)

### compute mean leaf thickness value
Traitdata$Leafthickness <- (Traitdata$`Leaf thick 1 (mm)`+ Traitdata$`Leaf thick 2 (mm)`+ Traitdata$`Leaf thick 3 (mm)`)/3

### Prepare data to take certain species out, if needed
Tetbal <- Traitdata[ which(Traitdata$Species=='Tetbal'), ]
Manbid <- Traitdata[ which(Traitdata$Species=='Manbid'), ]
Guagui <- Traitdata[ which(Traitdata$Species=='Guagui'), ]
Inlau <- Traitdata[ which(Traitdata$Species=='Inlau'), ]
Schmor <- Traitdata[ which(Traitdata$Species=='Schmor'), ]
Premon <- Traitdata[ which(Traitdata$Species=='Premon'), ]
Urebac <- Traitdata[ which(Traitdata$Species=='Urebac'), ]
Cecsch <- Traitdata[ which(Traitdata$Species=='Cecsch'), ]

Traitdata1 <- rbind(Tetbal, Manbid, Guagui, Inlau, Schmor, Premon)
Traitdata2 <- rbind(Tetbal, Manbid, Guagui, Inlau, Schmor)
Traitdata3 <- rbind (Tetbal, Manbid, Guagui, Inlau, Schmor, Premon, Urebac)
Traitdata4 <- rbind(Tetbal, Manbid, Guagui, Inlau, Schmor, Urebac)

## Boxplots leaf and stem traits

library(ggplot2)

### take CECSCH away in all trait boxplots
### So use traitdata3 and Color3

### LMA triat
LMA <- ggplot(Traitdata3, aes(fill = Species , y=`LMA (g/m2)`, x=Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "LMA (g/m2)")+ theme(text = element_text(size=20))
LMA

### LA trait
LA <- ggplot(Traitdata3, aes(fill = Species , y=`Leaf area (cm2)`, x=Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "LA (cm2)")+ theme(text = element_text(size=20))
LA

### Leaf thickness trait
LT <- ggplot(Traitdata3, aes(fill = Species , y=Leafthickness, x=Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "Leaf thickness (mm)")+ theme(text = element_text(size=20)) +  expand_limits(y = 0)
LT

### LDMC trait
LDMC <- ggplot(Traitdata3, aes(fill = Species , y=`LDMC (%)`, x=Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "LDMC (%)")+ theme(text = element_text(size=20))
LDMC

### WD trait
Traitdata4$`Stem density (g/cm3)` <- as.numeric(Traitdata4$`Stem density (g/cm3)`)
WD <- ggplot(Traitdata4, aes(fill = Species , y=`Stem density (g/cm3)` , x=Species)) + geom_boxplot() + scale_fill_manual(values = Color4) + labs(x = "" , y = "Stem density (g/cm3")+ theme(text = element_text(size=20))
WD

### WDMC trait
## take premon out, has no stem in seedling stage
## Urebac data varies so much to a low amount of datapoints, can be taken out as well
WDMC <- ggplot(Traitdata2, aes(fill = Species , y=`WDMC (%)`, x=Species)) + geom_boxplot() + scale_fill_manual(values = Color2) + labs(x = "" , y = "WDMC (%)")+ theme(text = element_text(size=20))
WDMC


### arrange muliple plots in one grid
# Link: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

library(ggpubr)
ggarrange(LMA, LA, LT, LDMC, WDMC, WD, common.legend = TRUE, legend = "bottom")

### Lineair regression leaf and stem traits

# Convert species names to a factor
Traitdata$Species <- as.factor(Traitdata$Species)

### plot traits
cols <- brewer.pal(8, "Dark2")

par(mfcol=c(2,3), mar=c(4,4,1,1))
for(k in 1:6){
  yvar <- c("LMA (g/m2)", "LDMC (%)", "Leaf area (cm2)", "Leafthickness","Stem density (g/cm3)", "WDMC (%)")[k]
  y <- unlist(Traitdata[,yvar])
  x <- Traitdata$Moisture
  plot(x, y, bg=cols[as.factor(Traitdata$Species)], pch=21, 
       xlab="Soil_moisture %",
       ylab=c("LMA (g/m2)",
              "LDMC (%)",
              "Leaf area (cm2)",
              "Leafthickness",
              "Stem density (g/cm3)",
              "WDMC (%)")[k])
  for(i in seq_along(levels(Traitdata$Species))){
    tmpy <- y[Traitdata$Species == levels(Traitdata$Species)[i]]
    tmpx <- x[Traitdata$Species == levels(Traitdata$Species)[i]]
    if(sum(!is.na(tmpy))>2){
      fit <- lm(tmpy ~ tmpx)
      lty <- ifelse(summary(fit)$coefficients[2,4] < 0.05, 1, 2)
      lwd <- ifelse(summary(fit)$coefficients[2,4] < 0.05, 3, 1)
      abline(fit, col=cols[i], lwd=lwd, lty=lty)
    }
    if(k==1){ 
      legend('topleft', legend=levels(Traitdata$Species),
             pt.bg=cols[seq_along(levels(Traitdata$Species))],
             pch=21, cex=0.9, pt.cex=2.5, bg='white')}
  }
}


### Boxplots of root traits

### Root mass fraction
RMF <- ggplot(Traitdata3, aes(fill = Species , y= root_fraction, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "Root mass fraction (%)")+ theme(text = element_text(size=20))
RMF

### Rooting depth
Rdepth <- ggplot(Traitdata3, aes(fill = Species , y=`Root depth(cm)`, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "Root depth(cm)")+ theme(text = element_text(size=20))
Rdepth

### Root length
Rlength <- ggplot(Traitdata3, aes(fill = Species , y= `Length(cm)`, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "Root length(cm)")+ theme(text = element_text(size=20))
Rlength

### Average Diameter
Averagediameter<- ggplot(Traitdata3, aes(fill = Species , y=`AvgDiam(mm)`, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "AvgDiam(mm)")+ theme(text = element_text(size=20)) +  expand_limits(y = 0)
Averagediameter

### Tips
Tips <- ggplot(Traitdata3, aes(fill = Species , y= Tips, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "Tips")+ theme(text = element_text(size=20))
Tips

### SRL
SRL <- ggplot(Traitdata3, aes(fill = Species , y=`SRL(cm/g)`, x = Species)) + geom_boxplot() + scale_fill_manual(values = Color3) + labs(x = "" , y = "SRL(cm/g)")+ theme(text = element_text(size=20))
SRL

ggarrange(RMF, Rdepth, Rlength, Averagediameter, Tips, SRL, common.legend = TRUE, legend = "bottom")

### Lineair regression root traits

par(mfcol=c(2,3), mar=c(4,4,1,1))
for(k in 1:6){
  yvar <- c("root_fraction", "Root depth(cm)", "Length(cm)","AvgDiam(mm)", "Tips", "SRL(cm/g)")[k]
  y <- unlist(Traitdata[,yvar])
  x <- Traitdata$Moisture
  plot(x, y, bg=cols[as.factor(Traitdata$Species)], pch=21, 
       xlab="Soil_moisture %",
       ylab=c("Root mass fration(%)",
              "Root depth(cm)",
              "Root length (cm)",
              "Average Diameter(mm)",
              "Tips",
              "SRL(cm/g)")[k])
  for(i in seq_along(levels(Traitdata$Species))){
    tmpy <- y[Traitdata$Species == levels(Traitdata$Species)[i]]
    tmpx <- x[Traitdata$Species == levels(Traitdata$Species)[i]]
    if(sum(!is.na(tmpy))>2){
      fit <- lm(tmpy ~ tmpx)
      lty <- ifelse(summary(fit)$coefficients[2,4] < 0.05, 1, 2)
      lwd <- ifelse(summary(fit)$coefficients[2,4] < 0.05, 3, 1)
      abline(fit, col=cols[i], lwd=lwd, lty=lty)
    }
    if(k==1){
      legend('topleft', legend=levels(Traitdata$Species),
             pt.bg=cols[seq_along(levels(Traitdata$Species))],
             pch=21, cex=0.75, pt.cex=1.5, bg='white')}
  }
}


### PCA with leaf and root traits
### remove NA values and 2 urebac outliers (56,94)
PCA_data <- Traitdata[-c(56,72,94,194,195,324,339,401,497),]

### PCA with only leaf traits
library(vegan)

Leafarea <- PCA_data$`Leaf area (cm2)`
LMA <- PCA_data$`LMA (g/m2)`
LDMC <- PCA_data$`LDMC (%)`
Leafthickness <- PCA_data$Leafthickness
Rootmassfraction <- PCA_data$root_fraction
Rootdepth <- PCA_data$`Root depth(cm)`
Rootlength <- PCA_data$`Length(cm)`
AvgRootDiam <- PCA_data$`AvgDiam(mm)`
Tips <- PCA_data$Tips
SRL <- PCA_data$`SRL(cm/g)`

ord <- prcomp( ~ Leafarea  +
                 LMA+
                 LDMC +
                 Leafthickness +
                 Rootmassfraction +
                 Rootlength +
                 Rootdepth +
                 AvgRootDiam+
                 Tips+
                 SRL ,
               center = TRUE, scale = TRUE)
summary(ord)

library("factoextra")
### visualize the variables used in the PCA
A1 <- fviz_pca_var(ord, col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE # Avoid text overlapping)
)

###  visualize results for individuals
B1 <- fviz_pca_ind(ord,
                   label = "none", # hide individual labels
                   habillage = PCA_data$Species, # color by groups
                   palette = c("#1B9E77","#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
                   addEllipses = TRUE # Concentration ellipses
)


### Put the plots together
PCAplot1 <- ggarrange(A1, B1, ncol = 2, nrow = 1)
PCAplot1

### Contribution of variables
tres.var <- get_pca_var(ord)
contrib1 <- tres.var$contrib
contrib1

### Coordinates for individuals and plotting it against soil moisture
res.ind <- get_pca_ind(ord)
Coord <- res.ind$coord

Dim1 <- Coord[,-c(2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
Dim2 <- Coord[,-c(1,3,4,5,6,7,8,9,10,11,12,13,14,15)]

PCA_data$Dim1 <- Dim1
PCA_data$Dim2 <- Dim2

Dimension1 = lm(Dim1 ~ Moisture, data = PCA_data)
summary(Dimension1)
Dimension2 = lm(Dim2 ~ Moisture, data = PCA_data)
summary(Dimension2)

ggplot(data = PCA_data, aes(x = Moisture, y = Dim1))+ 
  geom_point(color='black') +
  geom_smooth(method = "lm",color='red', se = FALSE)+
  labs(y= " PCA dimension 1", x = "Soil Moisture")


ggplot(data = PCA_data, aes(x = Moisture, y = Dim2))+ 
  geom_point(color='black') +
  geom_smooth(method = "lm",color='red', se = FALSE)+
  labs(y= "PCA dimension 2", x = "Soil Moisture")
