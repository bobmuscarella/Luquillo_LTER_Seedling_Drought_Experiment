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

