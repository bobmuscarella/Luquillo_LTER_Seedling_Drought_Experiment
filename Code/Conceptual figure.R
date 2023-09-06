# Conceptual Figure

pdf("~/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure1-20230832.pdf", width = 7, height=5)

library(pBrackets)

# Conceptual Figure (2023-08-31)

# pdf("/Users/au529793/Desktop/Figure1-alt3-20230832.pdf", width = 7, height=5)

par(mar=c(5,4,1,1))

plot(seq(0, 50, by = 0.1), plogis(seq(-5,5,length.out=501)), type='l', col=4, lwd=4,
     xlab="Soil Moisture (%)", ylab="Probability of Survival")
points(c(seq(0, 20, by = 0.1), 50), c(plogis(seq(-5,5,length.out=201)),1), xlab = 'probabilities', type='l', col=2, lwd=4)
# points(c(seq(0, 100, by = 0.1)), c(plogis(seq(-5,5,length.out=1001))), xlab = 'probabilities', type='l', col=7, lwd=4)

abline(v=moist_percentiles, lty=3)
arrows(moist_percentiles[1], 0.05, moist_percentiles[2], 0.05, 
       len=0.1, code=3, lwd=1.5)
text(mean(moist_percentiles), 0.01, 'Experimental Soil Moisure
(10% and 90% percentiles)', cex=0.75)

points(rep(moist_percentiles[1], 2), c(0.091, 0.855), bg=c(4,2), pch=21, cex=2)
points(rep(moist_percentiles[2], 2), c(0.945, 1), bg=c(4,2), pch=23, cex=2)

text(5, 0.8, "Tolerance")

arrows(8.25, 0.79, 13, 0.12,
       len=0.1, lwd=1.5, col=4)

arrows(8.25, 0.79, 12.8, 0.85,
       len=0.1, lwd=1.5, col=2)


text(45.2, mean(c(0.945, 0.091)), "Sensitivity", col=4)
text(45.2, mean(c(1, 0.855)), "Sensitivity", col=2)

brackets(40, 0.945, 40, 0.091, h=2, col=4, lwd=2)
brackets(40, 1, 40, 0.855, h=2, col=2, lwd=2)
segments(moist_percentiles[1], 0.855, 40, 0.855, col=2, lty=3)
segments(moist_percentiles[1], 0.091, 40, 0.091, col=4, lty=3)

dev.off()



# 
# 
# par(mar=c(5,4,1,1))
# 
# set.seed(42)
# x <- rnorm(100)
# y <- ifelse(x < 0, 0, 1)
# 
# y[x>0] <- sample(0:1, sum(x>0), replace=T, prob=c(0.1, 0.9))
# y[x<0] <- sample(0:1, sum(x<0), replace=T, prob=c(0.9, 0.1))
# 
# plot(x, y, ylab="Probability of survival", 
#      xlab="Soil moisture (%)", axes=F, pch=NA, xlim=c(-2.5, 2.5))
# axis(1, label=seq(0,100,10), at=seq(-3,7))
# axis(2, las=2, labels=seq(0,100,25), at=seq(0,1,0.25))
# box()
# mod <- glm(y ~ x, family='binomial')
# 
# nd <- data.frame(x=seq(-2.5, 2.5, 0.01))
# ypred <- predict(mod, newdata = nd, type='response')
# 
# lines(nd$x, ypred, type='l', lwd=3)
# 
# # arrows(nd[which(round(ypred, 2)==0.50)[1],1], 0.5, 
# #        nd[which(round(ypred, 2)==0.50)[1],1], 0, len=0.1, lwd=1.5)
# 
# # segments(-3.5, 0.5, 
# #        nd[which(round(ypred, 2)==0.50)[1],1], 0.5, len=0.1, lwd=1.5, lty=2)
# 
# abline(v=c(1.35175, 39.1400), lwd=1.5, lty=2)
# 
# 
# # points(nd[which(round(ypred, 2)==0.50)[1],1], 0.5, pch=21, cex=2, lwd=2)
# 
# text(0.65, 0.8, "β Demographic Rate", #"Sensitivity",
#      col=2, pos=4, font=2)
# 
# text(0.65, 0.65, "Slope of relationship between 
# a given demographic rate 
# and soil moisture.", cex=0.8, col=2, pos=4)
# 
# 
# text(0.1, 0.25, "Sur50", #"Tolerance", 
#      pos=4, font=2, col=4)
# 
# text(0.1, 0.1, "100 - Soil moisture 
# that corresponds to a 
# 50% survival probability.", cex=0.8, col=4, pos=4)
# 
# text(-1.5, 0.53, "50% Predicted Survival", 
#      col=4, font=3, cex=0.7)
# 
# dev.off()