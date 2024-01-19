# Conceptual Figure revised

pdf("~/Projects/GIT/Luquillo_LTER_Seedling_Drought_Experiment/Figures/Figure1-20231130.pdf", width=8)

library(pBrackets)

par(mar=c(5,4,1,1))

sp1x <- seq(0, 60, by = 0.1)
sp1y <- plogis(seq(-5,5,length.out=601))

sp2x <- c(seq(0, 20, by = 0.1), 60)
sp2y <- c(plogis(seq(-5,5,length.out=201)),1)


plot(sp1x, sp1y, type='l', col=4, lwd=4,
     xlab="Soil Moisture (%)", ylab="Survival probability", xlim=c(3,50))
points(sp2x, sp2y, xlab = 'probabilities', type='l', col=2, lwd=4)

abline(v=moist_percentiles, lty=3)

text(mean(moist_percentiles), 0.01, 'Experimental soil moisure range
(10% and 90% percentiles)', cex=0.9)

polygon(x=c(-10, moist_percentiles[1], moist_percentiles[1], -10),
        y=c(-1,-1,2,2), col=rgb(0,0,0,0.2), border=NA)

polygon(x=c(100, moist_percentiles[2], moist_percentiles[2], 100),
        y=c(-1,-1,2,2), col=rgb(0,0,0,0.2), border=NA)

y1pt1 <- sp1y[which(round(sp1x,1) == round(moist_percentiles[1],1))]
y2pt1 <- sp2y[which(round(sp2x,1) == round(moist_percentiles[1],1))]

y1pt2 <- sp1y[which(round(sp1x,1) == round(moist_percentiles[2],1))]
y2pt2 <- sp2y[which(round(sp2x,1) == round(moist_percentiles[2],1))]

points(rep(moist_percentiles[1], 2), c(y1pt1, y2pt1), bg=c(4,2), pch=21, cex=2)
points(rep(moist_percentiles[2], 2), c(y1pt2, 1), bg=c(4,2), pch=21, cex=2)

text(18, 0.6, "Tolerance")

arrows(18, 0.58, moist_percentiles[1]+0.75, y1pt1+0.02,
       len=0.1, lwd=1.5, col=4)
arrows(18, 0.62, moist_percentiles[1]+0.75, y2pt1-0.02,
       len=0.1, lwd=1.5, col=2)

segments(moist_percentiles[1], y2pt1, moist_percentiles[2], y2pt1, col=2, lty=3)
segments(moist_percentiles[1], y1pt1, moist_percentiles[2], y1pt1, col=4, lty=3)

brackets(38.8, y2pt1+0.005, 38.8, 1-0.015, h=1, col=2, lwd=1.5)
brackets(38.8, y1pt1+0.005, 38.8, y1pt2-0.015, h=1, col=4, lwd=1.5)

# arrows(38.8, y2pt1+0.005, 38.8, 1-0.015, code=3, col=2, lwd=1.5, len=0.1)
# arrows(38.8, y1pt1+0.005, 38.8, y1pt2-0.015, code=3, col=4, lwd=1.5, len=0.1)

text(34.75, 0.44, "Sensitivity", col=4)
text(34.75, 0.925, "Sensitivity", col=2)

arrows(34.5, 0.01, moist_percentiles[2], len=0.075, lwd=1.5)
arrows(18.5, 0.01, moist_percentiles[1], len=0.075, lwd=1.5)

dev.off()
