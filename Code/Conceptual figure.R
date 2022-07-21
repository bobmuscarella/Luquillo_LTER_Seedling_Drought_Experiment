# Conceptual Figure

pdf("/Users/au529793/Desktop/Figure1-alt2.pdf", width = 7, height=5)

par(mar=c(5,4,1,1))

set.seed(42)
x <- rnorm(100)
y <- ifelse(x < 0, 0, 1)

y[x>0] <- sample(0:1, sum(x>0), replace=T, prob=c(0.1, 0.9))
y[x<0] <- sample(0:1, sum(x<0), replace=T, prob=c(0.9, 0.1))

plot(x, y, ylab="Probability of survival", 
     xlab="Soil moisture (%)", axes=F, pch=NA, xlim=c(-2.5, 2.5))
axis(1, label=seq(0,100,10), at=seq(-3,7))
axis(2, las=2, labels=seq(0,100,25), at=seq(0,1,0.25))
box()
mod <- glm(y ~ x, family='binomial')

nd <- data.frame(x=seq(-2.5, 2.5, 0.01))
ypred <- predict(mod, newdata = nd, type='response')

lines(nd$x, ypred, type='l', lwd=3)

arrows(nd[which(round(ypred, 2)==0.50)[1],1], 0.5, 
       nd[which(round(ypred, 2)==0.50)[1],1], 0, len=0.1, lwd=1.5)

segments(-3.5, 0.5, 
       nd[which(round(ypred, 2)==0.50)[1],1], 0.5, len=0.1, lwd=1.5, lty=2)

# points(nd[which(round(ypred, 2)==0.50)[1],1], 0.5, pch=21, cex=2, lwd=2)

text(0.65, 0.8, "β Demographic Rate", #"Sensitivity",
     col=2, pos=4, font=2)

text(0.65, 0.65, "Slope of relationship between 
a given demographic rate 
and soil moisture.", cex=0.8, col=2, pos=4)


text(0.1, 0.25, "Sur50", #"Tolerance", 
     pos=4, font=2, col=4)

text(0.1, 0.1, "100 - Soil moisture 
that corresponds to a 
50% survival probability.", cex=0.8, col=4, pos=4)

text(-1.5, 0.53, "50% Predicted Survival", 
     col=4, font=3, cex=0.7)

dev.off()