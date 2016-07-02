library(survival)

control_survival = data.frame(type = "Control", survival = c(14, 13, 13, 15, 13, 13))
niclo_survival =  data.frame(type = "Niclosamide",  survival =c(13, 14, 14, 13, 15, 14))
nen_survival =  data.frame(type = "NEN",  survival =c(16, 18, 15, 17, 17, 18))

sample_survival = rbind(control_survival, nen_survival) #niclo_survival, 

fit <- survival::coxph(Surv(survival) ~ type, data=sample_survival)
summary(fit)

pdf("~/Documents/stanford/hcc/data/LIHC/nen_survival.pdf")
  my.fit1 <- survfit(Surv(survival) ~ type, data = sample_survival) 
  plot(my.fit1, xlab="weeks", ylab = "Cumulative Survival",  col=c("green", "red"), cex.lab=1.5, cex.axis=1.5)
  text(3, .03, labels = "Log rank P=0.0009")
  legend(1, .3, c("Control", "NEN"), lty=1, col = c("green", "red"))
dev.off()

