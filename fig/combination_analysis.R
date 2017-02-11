require(xlsx)
library(plyr)
library(ggplot2)

combination = read.xlsx("fig/fig_data.xlsx", sheetName = "combination")

cdata <- ddply(combination, .(drug), summarise, 
               N    = length(volume),
               mean = mean(volume),
               sd   = sd(volume, na.rm =T),
               se   = sd(volume, na.rm =T) / sqrt(length(volume))
            )

ylim_max = (max(cdata$mean) * 2)

cdata$star <- ""
drugs = c("NEN", "sorafenib", "NEN & sorafenib")
for (drug in drugs){
  p = t.test(combination$volume[combination$drug == "drinking water"], combination$volume[combination$drug == drug], "two.sided")$p.value
  if (p < 0.05){
    cdata$star[cdata$drug == drug]  <- "*"
  }else if (p < 0.01){
    cdata$star[cdata$drug == drug]  <- "**"
  }else if (p < 0.001){
    cdata$star[cdata$drug == drug]  <- "***"
  }
}

cdata$drug[cdata$drug == "drinking water"] = "Vehicle"

pdf(paste("fig/", "combination", ".pdf", sep=""))
q =ggplot(cdata, aes(x=drug, y=mean)) + coord_cartesian(ylim=c(0,ylim_max))  +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  geom_text(aes(label=star, y= mean + se + .01), colour="black", vjust= 0, size=10) +
  xlab("") +
  ylab("RQ") +
  ggtitle(gene) 

print(q)
dev.off()
