rdata <- read.csv("D:/JHU/lesson/BLAM/code/rdata.csv")
rdata

#fit repeated measures ANOVA model
sigma_p <- aov(Sigma~factor(MR)*factor(SWM)*factor(Bi)+Error(factor(Subject)), data = rdata)
alpha_p <- aov(Alpha~factor(MR)*factor(SWM)*factor(Bi)+Error(factor(Subject)), data = rdata)

#view model summary
summary(sigma_p)
summary(alpha_p)