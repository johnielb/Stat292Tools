library(car)
library(ggplot2)
library(tidyr)

rimu <- tibble(c(6.5,8.8,5.3,8.2,9.2,6.6,6.3,8.0,5.8,10.9,8.7,6.1,12.3,5.8,4.4,2.6,3.6,3.7),
               as.factor(c(1,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4)))
colnames(rimu) <- c("diam_breast_height", "site")
rimu$logDBH <- log(rimu$diam_breast_height)

ggplot(rimu, aes(x=site, y=diam_breast_height)) + 
  geom_boxplot() +
  stat_summary(geom="point", fun=mean, size=3, shape=1)

cat("ANOVA summary:\n")
anova <- aov(diam_breast_height ~ site, data = rimu)
summary(anova)

leveneTest(diam_breast_height ~ site, data = rimu)

kruskal.test(diam_breast_height ~ site, data = rimu)

TukeyHSD(anova)

plot(anova, 1)
plot(anova, 2)

cat("log ANOVA summary:\n")
log_anova <- aov(logDBH ~ site, data = rimu)
summary(log_anova)

leveneTest(logDBH ~ site, data = rimu)

kruskal.test(logDBH ~ site, data = rimu)

TukeyHSD(log_anova)

plot(log_anova, 1)
plot(log_anova, 2)