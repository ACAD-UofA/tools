# mammoth data from Pecnerova et al. 2017

library(dplyr)
library(reshape2)

source('explore.R')

d <- read.table("data/mammoths.csv", sep="\t", header=T)
d$material <- addNA(d$material)
d$material2 <- addNA(d$material2)

mres <- tibble(data=character(), model=character(),`coeff p`=numeric(), `LRT`=numeric())
i <- 1

# are sex ratios significantly different?
mf <- table(d$sex)
mf
mf["M"]/sum(mf)
binom.test(mf["M"], mf["M"]+mf["F"], alternative="greater")
m.null <- glm(sex~1, data=d, family=binomial)
summary(m.null)
mres[i,] <- list("mammoths", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

#explore.categorical(d, "material", "sex")
#mres[i,] <- list("mammoths", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
#i <- i + 1
m <- explore.categorical(d, "material2", "sex")
mres[i,] <- list("mammoths", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "loc", "sex")
mres[i,] <- list("mammoths", "loc", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "loc2", "sex")
mres[i,] <- list("mammoths", "loc2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "age", "sex")
mres[i,] <- list("mammoths", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

mres.mammoths <- mres
summary.mammoths <- dcast(mres, model~data, value.var='LRT')
summary.mammoths
