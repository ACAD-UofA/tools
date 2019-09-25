# mammoth data from Pecnerova et al. 2017

library(dplyr)
library(reshape2)

source('explore.R')

d <- read.table("data/mammoths.csv", sep="\t", header=T)
d$material <- addNA(d$material)
d$material2 <- addNA(d$material2)

mres <- tibble(`data`=character(), `model`=character(), `coeff`=character(),
                `Estimate`=numeric(),
                `Std. Error`=numeric(),
                `z value`=numeric(),
                `Pr(>|z|)`=numeric())

lres <- tibble(`data`=character(), `model`=character(),
                `Df`=numeric(),
                `Deviance`=numeric(),
                `Resid. Df`=numeric(),
                `Resid. Dev`=numeric(),
                `Pr(>Chi)`=numeric())

mappend <- function(m1, m2, data, model) {
  m <- as_tibble(coefficients(summary(m2)), rownames="coeff")
  if (model %in% colnames(m2$data) && is.factor(m2$data[model][,])) {
    fn <- function(x) gsub(paste("^", model, sep=""), "", x)
    m["coeff"] <- lapply(m["coeff"], fn)
  }
  rbind(m1, cbind(data, model, m))
}

lrtappend <- function(l, m2, data, model) {
  m <- anova(m2,test="LRT")[2,]
  row.names(m) = NULL
  rbind(l, cbind(data, model, m))
}

# are sex ratios significantly different?
mf <- table(d$sex)
mf
mf["M"]/sum(mf)
binom.test(mf["M"], mf["M"]+mf["F"], alternative="greater")
m.null <- glm(sex~1, data=d, family=binomial)
summary(m.null)
mres <- mappend(mres, m.null, "mammoths", "intercept-only")

m <- explore.categorical(d, "material", "sex")
mres <- mappend(mres, m, "mammoths", "material")
lres <- lrtappend(lres, m, "mammoths", "material")

m <- explore.categorical(d, "loc", "sex")
mres <- mappend(mres, m, "mammoths", "loc")
lres <- lrtappend(lres, m, "mammoths", "loc")

m <- explore.categorical(d, "loc2", "sex")
mres <- mappend(mres, m, "mammoths", "loc2")
lres <- lrtappend(lres, m, "mammoths", "loc2")

m <- explore.continuous(d, "age", "sex")
mres <- mappend(mres, m, "mammoths", "age")
lres <- lrtappend(lres, m, "mammoths", "age")

mres.mammoths <- mres
lres.mammoths <- lres
