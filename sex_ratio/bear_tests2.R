
library(dplyr)
library(reshape2)
library(car)

source('explore.R')
source('mapstuff.R')
source('two-sample-tests.R')

d <- read.table("data/brownbears.csv", sep="\t", header=TRUE)
d$site_type <- factor(d$site_type, levels=rev(levels(d$site_type)))

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
mres <- mappend(mres, m.null, "brown bears", "intercept-only")

m <- explore.categorical(d, "site_type", "sex")
mres <- mappend(mres, m, "brown bears", "site_type")
lres <- lrtappend(lres, m, "brown bears", "site_type")

m <- explore.categorical(d, "material", "sex")
mres <- mappend(mres, m, "brown bears", "material")
lres <- lrtappend(lres, m, "brown bears", "material")

m <- explore.categorical(d, "material2", "sex")
mres <- mappend(mres, m, "brown bears", "material2")
lres <- lrtappend(lres, m, "brown bears", "material2")


m <- explore.continuous(d, "lat", "sex")
mres <- mappend(mres, m, "brown bears", "lat")
lres <- lrtappend(lres, m, "brown bears", "lat")

m <- explore.continuous(d, "lon", "sex")
mres <- mappend(mres, m, "brown bears", "lon")
lres <- lrtappend(lres, m, "brown bears", "lon")

m <- explore.continuous(d, "alt", "sex") # *
mres <- mappend(mres, m, "brown bears", "alt")
lres <- lrtappend(lres, m, "brown bears", "alt")

m <- explore.continuous(d, "age", "sex") # *
mres <- mappend(mres, m, "brown bears", "age")
lres <- lrtappend(lres, m, "brown bears", "age")

m <- explore.continuous(d, "endog", "sex")
mres <- mappend(mres, m, "brown bears", "endog")
lres <- lrtappend(lres, m, "brown bears", "endog")

m <- explore.continuous(d, "GC", "sex")
mres <- mappend(mres, m, "brown bears", "GC")
lres <- lrtappend(lres, m, "brown bears", "GC")

m <- explore.continuous(d, "readlen", "sex")
mres <- mappend(mres, m, "brown bears", "readlen")
lres <- lrtappend(lres, m, "brown bears", "readlen")

m <- explore.continuous(d, "C2T", "sex")
mres <- mappend(mres, m, "brown bears", "C2T")
lres <- lrtappend(lres, m, "brown bears", "C2T")


m <- explore.categorical(d, "alpine", "sex") # *
mres <- mappend(mres, m, "brown bears", "alpine")
lres <- lrtappend(lres, m, "brown bears", "alpine")

d2 <- subset(d, alpine!="Alps")
d2 <- droplevels(d2)
d3 <- subset(d, alpine=="Alps")
d3 <- droplevels(d3)

m.null <- glm(sex~1, data=d2, family=binomial)
summary(m.null)
mres <- mappend(mres, m.null, "brown bears (non Alps)", "intercept-only")

m <- explore.categorical(d2, "material", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "material")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "material")

m <- explore.categorical(d2, "material2", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "material2")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "material2")

m <- explore.categorical(d2, "site_type", "sex") # *
mres <- mappend(mres, m, "brown bears (non Alps)", "site_type")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "site_type")

m <- explore.continuous(d2, "lat", "sex") # *
mres <- mappend(mres, m, "brown bears (non Alps)", "lat")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "lat")

m <- explore.continuous(d2, "lon", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "lon")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "lon")

m <- explore.continuous(d2, "alt", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "alt")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "alt")

m <- explore.continuous(d2, "age", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "age")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "age")

m <- explore.continuous(d2, "endog", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "endog")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "endog")

m <- explore.continuous(d2, "GC", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "GC")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "GC")

m <- explore.continuous(d2, "readlen", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "readlen")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "readlen")

m <- explore.continuous(d2, "C2T", "sex")
mres <- mappend(mres, m, "brown bears (non Alps)", "C2T")
lres <- lrtappend(lres, m, "brown bears (non Alps)", "C2T")


m.null <- glm(sex~1, data=d3, family=binomial)
summary(m.null)
mres <- mappend(mres, m.null, "brown bears (Alps)", "intercept-only")

m <- explore.categorical(d3, "material", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "material")
lres <- lrtappend(lres, m, "brown bears (Alps)", "material")

m <- explore.categorical(d3, "material2", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "material2")
lres <- lrtappend(lres, m, "brown bears (Alps)", "material2")

m <- explore.continuous(d3, "alt", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "alt")
lres <- lrtappend(lres, m, "brown bears (Alps)", "alt")

m <- explore.continuous(d3, "age", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "age")
lres <- lrtappend(lres, m, "brown bears (Alps)", "age")

m <- explore.continuous(d3, "endog", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "endog")
lres <- lrtappend(lres, m, "brown bears (Alps)", "endog")

m <- explore.continuous(d3, "lat", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "lat")
lres <- lrtappend(lres, m, "brown bears (Alps)", "lat")

m <- explore.continuous(d3, "lon", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "lon")
lres <- lrtappend(lres, m, "brown bears (Alps)", "lon")

m <- explore.continuous(d3, "GC", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "GC")
lres <- lrtappend(lres, m, "brown bears (Alps)", "GC")

m <- explore.continuous(d3, "readlen", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "readlen")
lres <- lrtappend(lres, m, "brown bears (Alps)", "readlen")

m <- explore.continuous(d3, "C2T", "sex")
mres <- mappend(mres, m, "brown bears (Alps)", "C2T")
lres <- lrtappend(lres, m, "brown bears (Alps)", "C2T")

mres.bears <- mres
lres.bears <- lres
