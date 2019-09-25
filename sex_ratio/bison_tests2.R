
library(dplyr)
library(reshape2)

source('explore.R')
source('mapstuff.R')
source('two-sample-tests.R')

d = read.table("data/bison.csv", sep="\t", header=TRUE)

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

d$material = addNA(d$material)
d$material2 = addNA(d$material2)

# reorder factors to make plots more aesthetically pleasing
d$site_type <- factor(d$site_type, levels=rev(levels(d$site_type)))

# divide into European and American, at 60 degrees
d$loncat = cut(d$lon, breaks=c(-20,60,Inf), labels=c("EU","AM"))

# are sex ratios significantly different?
mf = table(d$sex)
mf
mf["M"]/sum(mf)
binom.test(mf["M"], mf["M"]+mf["F"], alternative="greater")
m.null = glm(sex~1, data=d, family=binomial)
summary(m.null)
mres <- mappend(mres, m.null, "bison", "intercept-only")

m <- explore.categorical(d, "loncat", "sex")
mres <- mappend(mres, m, "bison", "loncat")
lres <- lrtappend(lres, m, "bison", "loncat")

m <- explore.categorical(d, "site_type", "sex") # *
mres <- mappend(mres, m, "bison", "site_type")
lres <- lrtappend(lres, m, "bison", "site_type")

m <- explore.categorical(d, "material", "sex")
mres <- mappend(mres, m, "bison", "material")
lres <- lrtappend(lres, m, "bison", "material")

m <- explore.categorical(d, "material2", "sex")
mres <- mappend(mres, m, "bison", "material2")
lres <- lrtappend(lres, m, "bison", "material2")

# check that the 'reference' level isn't fucking us
d.sc2.ref.post = transform(d, material2=addNA(relevel(material2, ref="Postcrania")))
explore.categorical(d.sc2.ref.post, "material2", "sex")

m <- explore.continuous(d, "lat", "sex")
mres <- mappend(mres, m, "bison", "lat")
lres <- lrtappend(lres, m, "bison", "lat")

m <- explore.continuous(d, "lon", "sex")
mres <- mappend(mres, m, "bison", "lon")
lres <- lrtappend(lres, m, "bison", "lon")

m <- explore.continuous(d, "alt", "sex")
mres <- mappend(mres, m, "bison", "alt")
lres <- lrtappend(lres, m, "bison", "alt")

m <- explore.continuous(d, "age", "sex")
mres <- mappend(mres, m, "bison", "age")
lres <- lrtappend(lres, m, "bison", "age")

m <- explore.continuous(d, "endog", "sex")
mres <- mappend(mres, m, "bison", "endog")
lres <- lrtappend(lres, m, "bison", "endog")

m <- explore.continuous(d, "GC", "sex")
mres <- mappend(mres, m, "bison", "GC")
lres <- lrtappend(lres, m, "bison", "GC")

m <- explore.continuous(d, "readlen", "sex")
mres <- mappend(mres, m, "bison", "readlen")
lres <- lrtappend(lres, m, "bison", "readlen")

m <- explore.continuous(d, "C2T", "sex")
mres <- mappend(mres, m, "bison", "C2T")
lres <- lrtappend(lres, m, "bison", "C2T")



# test for differences in spatial distribution between sexes
#d.sex.sorted = subset(d[order(d$sex),], !is.na(lat) & !is.na(lon))
#gcm = greatcircle.dmatrix(d.sex.sorted$lat, d.sex.sorted$lon)
#kern1 = kernel.test(gcm, table(d.sex.sorted$sex))
#kern1
#hist(kern1$permutations,100)
#abline(v=kern1$statistic, col="red")

# evaluate collinearity of indep. vars.

# the contingency table is too sparse to check for any confounders of
# site type, except site_type=="Sediment"

explore.categorical(d, "material", "site_type")
explore.categorical(d, "material2", "site_type")
explore.categorical(d.sc2.ref.post, "material2", "site_type")

# compare 2 predictor model to simple model
m = glm(sex~site_type+material, family=binomial, data=d)
m.alt1 = update(m, sex~site_type, data=m$model)
summary(m)
anova(m.alt1,m,test="LRT")

m = glm(sex~site_type+material2, family=binomial, data=d)
m.alt1 = update(m, sex~site_type, data=m$model)
summary(m)
anova(m.alt1,m,test="LRT")

site.sediment <- subset(d, site_type=="Sediment")
site.sediment <- droplevels(site.sediment)

m.null = glm(sex~1, data=site.sediment, family=binomial)
mres <- mappend(mres, m.null, "bison (sediment)", "intercept-only")

m <- explore.categorical(site.sediment, "material", "sex")
mres <- mappend(mres, m, "bison (sediment)", "material")
lres <- lrtappend(lres, m, "bison (sediment)", "material")

m <- explore.categorical(site.sediment, "material2", "sex")
mres <- mappend(mres, m, "bison (sediment)", "material2")
lres <- lrtappend(lres, m, "bison (sediment)", "material2")

m <- explore.categorical(site.sediment, "loncat", "sex")
mres <- mappend(mres, m, "bison (sediment)", "loncat")
lres <- lrtappend(lres, m, "bison (sediment)", "loncat")

m <- explore.continuous(site.sediment, "lat", "sex")
mres <- mappend(mres, m, "bison (sediment)", "lat")
lres <- lrtappend(lres, m, "bison (sediment)", "lat")

m <- explore.continuous(site.sediment, "lon", "sex")
mres <- mappend(mres, m, "bison (sediment)", "lon")
lres <- lrtappend(lres, m, "bison (sediment)", "lon")

m <- explore.continuous(site.sediment, "alt", "sex")
mres <- mappend(mres, m, "bison (sediment)", "alt")
lres <- lrtappend(lres, m, "bison (sediment)", "alt")

m <- explore.continuous(site.sediment, "age", "sex")
mres <- mappend(mres, m, "bison (sediment)", "age")
lres <- lrtappend(lres, m, "bison (sediment)", "age")

m <- explore.continuous(site.sediment, "endog", "sex")
mres <- mappend(mres, m, "bison (sediment)", "endog")
lres <- lrtappend(lres, m, "bison (sediment)", "endog")

m <- explore.continuous(site.sediment, "GC", "sex")
mres <- mappend(mres, m, "bison (sediment)", "GC")
lres <- lrtappend(lres, m, "bison (sediment)", "GC")

m <- explore.continuous(site.sediment, "readlen", "sex")
mres <- mappend(mres, m, "bison (sediment)", "readlen")
lres <- lrtappend(lres, m, "bison (sediment)", "readlen")

m <- explore.continuous(site.sediment, "C2T", "sex")
mres <- mappend(mres, m, "bison (sediment)", "C2T")
lres <- lrtappend(lres, m, "bison (sediment)", "C2T")


postcrania = subset(d, material2=="Postcrania")
postcrania = droplevels(postcrania)
m.null = glm(sex~1, data=postcrania, family=binomial)
mres <- mappend(mres, m.null, "bison (postcrania)", "intercept-only")

m <- explore.categorical(postcrania, "material", "sex") # *
mres <- mappend(mres, m, "bison (postcrania)", "material")
lres <- lrtappend(lres, m, "bison (postcrania)", "material")

m <- explore.categorical(postcrania, "site_type", "sex") # *
mres <- mappend(mres, m, "bison (postcrania)", "site_type")
lres <- lrtappend(lres, m, "bison (postcrania)", "site_type")

m <- explore.categorical(postcrania, "loncat", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "loncat")
lres <- lrtappend(lres, m, "bison (postcrania)", "loncat")

m <- explore.continuous(postcrania, "lat", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "lat")
lres <- lrtappend(lres, m, "bison (postcrania)", "lat")

m <- explore.continuous(postcrania, "lon", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "lon")
lres <- lrtappend(lres, m, "bison (postcrania)", "lon")

m <- explore.continuous(postcrania, "alt", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "alt")
lres <- lrtappend(lres, m, "bison (postcrania)", "alt")

m <- explore.continuous(postcrania, "age", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "age")
lres <- lrtappend(lres, m, "bison (postcrania)", "age")

m <- explore.continuous(postcrania, "endog", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "endog")
lres <- lrtappend(lres, m, "bison (postcrania)", "endog")

m <- explore.continuous(postcrania, "GC", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "GC")
lres <- lrtappend(lres, m, "bison (postcrania)", "GC")

m <- explore.continuous(postcrania, "readlen", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "readlen")
lres <- lrtappend(lres, m, "bison (postcrania)", "readlen")

m <- explore.continuous(postcrania, "C2T", "sex")
mres <- mappend(mres, m, "bison (postcrania)", "C2T")
lres <- lrtappend(lres, m, "bison (postcrania)", "C2T")


mres.bison <- mres
lres.bison <- lres
