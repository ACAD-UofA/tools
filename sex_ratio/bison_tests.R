
library(dplyr)
library(reshape2)

source('explore.R')
source('mapstuff.R')
source('two-sample-tests.R')

d = read.table("data/bison.csv", sep="\t", header=TRUE)
mres = tibble(data=character(), model=character(),`coeff p`=numeric(), `LRT`=numeric())
i = 1

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
mres[i,] <- list("bison", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(d, "loncat", "sex")
mres[i,] <- list("bison", "loncat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "site_type", "sex") # *
mres[i,] <- list("bison", "site_type", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "material", "sex")
mres[i,] <- list("bison", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "material2", "sex")
mres[i,] <- list("bison", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

# check that the 'reference' level isn't fucking us
d.sc2.ref.post = transform(d, material2=addNA(relevel(material2, ref="Postcrania")))
explore.categorical(d.sc2.ref.post, "material2", "sex")

m <- explore.continuous(d, "lat", "sex")
mres[i,] <- list("bison", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "lon", "sex")
mres[i,] <- list("bison", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "alt", "sex")
mres[i,] <- list("bison", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "age", "sex")
mres[i,] <- list("bison", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "endog", "sex")
mres[i,] <- list("bison", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "GC", "sex")
mres[i,] <- list("bison", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "readlen", "sex")
mres[i,] <- list("bison", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "C2T", "sex")
mres[i,] <- list("bison", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1


ggplot(d, aes(x=lon, y=lat, colour=sex, shape=sex)) +
  geom_jitter(width=.2, height=.2, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue")) +
  geom_density_2d()

#nh = get_map(location=c(-180,20,180,75), maptype="watercolor")
#ggmap(nh) +
#  geom_jitter(data=d, aes(x=lon,y=lat,colour=sex,shape=sex), width=.2, height=.2, size=2) +
#  scale_shape_manual(values=c(1,4)) +
#  scale_color_manual(values=c("red", "blue"))

world=world.panned()
ggplot(world, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill="white", colour="black") +
  geom_jitter(data=d, inherit.aes=F, aes(x=lon,y=lat,colour=sex,shape=sex), width=.5, height=.5, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue")) +
  coord_map() #+
#  ylim(0,80)

# test for differences in spatial distribution between sexes
d.sex.sorted = subset(d[order(d$sex),], !is.na(lat) & !is.na(lon))
gcm = greatcircle.dmatrix(d.sex.sorted$lat, d.sex.sorted$lon)
kern1 = kernel.test(gcm, table(d.sex.sorted$sex))
kern1
hist(kern1$permutations,100)
abline(v=kern1$statistic, col="red")

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
mres[i,] <- list("bison (sediment)", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(site.sediment, "material", "sex")
mres[i,] <- list("bison (sediment)", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(site.sediment, "material2", "sex")
mres[i,] <- list("bison (sediment)", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(site.sediment, "loncat", "sex")
mres[i,] <- list("bison (sediment)", "loncat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

m <- explore.continuous(site.sediment, "lat", "sex")
mres[i,] <- list("bison (sediment)", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "lon", "sex")
mres[i,] <- list("bison (sediment)", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "alt", "sex")
mres[i,] <- list("bison (sediment)", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "age", "sex")
mres[i,] <- list("bison (sediment)", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "endog", "sex")
mres[i,] <- list("bison (sediment)", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "GC", "sex")
mres[i,] <- list("bison (sediment)", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "readlen", "sex")
mres[i,] <- list("bison (sediment)", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(site.sediment, "C2T", "sex")
mres[i,] <- list("bison (sediment)", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

postcrania = subset(d, material2=="Postcrania")
postcrania = droplevels(postcrania)
m.null = glm(sex~1, data=postcrania, family=binomial)
mres[i,] <- list("bison (postcrania)", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(postcrania, "material", "sex") # *
mres[i,] <- list("bison (postcrania)", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(postcrania, "site_type", "sex") # *
mres[i,] <- list("bison (postcrania)", "site_type", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(postcrania, "loncat", "sex")
mres[i,] <- list("bison (postcrania)", "loncat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

m <- explore.continuous(postcrania, "lat", "sex")
mres[i,] <- list("bison (postcrania)", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "lon", "sex")
mres[i,] <- list("bison (postcrania)", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "alt", "sex")
mres[i,] <- list("bison (postcrania)", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "age", "sex")
mres[i,] <- list("bison (postcrania)", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "endog", "sex")
mres[i,] <- list("bison (postcrania)", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "GC", "sex")
mres[i,] <- list("bison (postcrania)", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "readlen", "sex")
mres[i,] <- list("bison (postcrania)", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(postcrania, "C2T", "sex")
mres[i,] <- list("bison (postcrania)", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

# summarise

mres.bison <- mres
summary.bison <- dcast(mres, model~data, value.var='LRT')
summary.bison
