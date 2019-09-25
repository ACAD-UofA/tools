
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

mres <- tibble(data=character(), model=character(),`coeff p`=numeric(), `LRT`=numeric())
i <- 1

# are sex ratios significantly different?
mf <- table(d$sex)
mf
mf["M"]/sum(mf)
binom.test(mf["M"], mf["M"]+mf["F"], alternative="greater")
m.null <- glm(sex~1, data=d, family=binomial)
summary(m.null)
mres[i,] <- list("brown bears", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(d, "site_type", "sex")
mres[i,] <- list("brown bears", "site_type", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "material", "sex")
mres[i,] <- list("brown bears", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d, "material2", "sex")
mres[i,] <- list("brown bears", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
#m <- explore.categorical(d, "collection", "sex")
#mres[i,] <- list("brown bears", "collection", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
#i <- i + 1

m <- explore.continuous(d, "lat", "sex")
mres[i,] <- list("brown bears", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "lon", "sex")
mres[i,] <- list("brown bears", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "alt", "sex") # *
mres[i,] <- list("brown bears", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "age", "sex") # *
mres[i,] <- list("brown bears", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "endog", "sex")
mres[i,] <- list("brown bears", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "GC", "sex")
mres[i,] <- list("brown bears", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "readlen", "sex")
mres[i,] <- list("brown bears", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d, "C2T", "sex")
mres[i,] <- list("brown bears", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

# evaluate (multi)collinearity of indep. vars.
car::scatterplotMatrix(d[c("lat", "lon", "alt", "age", "site_type", "material")])

# look at Variance Inflation Factors for full model
m <- glm(sex~lat+lon+alt+age+endog+material+site_type, family=binomial, data=d)
car::vif(m)

# remove site_type
#m = glm(sex~lat+lon+alt+age+endog+material, family=binomial, data=d)
#car::vif(m)


ggplot(d, aes(x=lon, y=alt, colour=sex)) + geom_jitter(width=2, height=30) + geom_density_2d()
ggplot(d, aes(x=age, y=alt, colour=sex)) + geom_jitter(width=1000, height=30) + geom_density_2d()

# strong geographic clustering is observable
ggplot(d, aes(x=lon, y=lat, colour=sex, shape=sex)) +
  geom_jitter(width=.2, height=.2, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue")) +
  geom_density_2d()

eu <- get_map(location=c(0,40,30,55), maptype="terrain")
ggmap(eu) +
  geom_jitter(data=d, aes(x=lon,y=lat,colour=sex,shape=sex), width=.1, height=.1, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue"))

nh <- get_map(location=c(-180,20,180,75), maptype="watercolor")
ggmap(nh) +
  geom_jitter(data=d, aes(x=lon,y=lat,colour=sex,shape=sex), width=.2, height=.2, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue"))

world <- world.panned()
ggplot(world, aes(x=long,y=lat,group=group)) +
  geom_polygon(fill="white", colour="black") +
  geom_jitter(data=d, inherit.aes=F, aes(x=lon,y=lat,colour=sex,shape=sex), width=.5, height=.5, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue")) +
  coord_map() +
  ylim(0,80)

# test for differences in spatial distribution between sexes
d.sex.sorted <- subset(d[order(d$sex),], !is.na(lat) & !is.na(lon))
gcm <- greatcircle.dmatrix(d.sex.sorted$lat, d.sex.sorted$lon)
kern1 <- kernel.test(gcm, table(d.sex.sorted$sex))
kern1
hist(kern1$permutations,100)
abline(v=kern1$statistic, col="red")


ggplot(d, aes(x=lon, y=lat, colour=sex, shape=alpine)) +
  geom_jitter(width=.2, height=.2, size=2) +
  scale_shape_manual(values=c(1,4)) +
  scale_color_manual(values=c("red", "blue")) +
  geom_density_2d()

m <- explore.categorical(d, "alpine", "sex") # *
mres[i,] <- list("brown bears", "alpine", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

explore.categorical(d, "material", "alpine")
explore.categorical(d, "material2", "alpine")

explore.continuous(d, "age", "alpine") # *
explore.continuous(d, "alt", "alpine") # *
explore.continuous(d, "lat", "alpine") # *
explore.continuous(d, "lon", "alpine")

explore.categorical(d, "alpine", "site_type")

d2 <- subset(d, alpine!="Alps")
d2 <- droplevels(d2)
d3 <- subset(d, alpine=="Alps")
d3 <- droplevels(d3)

m.null <- glm(sex~1, data=d2, family=binomial)
summary(m.null)
mres[i,] <- list("brown bears (non Alps)", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(d2, "material", "sex")
mres[i,] <- list("brown bears (non Alps)", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d2, "material2", "sex")
mres[i,] <- list("brown bears (non Alps)", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

m <- explore.categorical(d2, "site_type", "sex") # *
mres[i,] <- list("brown bears (non Alps)", "site_type", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
explore.categorical(d2, "material", "site_type")

m <- explore.continuous(d2, "lat", "sex") # *
mres[i,] <- list("brown bears (non Alps)", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "lon", "sex")
mres[i,] <- list("brown bears (non Alps)", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "alt", "sex")
mres[i,] <- list("brown bears (non Alps)", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "age", "sex")
mres[i,] <- list("brown bears (non Alps)", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "endog", "sex")
mres[i,] <- list("brown bears (non Alps)", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "GC", "sex")
mres[i,] <- list("brown bears (non Alps)", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "readlen", "sex")
mres[i,] <- list("brown bears (non Alps)", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d2, "C2T", "sex")
mres[i,] <- list("brown bears (non Alps)", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1


m.null <- glm(sex~1, data=d3, family=binomial)
summary(m.null)
mres[i,] <- list("brown bears (Alps)", "intercept-only", coefficients(summary(m.null))[1,4], NA)
i <- i + 1

m <- explore.categorical(d3, "material", "sex")
mres[i,] <- list("brown bears (Alps)", "material", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.categorical(d3, "material2", "sex")
mres[i,] <- list("brown bears (Alps)", "material2", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "alt", "sex")
mres[i,] <- list("brown bears (Alps)", "alt", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "age", "sex")
mres[i,] <- list("brown bears (Alps)", "age", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "endog", "sex")
mres[i,] <- list("brown bears (Alps)", "endog", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "lat", "sex")
mres[i,] <- list("brown bears (Alps)", "lat", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "lon", "sex")
mres[i,] <- list("brown bears (Alps)", "lon", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "GC", "sex")
mres[i,] <- list("brown bears (Alps)", "GC", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "readlen", "sex")
mres[i,] <- list("brown bears (Alps)", "readlen", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1
m <- explore.continuous(d3, "C2T", "sex")
mres[i,] <- list("brown bears (Alps)", "C2T", coefficients(summary(m))[2,4], anova(m,test="LRT")$"Pr(>Chi)"[[2]])
i <- i + 1

mres.bears <- mres
summary.bears <- dcast(mres, model~data, value.var='LRT')
summary.bears

# look at Variance Inflation Factors for full model
m <- glm(sex~lat+lon+alt+age+endog+material+site_type, family=binomial, data=d2)
car::vif(m)

m <- glm(sex~site_type+lat, data=d2, family=binomial)
anova(m, test="LRT")
car::vif(m)

m <- glm(sex~lat+site_type, data=d2, family=binomial)
anova(m, test="LRT")

m.null = update(m, sex~1, data=m$model)
print(anova(m.null, m, test="LRT"))

d4 <- subset(d2, !is.na(site_type) && !is.na(lat))
m1 <- glm(sex~site_type, data=d4, family=binomial)
m2 <- glm(sex~lat, data=d4, family=binomial)
anova(m1,m2, test="LRT")
anova(m2,m1, test="LRT")

explore.continuous(d2, "lat", "site_type")


# test for differences in spatial distribution between sexes
d2.sex.sorted <- subset(d2[order(d2$sex),], !is.na(lat) & !is.na(lon))
gcm2 <- greatcircle.dmatrix(d2.sex.sorted$lat, d2.sex.sorted$lon)
kern2 <- kernel.test(gcm2, table(d2.sex.sorted$sex))
kern2
hist(kern2$permutations,100)
abline(v=kern2$statistic, col="red")
