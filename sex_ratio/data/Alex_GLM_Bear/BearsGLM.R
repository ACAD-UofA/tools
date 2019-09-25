
#Remove anything left over from earlier sessions
rm(list=ls()) 

dir()


#load required packages. 
library(MuMIn)
library(heplots)
library(MASS)
library(car)


#The MuMIn package can be a bit finicky -, it seems to require that we change the "options" setting in our workspace to "na.fail".
#have to remove variables with NA
options(na.action = "na.fail")

#Read in the brown bear meta data file,
dat <- read.csv("SexRatioBrownBear.csv")


#plotting variables checking for normality in continuous variables 
plot(data.frame(dat$Region, dat$Latitude, dat$Site.type, dat$Bone.Type, dat$Collected.From, dat$Date, dat$Endogenous, dat$Clonal, dat$RatioMtNu))
plot(dat$Endogenous~dat$Sex.Assignment) #slightly skewed
plot(dat$Date~dat$Sex.Assignment, na.action="na.exclude") #skewed
plot(dat$Latitude~dat$Sex.Assignment) #seems fine
plot(dat$Clonal~dat$Sex.Assignment) #very skewed
plot(dat$RatioMtNu~dat$Sex.Assignment) #very skewed


#transformation to try and force normality #probably don't need to do this for glm
dat$logDate <- log(dat$Date+1) #+1 as historic samples and counted as 0
dat$logEndogenous <- log(dat$Endogenous)
dat$logClonal <- log(dat$Clonal)
dat$logRatioMtNu <- log(dat$RatioMtNu)

#plots of transformed variables
plot(dat$logEndogenous~dat$Sex.Assignment, na.action="na.exclude") #seems better
plot(dat$logDate~dat$Sex.Assignment, na.action="na.exclude") #questionable whether this is an improvement
plot(dat$logClonal~dat$Sex.Assignment, na.action="na.exclude") #seems better
plot(dat$logRatioMtNu~dat$Sex.Assignment, na.action="na.exclude") #seems better


#centering (makes interpretating interactions easier)
dat$lE <- as.numeric(scale(dat$logEndogenous, center = TRUE, scale = FALSE))
dat$lD <- as.numeric(scale(dat$logDate, center = TRUE, scale = FALSE))
dat$lC <- as.numeric(scale(dat$logClonal, center = TRUE, scale = FALSE))
dat$lR <- as.numeric(scale(dat$logRatioMtNu, center = TRUE, scale = FALSE))
dat$Lat <- as.numeric(scale(dat$Latitude, center = TRUE, scale = FALSE))


###Model averaging for generalised linear models


#Glm with all variables
MSex <- glm(Sex.Assignment ~ Site.type + Region + Bone.Type + lE + lD + lC + lR + Lat, family = binomial(link = "logit"), data = dat, na.action = "na.omit")
anova(MSex, test = "LRT")

#without variables with NAs/no interactions
MSex2 <- glm(Sex.Assignment ~ Site.type + Region + Bone.Type + lE + lC + lR + Lat, family = binomial(link = "logit"), data = dat, na.action = "na.fail")
anova(MSex2, test = "LRT")

#without variables with NAs/ with interactions, takes very long to compute
MSex3 <- glm(Sex.Assignment ~ (Site.type + Region + Bone.Type + lE + lC + lR + Lat)^2, family = binomial(link = "logit"), data = dat, na.action = "na.fail")



#Use AIC first to compare the models and find the one with the best AIC
extractAIC(MSex)
extractAIC(MSex2)
extractAIC(MSex3)
MSex2Step <- stepAIC(MSex2)

#The model averaging procedure using dredge (can't handle NAs)
MSex2Full <- dredge(MSex3, rank = "AIC", extra = "R^2")
MSex2Full
#other selection criteria
MSex2Avg <- model.avg(MSex2Full, subset = cumsum(weight) <= 0.95)
MSex2Avg$coefficients
MSex2Avg$importance

#Lindsey seems right in that I probably don't have to transform my predictors
#therefore:
#centering continuous predictors, shouldn't need to standardise, but can if needed (scale=TRUE)
dat$E <- as.numeric(scale(dat$Endogenous, center = TRUE, scale = FALSE))
dat$C <- as.numeric(scale(dat$Clonal, center = TRUE, scale = FALSE))
dat$R <- as.numeric(scale(dat$RatioMtNu, center = TRUE, scale = FALSE))
dat$D <- as.numeric(scale(dat$logDate, center = TRUE, scale = FALSE))

###
MSex4 <- glm(Sex.Assignment ~ (Site.type + Region + Bone.Type + E + C + R + Lat)^2, family = binomial(link = "logit"), data = dat, na.action = "na.fail")  
MSex4Full <- dredge(MSex4, rank = "AIC", extra = "R^2") #Takes impossibly long to run. Redo with just interactions interested for now
MSex4Full

MSex5 <- glm(Sex.Assignment ~ Site.type + Region + Bone.Type + E + C + R + Lat + Site.type:Region+Site.type:Bone.Type+Site.type:E+Site.type:Lat+Region:Bone.Type+Region:E+Region:Lat+Bone.Type:Lat+E:Lat, family = binomial(link = "logit"), data = dat, na.action = "na.fail") #only with interactions deemed relevant included to decrease compute time
MSex5Full <- dredge(MSex5, rank = "AIC", extra = "R^2")
MSex5Full 

MSex5Avg <- model.avg(MSex5Full, subset = cumsum(weight) <= 0.95)
MSex5Avg$coefficients
MSex5Avg$importance #4 variables (Site.type, Region, E, and Lat) have similar importance

MSexFinal <- glm(Sex.Assignment ~ Region + Site.type + E + Lat , data = dat, family = binomial(link = "logit")) #included Lat because it there is no major difference in AIC between the two best models, also it is included in a similar amount of models as the other predictors included when looking at importance

#This gives us a single model, on which we can look at parameter estimates, effect sizes (eta-squared), diagnostic plots, etc. 
coef(MSexFinal)
anova(MSexFinal, test = "LRT")
etasq(MSexFinal)
hist(resid(MSexFinal))
glm.diag(MSexFinal)
glm.diag.plots(MSexFinal) #not looking great
plot(fitted(MSexFinal), rstandard(MSexFinal)) # fitted vs residual. this looks weird  


#####PLOTS#######
counts1 <- table(dat$Sex.Assignment)
pcounts1 <- prop.table (counts1)
barplot(counts1, beside=T, legend= rownames(counts1), col=c("brown3","cadetblue"), main="Brown bear", ylim=c(0,80))
barplot(pcounts1, beside=T, legend= rownames(pcounts1), col=c("brown3","cadetblue"), main="Brown bear",ylim=c(0,1))

counts2 <- table(dat$Sex.Assignment,dat$Site.type)
pcounts2 <- prop.table (counts2,2)
barplot(counts2, beside=T, legend= rownames(counts2), col=c("brown3","cadetblue"), main="Brown bear", ylim=c(0,50))
barplot(pcounts2, beside=T, legend= rownames(pcounts2), col=c("brown3","cadetblue"), main="Brown bear",ylim=c(0,1.0))

counts3 <- table(dat$Sex.Assignment,dat$Region)
pcounts3 <- prop.table (counts3,2)
barplot(counts3, beside=T, legend= rownames(counts3), col=c("brown3","cadetblue"), main="Brown bear", ylim=c(0,50))
barplot(pcounts3, beside=T, legend= rownames(pcounts3), col=c("brown3","cadetblue"), main="Brown bear",ylim=c(0,1.0))
