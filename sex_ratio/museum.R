
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)

binom.split <- function(binom, delim="_") {
  fields <- strsplit(as.character(binom), delim)
  g = unlist(lapply(fields, function(x) x[[1]]))
  s = unlist(lapply(fields, function(x) x[[2]]))
  return(list(Genus=g, Species=s))
}

mammals <- read.table("data/museum.tab", sep="\t", col.names=c("taxid", "M", "F", "mratio"), stringsAsFactors=FALSE)
mammals <- subset(mammals, M+F>100)
mammals <- transform(mammals, mratio=M/(M+F))
mammals <- cbind(mammals, binom.split(mammals$taxid, " "), stringsAsFactors=FALSE)
mammals$taxid <- NULL

lh <- read.table("data/LH\ paper\ dataset.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
binom <- binom.split(lh$Species)
lh$Species <- NULL
lh <- cbind(lh, binom, stringsAsFactors=FALSE)
lh$Order <- NULL

sh <- read.table("data/Stockley&Hobson.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
binom <- binom.split(sh$Species)
sh$Species <- NULL
sh <- cbind(sh, binom, stringsAsFactors=FALSE)
sh$Order <- NULL

mam2 <- full_join(lh, sh, by=c("Genus", "Species"))
mammals <- inner_join(mam2, mammals, by=c("Genus", "Species"))

anage <- read.table("data/anage_data.txt", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
anage$Order <- NULL
mammals <- left_join(mammals, anage, by=c("Genus","Species"))

mom <- read.table("data/MOM_v4.1.csv", sep="\t", quote="", header=TRUE, na.strings=c("NA","-999.000"), stringsAsFactors=FALSE)
mom <- summarise(group_by(mom, Genus, Species), mom.mass=median(Mass,na.rm=TRUE))
mom$Order <- NULL
mammals <- left_join(mammals, mom, by=c("Genus","Species"))

pt <- read.table("data/PanTHERIA_1-0_WR05_Aug2008.txt", sep="\t", quote="", header=TRUE, na.strings=c("NA","-999.00","-999"), stringsAsFactors=FALSE)
mammals <- left_join(mammals, pt, by=c("Genus"="MSW05_Genus","Species"="MSW05_Species"))

# average (median) the body mass across the data sources
mammals$lbm <- log(apply(mammals[c("BodyMass", "BMf", "Adult.weight..g.", "mom.mass", "X5.1_AdultBodyMass_g")], 1, median, na.rm=TRUE))

# Beta CI
mammals$mr.lwr <- qbeta(0.025, mammals$M+0.5, mammals$F+0.5)
mammals$mr.upr <- qbeta(0.975, mammals$M+0.5, mammals$F+0.5)

names(mammals)[names(mammals)=="MSW05_Order"] <- "Order"
names(mammals)[names(mammals)=="MSW05_Family"] <- "Family"
names(mammals)[names(mammals)=="MSW05_Genus"] <- "Genus"
names(mammals)[names(mammals)=="MSW05_Species"] <- "Species"
names(mammals)[names(mammals)=="MSW05_Binomial"] <- "Binom"

mammals$PaternalCare <- ifelse(is.na(mammals$MCAll), mammals$MaleCare, mammals$MCAll)
# Do the datasets agree?
mammals[which(xor(mammals$MCAll, mammals$MaleCare)),
        c("Binom","Common.name","Order","MCAll","MaleCare","M","F")]
# set to NA where the datasets disagree
mammals$PaternalCare[which(xor(mammals$MCAll, mammals$MaleCare))] <- NA
mammals <- transform(mammals,
                Genus=factor(Genus),
                Species=factor(Species),
                Order=factor(Order),
                PaternalCare=factor(PaternalCare, labels=c("yes","no")),
                MCAll=factor(MCAll),
                MCc=factor(MCc),
                MCp=factor(MCp),
                MCpf=factor(MCpf),
                MCg=factor(MCg),
                MCh=factor(MCh),
                AlloAll=factor(AlloAll),
                SM=factor(SM),
                MaleCare=factor(MaleCare),
                MaleProv=factor(MaleProv),
                Monotocous=factor(Monotocous),
                X1.1_ActivityCycle=factor(X1.1_ActivityCycle, labels=c("nocturnal", "mid", "diurnal"), ordered=TRUE),
                X6.1_DietBreadth=factor(X6.1_DietBreadth, ordered=TRUE),
                X6.2_TrophicLevel=factor(X6.2_TrophicLevel, labels=c("herbivore","omnivore","carnivore"), ordered=TRUE),
                X12.1_HabitatBreadth=factor(X12.1_HabitatBreadth, ordered=TRUE),
                X12.2_Terrestriality=factor(X12.2_Terrestriality, labels=c("fossorial", "aboveground")),
                X24.1_TeatNumber=factor(X24.1_TeatNumber, ordered=TRUE)
                )

mammals <- transform(mammals, with_young=100*X25.1_WeaningAge_d/X14.1_InterbirthInterval_d)

mammals$Binom <- fct_reorder(mammals$Binom, mammals$mratio)
mammals <- droplevels(mammals)
mammals <- as_tibble(mammals)

to.rm <- c()
for (o in as.character(levels(mammals$Order))) {
  if (sum(mammals$Order==o) < 5) {
    to.rm = c(to.rm, which(mammals$Order==o))
  }
}
mammals2 <- mammals[-to.rm,]
mammals2 <- droplevels(mammals2)

ggplot(mammals, aes(y=Binom, x=mratio, colour=Order, shape=Order)) +
  geom_jitter(size=3, width=0.02, height=0) +
  geom_vline(xintercept=0.5, linetype=2, colour="red") +
  geom_errorbarh(aes(xmin=mr.lwr, xmax=mr.upr)) +
  xlim(0,1) +
  theme_bw() +
  xlab("M/(M+F)")

ggplot(mammals, aes(x=mratio)) +
  geom_density() +
  geom_rug() +
  geom_vline(xintercept=0.5, linetype=2, colour="red") +
  xlim(0,1) +
  theme_bw() +
  xlab("M/(M+F)")


m <- lm(mratio~MCAll, data=mammals)
summary(m)
m <- lm(mratio~MCp, data=mammals)
summary(m)
m <- lm(mratio~MCpf, data=mammals)
summary(m)
m <- lm(mratio~MCc, data=mammals)
summary(m)
m <- lm(mratio~MCg, data=mammals)
summary(m)
m <- lm(mratio~MCh, data=mammals)
summary(m)
m <- lm(mratio~AlloAll, data=mammals)
summary(m)
m <- lm(mratio~SM, data=mammals)
summary(m)

m <- lm(mratio~MaleCare, data=mammals)
summary(m)
m <- lm(mratio~MaleProv, data=mammals)
summary(m)
m <- lm(mratio~Monotocous, data=mammals)
summary(m)

m <- lm(mratio~PaternalCare, data=mammals)
summary(m)

m.full <- lm(mratio~Order+lbm+Order:lbm, data=mammals)
summary(m.full)

m2 = lm(mratio~lbm+Order, data=mammals)
summary(m2)
m3 = lm(mratio~lbm, data=mammals)
summary(m3)
m4 = lm(mratio~Order, data=mammals)
summary(m4)
TukeyHSD(aov(m4))

AIC(m.full,m2,m3,m4)
BIC(m.full,m2,m3,m4)
anova(m.full,m2,m3,m4,test="LRT")

m5 = lm(mratio~Order*PaternalCare, data=mammals)
summary(m5)

text_size=14
theme_set(
  theme_bw() +
    theme(legend.text=element_text(size=text_size, face="bold")) +
    theme(axis.title=element_text(size=text_size, face="bold")) +
    theme(axis.text=element_text(size=text_size, face="bold")) +
    theme(plot.title=element_text(hjust=0.5))
)

ggplot(mammals, aes(lbm, mratio)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  ylim(0,1) +
  xlab("log(Mass (g))") +
  ylab("M/(M+F)")

ggplot(mammals2, aes(lbm, mratio)) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  ylim(0,1) +
  xlab("log(Mass (g))") +
  ylab("M/(M+F)") +
  facet_grid(Order~.) +
  theme(strip.text = element_text(size=12, face="bold"))

ggplot(mammals2, aes(Order, mratio)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, height=0) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

ggplot(mammals2, aes(Order, mratio, colour=lbm)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1) +
  #guides(colour=guide_legend(title="log(Mass)")) + # This makes colours discrete. WTF?
  scale_colour_continuous(name="log(Mass)")

summary(lm(mratio~Order*PaternalCare, data=mammals))
ggplot(mammals2, aes(Order, mratio, colour=PaternalCare)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

summary(lm(mratio~Order*X1.1_ActivityCycle, data=mammals))
ggplot(mammals2, aes(Order, mratio, colour=X1.1_ActivityCycle)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

summary(lm(mratio~Order*X6.2_TrophicLevel, data=mammals))
ggplot(mammals2, aes(Order, mratio, colour=X6.2_TrophicLevel)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

summary(lm(mratio~Order*X12.1_HabitatBreadth, data=mammals))
ggplot(mammals2, aes(Order, mratio, colour=X12.1_HabitatBreadth)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

