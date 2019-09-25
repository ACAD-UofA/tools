
#library(dplyr)
#library(tibble)
#library(forcats)
#library(ggplot2)
library(tidyverse)

coll <- read.delim("data/collections-merged.tab", col.names=c("Genus", "Species", "collection.date", "collection", "sex"), stringsAsFactors=FALSE)
coll$binom <- paste(coll$Genus, coll$Species, sep=" ")
coll <- transform(coll, collection.date=as.integer(collection.date))
coll$collection.date <- ifelse(coll$collection.date<1500, NA, coll$collection.date)

# split samples into categories based on date
coll$date.class <- cut(coll$collection.date, c(0,1940,1960,2020), c("<1940","1940-1960",">1960"))

mammals <- coll %>%
#  group_by(binom, date.class, collection) %>%
#  group_by(binom, collection) %>%
  group_by(binom) %>%
  summarise(Genus=unique(Genus)[1],
            Species=unique(Species)[1],
            M=sum(sex=="M"),
            F=sum(sex=="F")) %>%
  filter(M+F>=100) %>%
  mutate(mratio=M/(M+F),
         # Beta CI
         mr.lwr=qbeta(0.025, M+0.5, F+0.5),
         mr.upr=qbeta(0.975, M+0.5, F+0.5))

pt <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt", quote="", header=TRUE, na.strings=c("NA","-999.00","-999"), stringsAsFactors=FALSE)
mammals <- left_join(mammals, pt, by=c("Genus"="MSW05_Genus","Species"="MSW05_Species"))

names(mammals)[names(mammals)=="MSW05_Order"] <- "Order"
names(mammals)[names(mammals)=="MSW05_Family"] <- "Family"
names(mammals)[names(mammals)=="MSW05_Genus"] <- "Genus"
names(mammals)[names(mammals)=="MSW05_Species"] <- "Species"
mammals$binom <- NULL
names(mammals)[names(mammals)=="MSW05_Binomial"] <- "Binom"

# Fill in missing Orders by matching to species of same genera that have info
for (genus in unique(mammals$Genus)) {
  x <- subset(mammals, Genus==genus)$Order
  order <- unique(as.character(na.omit(x)))
  if (length(order) == 1) {
    mammals[mammals$Genus==genus, "Order"] <- order[[1]]
  }
}
mammals[mammals$Genus=="Clethrionomys", "Order"] <- "Rodentia"
mammals[mammals$Genus=="Nanonycteris", "Order"] <- "Chiroptera"
mammals[mammals$Genus=="Neoplatymops", "Order"] <- "Chiroptera"
# remove bison beetle
mammals <- mammals[-which(mammals$Genus=="Prosopocoilus"), ]

# Natural ordering for the mammalian orders, from Wilson & Reeder 2005.
# Suggested by Kris.
orderorder = c("Monotremata", "Didelphimorphia", "Paucituberculata",
               "Microbiotheria", "Notoryctemorphia", "Dasyuromorphia",
               "Peramelemorphia", "Diprotodontia", "Afrosoricida",
               "Macroscelidea", "Tubulidentata", "Hyracoidea",
               "Proboscidea", "Sirenia", "Cingulata",
               "Pilosa", "Scandentia", "Dermoptera",
               "Primates", "Rodentia", "Lagomorpha", "Erinaceomorpha",
               "Soricomorpha", "Chiroptera", "Pholidota",
               "Carnivora", "Perissodactyla", "Artiodactyla",
               "Cetacea")

mammals <- transform(mammals,
                Genus=factor(Genus),
                Species=factor(Species),
                Order=factor(Order, levels=orderorder),
                Family=factor(Family),
                Binom=factor(Binom),
                X1.1_ActivityCycle=factor(X1.1_ActivityCycle, labels=c("nocturnal", "mid", "diurnal"), ordered=TRUE),
                X6.1_DietBreadth=factor(X6.1_DietBreadth, ordered=TRUE),
                X6.2_TrophicLevel=factor(X6.2_TrophicLevel, labels=c("herbivore","omnivore","carnivore"), ordered=TRUE),
                X12.1_HabitatBreadth=factor(X12.1_HabitatBreadth, ordered=TRUE),
                X12.2_Terrestriality=factor(X12.2_Terrestriality, labels=c("fossorial", "aboveground")),
                X24.1_TeatNumber=factor(X24.1_TeatNumber, ordered=TRUE)
                )

mammals$Binom <- fct_reorder(mammals$Binom, mammals$mratio)
mammals <- mammals[-which(is.na(mammals$Binom)),] # remove unknown taxa
mammals <- droplevels(mammals)
mammals <- as_tibble(mammals)
write.table(mammals, "mammal-sex-ratios.csv", sep="\t", quote=FALSE, row.names=FALSE)

#to.rm <- c()
#for (o in as.character(levels(mammals$Order))) {
#  if (sum(mammals$Order==o) < 5) {
#    to.rm = c(to.rm, which(mammals$Order==o))
#  }
#}
#mammals2 <- mammals[-to.rm,]
#mammals2 <- droplevels(mammals2)

theme_set(
  theme_bw() +
    theme(axis.title=element_text(size=12, face="bold")) +
    theme(axis.text=element_text(size=10, face="bold")) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(plot.title=element_text(size=12,face="bold",hjust=0.5))
)
scale <- 0.75*3
width <- scale*4
height <- scale*1.5

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

m2 <- mammals[mammals$mratio>0.8,]
summary(lm(mratio~Order, data=mammals))
ggplot(mammals, aes(Order, mratio)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_jitter(width=0.2, height=0, alpha=1, colour="black", shape=20, size=2) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  #geom_violin(alpha=0.5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("M/(M+F)") +
  #ylim(0,1) +
  ggtitle("Sex ratios in mammal collections") #+
  #annotate("text", x=m2$Order, y=m2$mratio+0.05, label=m2$Binom)

ggsave("museum-boxplot.pdf", width=width, height=height)

summary(lm(mratio~Order*date.class, data=mammals))
ggplot(mammals, aes(Order, mratio, colour=date.class)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_point(position=position_jitterdodge(jitter.width=0.1, jitter.height=0, dodge.width=0.75)) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("M/(M+F)") +
  ylim(0,1) +
  ggtitle("Sex ratios in mammal collections")

summary(lm(mratio~Order*collection, data=mammals))
ggplot(mammals, aes(Order, mratio, colour=collection)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_point(position=position_jitterdodge(jitter.width=0.1, jitter.height=0, dodge.width=0.75)) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab("M/(M+F)") +
  ylim(0,1) +
  ggtitle("Sex ratios in mammal collections")

summary(lm(mratio~Order*X1.1_ActivityCycle, data=mammals))
ggplot(mammals, aes(Order, mratio, colour=X1.1_ActivityCycle)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

summary(lm(mratio~Order*X6.2_TrophicLevel, data=mammals))
ggplot(mammals, aes(Order, mratio, colour=X6.2_TrophicLevel)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

summary(lm(mratio~Order*X12.1_HabitatBreadth, data=mammals))
ggplot(mammals, aes(Order, mratio, colour=X12.1_HabitatBreadth)) +
  geom_hline(yintercept=0.5, linetype=2, colour="gray") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, jitter.height=0, dodge.width=0.75)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("M/(M+F)") +
  ylim(0,1)

