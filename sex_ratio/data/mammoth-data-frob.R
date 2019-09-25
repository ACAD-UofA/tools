# mammoth data from Pecnerova et al. 2017

# Length of chrX and autosome, from LoxAfr4.
Lx = 120050768
La = 3050898245

d1 = read.table("P_Mammoth.txt", sep="\t", header=T)
d1 = subset(d1, sex!="U")
d1 = droplevels(d1)

d2 = read.table("Pecnerova_TableS1.csv", sep=",", header=T,
                na.strings=c("NA", "nd", ""),
                col.names=c("sample","age","material","loc","mappedreads","chrX","chr8"))

d = merge(d1,d2,by="sample")

# Two samples are beyond max C14 date, and are coerced to NA. 
d$age = as.integer(as.character(d$age))

# Collapse all bones into the same category.
d$material2 = d$material
levels(d$material2) <- list(bone=c("bone", "femur", "humerus", "tibia", "vertebra"),
                            hair="hair", tusk="tusk", tooth="tooth")

d$loc2 = ifelse(d$loc=="Wrangel Island", "Wrangel", "Other")
d$loc2 = factor(d$loc2)

out = data.frame(sample=d$sample, sex=d$sex, age=d$age,
                 material=d$material, material2=d$material2,
                 loc=d$loc, loc2=d$loc2,
                 Nx=d$chrX, Na=d$mappedreads, Lx=Lx, La=La)

write.table(out, "mammoths.csv", sep="\t", quote=FALSE, row.names=FALSE)
