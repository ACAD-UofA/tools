
d1 = read.table("meta.bears.txt", sep="\t", header=TRUE, na.strings = ("None"))
d2 = read.table("BrownBears.txt", sep="\t", header=TRUE, na.strings = ("None"))
d = merge(d1, d2, by="sample")

d3 = read.table("brownbear-gc.txt", col.names=c("sample","GC"))
d4 = read.table("brownbear-mean-len.txt", col.names=c("sample","readlen"))
d5 = read.table("brownbear-damage.txt", col.names=c("sample","C2T", "G2A"))
d = merge(merge(merge(d,d3),d4),d5)

# divide into Alps/Non-alps
d$alpine = factor((d$lat>45 & d$lat<50) & (d$lon>5 & d$lon<20), labels=c("nonAlps","Alps"))

write.table(d, "brownbears.csv", sep="\t", quote=FALSE, row.names=FALSE)
