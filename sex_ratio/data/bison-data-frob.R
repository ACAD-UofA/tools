
d1 = read.delim("bison.txt", sep="\t", header=TRUE, na.strings = ("None"),
                stringsAsFactors=FALSE)

# merge sample info for duplicates between Am and Eu subsets
for (i in which(duplicated(d1$sample))) {
  sample = as.character(d1[i,]$sample)
  x = which(d1$sample == d1[i,]$sample)[[1]] # master
  for (j in 1:ncol(d1)) {
    d1[x,j] = ifelse(is.na(d1[x,j]), d1[i,j], d1[x,j])
  }
}
d1 = d1[-which(duplicated(d1$sample)),]


d2 = read.table("bison-gc.txt", col.names=c("sample", "GC"))
d3 = read.table("bison-mean-len.txt", col.names=c("sample", "readlen"))
d4 = read.table("bison-damage.txt", col.names=c("sample", "C2T", "G2A"))

dd = merge(merge(merge(d1,d2), d3), d4)

write.table(dd, "bison-v2.csv", sep="\t", quote=FALSE, row.names=FALSE)
