# Accumulate results for copy/pasting into LaTeX.
# All the code here is unfathomably bizarre and unintuitive;
# maybe the `tables' package provides an improvement?

source('bear_tests2.R')
source('bison_tests2.R')
source('mammoth_tests2.R')

mres <- rbind(mres.bears, mres.bison, mres.mammoths)
lres <- rbind(lres.bears, lres.bison, lres.mammoths)
#write.table(mres, file="mres.txt", row.names=FALSE,quote=FALSE,sep="\t")
write.table(lres, file="lres.txt", row.names=FALSE,quote=FALSE,sep="\t")
