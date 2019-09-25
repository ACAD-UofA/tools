# Accumulate results for copy/pasting into LaTeX.
# All the code here is unfathomably bizarre and unintuitive;
# maybe the `tables' package provides an improvement?

source('bear_tests.R')
source('bison_tests.R')
source('mammoth_tests.R')

mres <- rbind(mres.bears, mres.bison, mres.mammoths)
res <- dcast(mres, model~data, value.var='LRT')
ires <- dcast(mres, model~data, value.var='coeff p')

# Take the intercept-only coefficient, and put it in with the LRTs.
intercept.p <- ires[which(ires$model=="intercept-only"),2:ncol(ires)]
res[which(res$model=="intercept-only"),2:ncol(res)] <- intercept.p

# I want three significant figures, but format() interacts poorly with
# write.table()'s na parameter.
res2 <- format(res, digits=3)
res2[is.na(res)] <- NA

# Format output in a sensible way.
# The authors of format() and formatC() were fucking high!
fn <- function(x) {
  xf = as.numeric(x)
  fmt <- function(z) ifelse(!is.na(z), sprintf("%#.03G", z), NA)
  # Highlight statistically significant results.
  x <- ifelse(xf<0.05, sprintf("\\bfseries{\\textit{%s}}",fmt(xf)), fmt(xf))
  return (x)
}
res2[,2:ncol(res2)] <- apply(res2[,2:ncol(res2)], 2, fn)

# Reorder rows, and only take those of interest.
ord <- c("intercept-only", "site_type", "material", "material2",
         "age", "lat", "lon", "alt", "alpine", "endog", "GC", "readlen", "C2T")
ord.i <- unlist(lapply(ord, function(x) which(res2$model==x)))
res2 <- res2[ord.i,]

# Rename some variables.
res2[which(res2[,1]=="site_type"),1] <- "Cave/non-cave"
res2[which(res2[,1]=="lat"),1] <- "Latitude"
res2[which(res2[,1]=="lon"),1] <- "Longitude"
res2[which(res2[,1]=="alt"),1] <- "Altitude"
res2[which(res2[,1]=="alpine"),1] <- "Alps/non-Alps"
res2[which(res2[,1]=="endog"),1] <- "Endogenous"
res2[which(res2[,1]=="GC"),1] <- "GC ratio"
res2[which(res2[,1]=="readlen"),1] <- "DNA fragment length"
res2[which(res2[,1]=="C2T"),1] <- "5' C->T"

# Well that was an ordeal!
write.table(res2, row.names=FALSE, quote=FALSE, na="", sep=" & ", eol=" \\\\\n")
