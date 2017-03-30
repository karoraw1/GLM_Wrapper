library(edgeR)
library(methods)
path <- commandArgs(TRUE)[1]
f_in <- commandArgs(TRUE)[2]

if (any(is.na(c(path, f_in)))){ 
	print("Improper commands provided")
	print(paste("Path:", path))
	print(paste("File In:", f_in))
	stop("Exiting")
}

infile = paste0(path, "/", f_in)
basename = sub("\\..*", "", f_in)
f_out = paste0(basename, "_vst.csv")
outfile = paste0(path, "/", f_out)

otuTable = read.csv(infile, header=T, row.names=1)
x = as(otuTable, "matrix")
y = edgeR::DGEList(counts = t(x), remove.zeros = TRUE)
z = edgeR::calcNormFactors(y, method = "TMM")
if (!all(is.finite(z$samples$norm.factors))) {
	stop("Something wrong with calcNormFactors, non-finite $norm.factors")
}
offset = t(z$counts) - 1
write.table(offset, outfile, row.names=T, na="", col.names=T, sep=",")



	

