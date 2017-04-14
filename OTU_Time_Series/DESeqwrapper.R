library("DESeq2")
library(methods)

# read in path and otu table and ancillary data

path <- commandArgs(TRUE)[1]
f_in <- commandArgs(TRUE)[2]
f_cond <- commandArgs(TRUE)[3]

if (any(is.na(c(path, f_in)))){ 
	print("Improper commands provided")
	print(paste("Path:", path))
	print(paste("File In:", f_in))
	stop("Exiting")
}

infile = paste0(path, "/", f_in)
infile2 = paste0(path, "/", f_cond)
basename = sub("\\..*", "", f_in)
f_out = paste0(basename, "_vst.csv")
outfile = paste0(path, "/", f_out)

countData = as.matrix(read.csv(infile, header=TRUE, row.names=1))
colData = read.csv(infile2, row.names=1, header=FALSE)
colnames(colData) = c("Date")
colData$Date = as.factor(colData$Date)

dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Date)
rld <- rlog(dds, blind=FALSE)
to_write = as.data.frame(assay(rld))
# outfile = '/Users/login/Documents/GLM_Wrapper/OTU_Time_Series/rlog_saved.csv'
write.table(to_write, outfile, row.names=T, na="", col.names=T, sep=",")



