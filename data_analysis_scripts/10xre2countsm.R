#!/user/bin/Rscript
#convert 10x result to counts matrix
#self.R 10xresultdir *_counts.txt
args=commandArgs(T)
library(Seurat)
data <- Read10X(data.dir = args[1])
write.table(data,paste0(args[2],"_counts.txt"),sep="\t",quote=FALSE,col.names=NA)
