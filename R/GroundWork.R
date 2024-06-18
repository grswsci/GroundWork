getGPL <- function(GPL_ID){
  data_file <- system.file("data", paste0(GPL_ID,".RDS"), package = "GroundWork")
  GPL <- readRDS(data_file)
  GPL <- GPL[GPL$`Gene Symbol`!="",]
  return(GPL)
}

getANN <- function(GEO_ID){
lines <- readLines(paste0(GEO_ID,"_series_matrix.txt"))  
lines_filtered <- lines[!grepl("!", lines)]  
data <- read.table(textConnection(paste(lines_filtered, collapse="\n")), header=TRUE)  
GPL$ID <- paste0("GPL",GPL$ID)
data$ID_REF <- paste0("GPL",data$ID_REF)
suppressPackageStartupMessages(library(tidyverse))
gene_mapping <- as.vector(GPL$`Gene Symbol`) 
names(gene_mapping) <- GPL$ID
head(gene_mapping)
data$symbol <- unlist(lapply(data$ID_REF, function(x) gene_mapping[x] %||% x))
data <- na.omit(data)
data$ID_REF <- data$symbol
data <- data[,-ncol(data)]
data <- as.matrix(data)
rownames(data) <- data[,1]
data <- data[,-1]
class(data) <- "numeric"
data <- na.omit(data)
suppressPackageStartupMessages(library(limma))
data <- avereps(data)
out=rbind(symbol=colnames(data),data)
write.table(out,file=paste0(GEOID,".Matrix.txt"),sep="\t",col.names=F,quote=F)
return(data)
}

getClinical <- function(GEO_ID){
clinicallines <- readLines(paste0(GEO_ID,"_series_matrix.txt"))  
clinicallines_filtered <- clinicallines[substr(clinicallines,1,7)=="!Sample"] 
clinicaldata <- read.table(textConnection(paste(clinicallines_filtered, collapse="\n")), header=TRUE)  
clinicaldata <- t(clinicaldata)
write.table(clinicaldata,file=paste0(GEOID,".clinical.txt"),sep="\t",col.names=F,row.names=T,quote=F)
return(clinicaldata)
}

zscore <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-") 
  rv <- sweep(rv, 1, rowsd, "/") 
  return(rv)
}
