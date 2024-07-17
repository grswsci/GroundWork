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
write.table(out,file=paste0(GEO_ID,".Matrix.txt"),sep="\t",col.names=F,quote=F)
return(data)
}

getClinical <- function(GEO_ID){
clinicallines <- readLines(paste0(GEO_ID,"_series_matrix.txt"))  
clinicallines_filtered <- clinicallines[substr(clinicallines,1,7)=="!Sample"] 
clinicaldata <- read.table(textConnection(paste(clinicallines_filtered, collapse="\n")), header=TRUE)  
clinicaldata <- t(clinicaldata)
write.table(clinicaldata,file=paste0(GEO_ID,".clinical.txt"),sep="\t",col.names=F,row.names=T,quote=F)
return(clinicaldata)
}

getZscore <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-") 
  rv <- sweep(rv, 1, rowsd, "/") 
  return(rv)
}

getWebsite <- function(){
browseURL("https://grswsci.top/", browser = getOption("browser"),encodeIfNeeded = FALSE)
}

getplot <- function(type = "pdf", ncol = 1) {
  library(ggplotify)
  library(cowplot)
  library(magick)
  library(pdftools)
  fnames <- Sys.glob(paste0("*.",type))
  if(type == "pdf"){
    p <- lapply(fnames,function(i){
      pn <- as.ggplot(image_read_pdf(i))
    })
  }else if(type %in% c("jpg","png","tiff")){
    p <- lapply(fnames,function(i){
      pn <- as.ggplot(image_read(i))
    })
  }
  
  plot_grid(plotlist = p, ncol = ncol)
  return(p)
}

getMerge <- function(Tumor,GEO_ID,clinicalFile){
library(limma)             
Matrix = paste0(GEO_ID,".Matrix.txt")     
data = read.table(Matrix,sep="\t",header=TRUE,check.names=FALSE)
data = as.matrix(data)
rownames(data) = data[,1]
data = data[,-1]
class(data) = "numeric"
data = avereps(data)
data = t(data)
clinical = read.table(clinicalFile, header=T, sep="\t", check.names=F, fileEncoding="GB18030",row.names=1)
clinical$CancerType = paste0(Tumor,".",GEO_ID)
clinical = clinical[,c(2,1)]
sameSample = intersect(row.names(data),row.names(clinical))
data = data[sameSample,]
clinical = clinical[sameSample,]
out = cbind(clinical,data)
out = cbind(id=row.names(out),out)
write.table(out,
            file=paste0(Tumor,".GEO.",GEO_ID,".",Tumor,".txt"),
            sep="\t",
            row.names=F,
            quote=F)
return(out)
}
getMerge2 <- function(Tumor,GEO_ID,clinicalFile){
  library(limma)             
  Matrix = paste0(GEO_ID,".Matrix.txt")     
  data = read.table(Matrix,sep="\t",header=TRUE,check.names=FALSE)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  class(data) = "numeric"
  data = avereps(data)
  data = t(data)
  clinical = read.table(clinicalFile, header=T, sep="\t", check.names=F, fileEncoding="GB18030",row.names=1)
  clinical$DiseaseType = paste0(Tumor,".",GEO_ID)
  clinical = clinical[,c(2,1)]
  sameSample = intersect(row.names(data),row.names(clinical))
  data = data[sameSample,]
  clinical = clinical[sameSample,]
  out = cbind(clinical,data)
  out = cbind(id=row.names(out),out)
  write.table(out,
              file=paste0(Tumor,".GEO.",GEO_ID,".",Tumor,".txt"),
              sep="\t",
              row.names=F,
              quote=F)
  return(out)
}
getGEO <- function(GEO_ID){
  options(timeout=1000000)
  if(nchar(GEO_ID) == 8){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,5),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
  }else if(nchar(GEO_ID) == 9){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,6),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
  }else if(nchar(GEO_ID) == 7){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,4),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
    }
}

getGEOs <- function(GEO_IDs){
  options(timeout=1000000)
  for (GEO_ID in GEO_IDs){
  if(nchar(GEO_ID) == 8){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,5),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
  }else if(nchar(GEO_ID) == 9){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,6),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
  }else if(nchar(GEO_ID) == 7){
    url = paste0("http://218.108.182.182:9000/pubmed/DownFile/GEO?fileName=/geo/series/",substr(GEO_ID,1,4),"nnn/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz")
    destfile <- paste0(GEO_ID,"_series_matrix.txt.gz")
    download.file(url, destfile, mode = "wb")  
  }
  }
}