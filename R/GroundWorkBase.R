row_merge <- function(data1,data2){
  samesample = intersect(rownames(data1),rownames(data2))
  data1 = data1[samesample,,drop=FALSE]
  data2 = data2[samesample,,drop=FALSE]
  data3 = cbind(data1,data2)
  return(data3)
}
row_same <- function(data1,data2){
  samesample = intersect(rownames(data1),rownames(data2))
  return(samesample)
}

`%!in%` <- Negate(`%in%`)

col_merge <- function(data1,data2){
  samesample = intersect(colnames(data1),colnames(data2))
  data1 = data1[,samesample,drop=FALSE]
  data2 = data2[,samesample,drop=FALSE]
  data3 = rbind(data1,data2)
  return(data3)
}

col_same <- function(data1,data2){
  samesample = intersect(colnames(data1),colnames(data2))
  return(samesample)
}

read_txt_numeric <- function(data){
library(limma)
data = read.table(data,sep="\t",header=TRUE,check.names=FALSE)
data = as.matrix(data)
rownames(data) = data[,1]
data = data[,-1]
class(data) = "numeric"
data = avereps(data)
return(data)
}

read_txt_charactor <- function(data){
  library(limma)
  data = read.table(data,sep="\t",header=TRUE,check.names=FALSE)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  return(data)
}
write_txt <- function(data,OutputName){
  out = rbind(symbol=colnames(data),data)
  write.table(out,file=paste0(OutputName,".txt"),
              sep="\t",
              col.names=F,
              quote=F)
}

write_csv <- function(data,OutputName,row_names = TRUE){
  write.csv(data,file=paste0(OutputName,".csv"),
              row.names = row_names,
              quote=F)
}