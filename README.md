#安装GroundWork包

if (!require("devtools", quietly = TRUE))

  install.packages("devtools") 
  
devtools::install_github("grswsci/GroundWork") 


#调用GroundWork包

library(GroundWork)


#获得GEO探针名

GPL <- getGPL(GPL_ID = "GPL570") 


#输入文件为GSE63885_series_matrix.txt(压缩包的解压原始文件)，对矩阵文件进行注释(需先运行getGPL函数)

data <- getANN(GEO_ID = "GSE63885") 

#输入文件为GSE63885_series_matrix.txt(压缩包的解压原始文件)，获得临床文件

Clinical <- getClinical(GEO_ID = "GSE63885")


#对矩阵进行按行zscore

data.scale <- getZscore(data)

mean(data.scale[1,])

#1.436377e-15

var(data.scale[1,])

#1
