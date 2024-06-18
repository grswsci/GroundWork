#安装GroundWork包

if (!require("devtools", quietly = TRUE))

  install.packages("devtools") 
  
devtools::install_github("grswsci/GroundWork") 


#调用GroundWork包

library(GroundWork)

#获得GEO探针名

GPL <- getGPL(GPL_ID = "GPL570") 

head(GPL)

#ID            Gene Symbol

#1007_s_at     DDR1 /// MIR4640

#1053_at       RFC2

#117_at        HSPA6

#121_at        PAX8

#1255_g_at     GUCA1A

#1294_at       MIR5193 /// UBA7

#输入文件为GSE63885_series_matrix.txt(压缩包的解压原始文件)

#对矩阵文件进行注释(需先运行getGPL函数)

data <- getANN(GEO_ID = "GSE63885") 

data[1:3,1:3]

#                 GSM1559299 GSM1559300 GSM1559301

#DDR1 /// MIR4640  11.255891  11.543291  11.766435

#RFC2               9.021955   8.478034   9.449631

#HSPA6              7.548759   8.075928   8.131847

#输入文件为GSE63885_series_matrix.txt(压缩包的解压原始文件)

#获得临床文件

Clinical <- getClinical(GEO_ID = "GSE63885")


#对矩阵进行按行zscore

data.scale <- getZscore(data)

mean(data.scale[1,])

#1.436377e-15

var(data.scale[1,])

#1
