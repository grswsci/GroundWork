`%!in%` <- function(x, table) {  
  !(x %in% table)  
} 

getPanRNA <- function(){
  library(data.table)
  panexpr <- fread("/mnt/ExpBig/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.grsw.impute.log2(Xsum1).tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
  panexpr <- as.data.frame(panexpr)
  rownames(panexpr) <- panexpr[,1]
  panexpr <- panexpr[,-1]
  colnames(panexpr) <- substr(colnames(panexpr),1,15)
  return(panexpr)
}
getPansample <- function(){
  rawAnno <- read.delim("/mnt/Merge/merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
  rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
  samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
  samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
  return(samAnno)
}
getPanGeneZscore <- function(geneName){
  data = read.table(paste0("/mnt/expTime/Expression.",geneName,".zscore.txt"),sep="\t",header=T,check.names=F, row.names=1)  
  return(data)
}
getPanGeneZscore_SingleCancer <- function(geneName,CancerType){
  data = read.table(paste0("/mnt/expTime/Expression.",geneName,".zscore.txt"), header=T, sep="\t", check.names=F,row.names = 1)
  data = data[which(data$CancerType == CancerType),]
  return(data)
}
getPanGeneTime <- function(geneName){
  data = read.table(paste0("/mnt/expTime/expTime.",geneName,".txt"), header=T, sep="\t", check.names=F,row.names = 1)
  return(data)
}
getPanGeneTime_SingleCancer <- function(geneName,CancerType){
  data = read.table(paste0("/mnt/expTime/expTime.",geneName,".txt"), header=T, sep="\t", check.names=F,row.names = 1)
  data = data[which(data$CancerType == CancerType),]
  return(data)
}
getPanGeneExp <- function(geneName){
  data = read.table(paste0("/mnt/expTime/Expression.",geneName,".txt"), header=T, sep="\t", check.names=F,row.names = 1)
  return(data)
}
getPanGeneExp_SingleCancer <- function(geneName,CancerType){
  data = read.table(paste0("/mnt/expTime/Expression.",geneName,".txt"), header=T, sep="\t", check.names=F,row.names = 1)
  data = data[which(data$CancerType == CancerType),]
  return(data)
}
getscRNA <- function(CancerName){
  scRNA = readRDS(paste0("/mnt/TISCH/",CancerName,".RDS"))
  #scRNA[["RNA"]] = as(scRNA[["RNA"]], "Assay")
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data)=="Celltype (major-lineage)"] = "CellType_MainLabel"
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data)=="Celltype (minor-lineage)"] = "CellType_FineLabel"
  Idents(scRNA) = scRNA$CellType_MainLabel
  return(scRNA)
}

getRNA <- function(CancerName){
  expFile = paste0("/mnt/ExpBig/Expression.Firehose.",CancerName,".txt")
  library(limma)
  library(data.table)
  data = fread(expFile,sep="\t",header=TRUE,check.names=FALSE)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  class(data) = "numeric"
  return(data)
}

getRNA_GEO_Tumor <- function(CancerName,DatasetName){
  load(paste0("/mnt/Best/",CancerName,"/symbol.rda"))
  library(limma)
  library(data.table)
  data = total_expr_list[[DatasetName]]
  clinical = total_clin_list[[DatasetName]]
  if("Tissue" %in% colnames(clinical)){
    clinical = clinical[which(clinical$Tissue == "Tumor"),]
  }
  data_Tumor = data[,rownames(clinical)]
  return(data_Tumor)
}

getstRNA <- function(CancerName){
  stRNA = readRDS(paste0("/mnt/ST/",CancerName,".stRNA.RDS"))
  return(stRNA)
}

getProtein <- function(CancerName){
  library(data.table)
  DataFile = paste0("/mnt/Proteomics/expression_",CancerName,"_protein.csv")
  data = fread(DataFile)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  return(data)
}
getProteinImpute <- function(CancerName){
  library(data.table)
  library(limma)
  library(impute)
  DataFile = paste0("/mnt/Proteomics/expression_",CancerName,"_protein.csv")
  data = fread(DataFile)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  class(data) = "numeric"
  mat=impute.knn(data)
  data = mat$data
  return(data)
}
getZscore <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #表达量-均值
  rv <- sweep(rv, 1, rowsd, "/")  #再除以标准差
  return(rv)
}
getZscore_col <- function(data,GeneName,FilterOut = TRUE) {
  data$Zscore = as.numeric(scale(data[,GeneName], center = TRUE, scale = TRUE))
  if(FilterOut == TRUE){
    data = data[which(data$Zscore <3 & data$Zscore > -3),]
    return(data)
  }else{
    return(data)
  }
}
getscTumor <- function(scRNA) {
  if("Tissue" %in% colnames(scRNA@meta.data)){
    scRNA@meta.data$Tissue <- gsub("tumor","Tumor",scRNA@meta.data$Tissue)
    if("Tumor" %in% unique(scRNA@meta.data$Tissue)){
      scRNA <- scRNA[,scRNA@meta.data$Tissue== "Tumor"]
    }
    gc()
  }else if("Source" %in% colnames(scRNA@meta.data)){
    scRNA@meta.data$Source <- gsub("tumor","Tumor",scRNA@meta.data$Source)
    if("Tumor" %in% unique(scRNA@meta.data$Source)){
      scRNA <- scRNA[,scRNA@meta.data$Source== "Tumor"]
    }
    gc()
  }else{
    scRNA <- scRNA
    print(c("Nothing"))
  }
  return(scRNA)
}
getTime <- function(){
  startTime <- Sys.time()
  startTime <- as.character(startTime)
  startTime <- gsub(" " ,"-",startTime)
  startTime <- gsub(":" ,"-",startTime)
  return(startTime)
}
getForest <- function(data){
  data$HR = round(unlist(as.numeric(data$HR)),3)
  data$HR.95L = round(unlist(as.numeric(data$HR.95L)),3)
  data$HR.95H = round(unlist(as.numeric(data$HR.95H)),3)
  data$pvalue = round(unlist(as.numeric(data$pvalue)),3)
  data$pvalue = ifelse(data$pvalue<0.001,"<0.001",data$pvalue)
  lineVec = nrow(data)+1
  data$' ' <- paste(rep(" ", 10), collapse = " ")
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(forestploter))
  tm <- forest_theme(base_size = 18,   #图形整体的大小
                     #可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                     ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2, 
                     # 参考线形状、宽度、颜色
                     refline_lty="dashed", refline_lwd=1, refline_col="grey20",
                     #x轴刻度字体的大小
                     xaxis_cex=0.8,
                     #脚注大小、颜色
                     footnote_cex = 0.6, footnote_col = "blue")
  #绘制图形
  if(max(data$HR.95H)<1){
    xlim = c(min(data$HR.95L),1.2)
  }else if(max(data$HR.95H)<2){
    xlim = c(min(data$HR.95L),2)
  }else if(max(data$HR.95H)<3){
    xlim = c(min(data$HR.95L),3)
  }else if(max(data$HR.95H)<4){
    xlim = c(min(data$HR.95L),4)
  }else{
    xlim = c(min(data$HR.95L),5)
  }
  plot <- forestploter::forest(data[,c(1:(ncol(data)-2),ncol(data),(ncol(data)-1))],
                               est = as.numeric(data$HR),
                               lower = as.numeric(data$HR.95L),
                               upper = as.numeric(data$HR.95H),
                               ci_column = 5,     #可信区间所在的列
                               ref_line = 1,      #参考线条的位置
                               xlim = xlim,    #X轴的范围
                               #ticks_at = c(0,0.5, 1, 5,10),
                               theme = tm,        #图形的参数
  )
  
  boxcolor = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")
  
  boxcolor =ifelse(data$HR>1,boxcolor[1],boxcolor[2])
  #boxcolor = boxcolor[as.numeric(as.factor(data$Method))]
  for(i in 1:nrow(data)){
    plot <- edit_plot(plot, col=5,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25)) # 改col，box的列
  }
  #设置pvalue的字体
  pos_bold_pval = which(as.numeric(gsub('<',"",data$pvalue))<0.05)
  if(length(pos_bold_pval)>0){
    for(i in pos_bold_pval){
      plot <- edit_plot(plot, col=c(1,2,3,4,6),row = i, which = "text", gp = gpar(fontface="bold"))  # 改col pvalue的列
    }
  }
  #在图形中增加线段
  plot <- add_border(plot, part = "header", row =1,where = "top",gp = gpar(lwd =2))
  plot <- add_border(plot, part = "header", row = c(lineVec), gp = gpar(lwd =1))
  #设置字体的大小, 并且将文字居中
  plot <- edit_plot(plot, col=1:ncol(data),row = 1:nrow(data), which = "text", gp = gpar(fontsize=12))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="header",
                    x = unit(0.5, "npc"))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),
                    x = unit(0.5, "npc"))
  #输出图形
  pdf("Forest.pdf", width=12, heigh=10)
  print(plot)
  dev.off()
}

getForest_rbind <- function(data,heigh=10){
  data$HR = round(unlist(as.numeric(data$HR)),3)
  data$HR.95L = round(unlist(as.numeric(data$HR.95L)),3)
  data$HR.95H = round(unlist(as.numeric(data$HR.95H)),3)
  data$pvalue = round(unlist(as.numeric(data$pvalue)),3)
  data$pvalue = ifelse(data$pvalue<0.001,"<0.001",data$pvalue)
  lineVec = nrow(data)+1
  data$' ' = paste(rep(" ", 10), collapse = " ")
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(forestploter))
  tm <- forest_theme(base_size = 18,   #图形整体的大小
                     #可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                     ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2, 
                     # 参考线形状、宽度、颜色
                     refline_lty="dashed", refline_lwd=1, refline_col="grey20",
                     #x轴刻度字体的大小
                     xaxis_cex=0.8,
                     #脚注大小、颜色
                     footnote_cex = 0.6, footnote_col = "blue")
  #绘制图形
  if(max(data$HR.95H)<1){
    xlim = c(min(data$HR.95L),1.2)
  }else if(max(data$HR.95H)<2){
    xlim = c(min(data$HR.95L),2)
  }else if(max(data$HR.95H)<3){
    xlim = c(min(data$HR.95L),3)
  }else if(max(data$HR.95H)<4){
    xlim = c(min(data$HR.95L),4)
  }else{
    xlim = c(min(data$HR.95L),5)
  }
  plot <- forestploter::forest(data[,c(1:(ncol(data)-3),ncol(data),(ncol(data)-2):(ncol(data)-1))],
                               est = as.numeric(data$HR),
                               lower = as.numeric(data$HR.95L),
                               upper = as.numeric(data$HR.95H),
                               ci_column = 5,     #可信区间所在的列
                               ref_line = 1,      #参考线条的位置
                               xlim = xlim,    #X轴的范围
                               #ticks_at = c(0,0.5, 1, 5,10),
                               theme = tm,        #图形的参数
  )
  
  boxcolor = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")
  
  boxcolor =ifelse(data$HR>1,boxcolor[1],boxcolor[2])
  #boxcolor = boxcolor[as.numeric(as.factor(data$Method))]
  for(i in 1:nrow(data)){
    plot <- edit_plot(plot, col=5,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25)) # 改col，box的列
  }
  #设置pvalue的字体
  pos_bold_pval = which(as.numeric(gsub('<',"",data$pvalue))<0.05)
  if(length(pos_bold_pval)>0){
    for(i in pos_bold_pval){
      plot <- edit_plot(plot, col=c(1,2,3,4,6),row = i, which = "text", gp = gpar(fontface="bold"))  # 改col pvalue的列
    }
  }
  #在图形中增加线段
  plot <- add_border(plot, part = "header", row =1,where = "top",gp = gpar(lwd =2))
  plot <- add_border(plot, part = "header", row = c(lineVec), gp = gpar(lwd =1))
  #设置字体的大小, 并且将文字居中
  plot <- edit_plot(plot, col=1:ncol(data),row = 1:nrow(data), which = "text", gp = gpar(fontsize=12))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="header",
                    x = unit(0.5, "npc"))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),
                    x = unit(0.5, "npc"))
  #输出图形
  pdf("Forest.pdf", width=12, heigh=heigh)
  print(plot)
  dev.off()
}

getunicox <- function(data){
  unicox_merge = data.frame()
  library(survival)
  for (variable in colnames(data)[3:ncol(data)]){
    data_subset = data[,c(colnames(data)[1:2],variable)]
    unicox = as.formula(paste0("Surv(",colnames(data_subset)[1],",",colnames(data_subset)[2],")","~",variable))
    myunicox = coxph(unicox, data = data_subset)
    unicoxResult = summary(myunicox)
    unicox_merge = rbind(unicox_merge,
                         cbind(id = variable,
                               HR = round(unicoxResult$conf.int[,"exp(coef)"],3),
                               HR.95L = round(unicoxResult$conf.int[,"lower .95"],3),
                               HR.95H = round(unicoxResult$conf.int[,"upper .95"],3),
                               pvalue = unicoxResult$coefficients[,"Pr(>|z|)"])
    )
  }
  return(unicox_merge)
}

getmulticox <- function(data){
  library(survival)
  multicoxTable = data.frame()
  multicox = as.formula(paste0("Surv(",colnames(data)[1],",",colnames(data)[2],")","~","."))
  multicox = coxph(multicox, data = data)
  multicoxsummary = summary(multicox)
  multicoxTable = rbind(multicoxTable,
                        cbind(id = rownames(multicoxsummary[["coefficients"]]),
                              HR = round(multicoxsummary$conf.int[,"exp(coef)"],3),
                              HR.95L = round(multicoxsummary$conf.int[,"lower .95"],3),
                              HR.95H = round(multicoxsummary$conf.int[,"upper .95"],3),
                              pvalue = multicoxsummary$coefficients[,"Pr(>|z|)"])
  )
  multicoxTable = as.data.frame(multicoxTable)
  return(multicoxTable)
}
getRNAmultigene <- function(CancerName,multigene){
  expFile = paste0("/mnt/ExpBig/Expression.Firehose.",CancerName,".txt")
  library(limma)
  library(data.table)
  data = fread(expFile,sep="\t",header=TRUE,check.names=FALSE)
  data = as.matrix(data)
  rownames(data) = data[,1]
  data = data[,-1]
  class(data) = "numeric"
  data = data[intersect(as.vector(multigene),rownames(data)),]
  return(data)
}
getsurvival <- function(){
  library(data.table)
  data = fread("/mnt/Merge/Survival_SupplementalTable_S1_20171025_xena_sp")
  data = as.data.frame(data)
  data = data[which(substr(data$sample,14,14) == "0"),]
  data = data[,c("sample","cancer type abbreviation","OS.time","OS","DSS.time","DSS","PFI.time","PFI","DFI.time","DFI")]
  colnames(data)[2] = "CancerType"
  rownames(data) = data[,1]
  return(data)
}

getsurvivalfilter <- function(data,SurvivalType){
  if(SurvivalType == "OS"){
    data$OS.time = data$OS.time/365
    data = as.data.frame(data)
    data = data[data$OS.time>0,]
    data = data[!(is.na(data$OS.time)),]
    data = data[!(is.na(data$OS)),]
    return(data)
  }else if(SurvivalType == "DSS"){
    data$DSS.time = data$DSS.time/365
    data = as.data.frame(data)
    data = data[data$DSS.time>0,]
    data = data[!(is.na(data$DSS.time)),]
    data = data[!(is.na(data$DSS)),]
    return(data)
  }else if(SurvivalType == "DFI"){
    data$DFI.time = data$DFI.time/365
    data = as.data.frame(data)
    data = data[data$DFI.time>0,]
    data = data[!(is.na(data$DFI.time)),]
    data = data[!(is.na(data$DFI)),]
    return(data)
  }else if(SurvivalType == "PFI"){
    data$PFI.time = data$PFI.time/365
    data = as.data.frame(data)
    data = data[data$PFI.time>0,]
    data = data[!(is.na(data$PFI.time)),]
    data = data[!(is.na(data$PFI)),]
    return(data)
  }
}

getmerge_row <- function(data1,data2){
  samesample = intersect(rownames(data1),rownames(data2))
  data1 = data1[samesample,,drop=FALSE]
  data2 = data2[samesample,,drop=FALSE]
  data3 = cbind(data1,data2)
  return(data3)
}

getwilcoxonggplot2_another_two <- function(data,
                                           Type = "Type",
                                           variable,
                                           leadin = "",
                                           levels = NULL,
                                           mycolor = c("#0072B5FF","#BC3C29FF","#E18727FF",
                                                       "#20854EFF","#7876B1FF","#6F99ADFF", 
                                                       "#FFDC91FF","#EE4C97FF","#E64B35FF",
                                                       "#4DBBD5FF","#00A087FF","#3C5488FF",
                                                       "#F39B7FFF","#8491B4FF","#91D1C2FF",
                                                       "#DC0000FF","#7E6148FF","#B09C85FF",
                                                       "#3B4992FF","#EE0000FF","#008B45FF",
                                                       "#631879FF","#008280FF","#BB0021FF",
                                                       "#5F559BFF","#A20056FF","#808180FF",
                                                       "#00468BFF","#ED0000FF","#42B540FF", 
                                                       "#0099B4FF","#925E9FFF","#FDAF91FF",
                                                       "#AD002AFF","#ADB6B6FF","#374E55FF",
                                                       "#DF8F44FF","#00A1D5FF","#B24745FF",
                                                       "#79AF97FF","#6A6599FF","#80796BFF",
                                                       "#1f77b4",  "#ff7f0e",  "#279e68",
                                                       "#d62728",  "#aa40fc",  "#8c564b",
                                                       "#e377c2",  "#b5bd61",  "#17becf",
                                                       "#aec7e8"),
                                           width = 6,
                                           height = 4,
                                           CancerName){
  library(beeswarm)
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  p = wilcox.test(expression ~ Type, data = data)$p.value
  if(p < 0.001){
    pval = signif(p,4)
    pval = format(pval, scientific = TRUE)
  }else{
    pval = sprintf("%.03f",p)
  }
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  box_plot = boxplot(expression ~ Type, data = data,outline = FALSE, plot=F)
  yMin = min(box_plot$stats)
  yMax = max(box_plot$stats/5 + box_plot$stats)
  n = ncol(box_plot$stats)
  pdf(file = paste0(leadin,"_",CancerName,"_",variable,"_wilcoxon.pdf"), width = width, height = height)
  par(mar = c(4.5,6,3,3))
  boxplot(expression ~ Type, 
          data = data,
          ylab = variable,
          xlab = leadin,
          main = paste0("Wilcoxon Rank Sum Tests (p=",pval,")"),
          cex.main = 1.4, 
          cex.lab = 1.4, 
          cex.axis = 1.3,
          ylim = c(yMin,yMax),
          outline = FALSE)
  beeswarm(expression ~ Type, 
           data = data, 
           col = mycolor, 
           lwd = 0.1,
           pch = 16, 
           add = TRUE, 
           corral="wrap")
  dev.off()
}  

getwilcoxonggplot2 <- function(data,
                               Type = "Type",
                               variable,
                               levels = NULL,
                               mycolor = c("#0072B5FF","#BC3C29FF"),
                               DatasetName){
  suppressPackageStartupMessages(library(ggplot2,quietly = TRUE))
  suppressPackageStartupMessages(library(ggpubr,quietly = TRUE))
  suppressPackageStartupMessages(library(limma,quietly = TRUE))
  suppressPackageStartupMessages(library(ggExtra,quietly = TRUE))
  suppressPackageStartupMessages(library(reshape2,quietly = TRUE))
  suppressPackageStartupMessages(library(data.table,quietly = TRUE))
  suppressPackageStartupMessages(library(aplot,quietly = TRUE))
  suppressPackageStartupMessages(library(dplyr,quietly = TRUE))
  suppressPackageStartupMessages(library(hrbrthemes,quietly = TRUE))
  suppressPackageStartupMessages(library(ggtext,quietly = TRUE))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  p <- wilcox.test(expression ~ Type, data = data)$p.value
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  # 使用ggplot2函数进行绘图，数据来源于tmp数据框，x轴为Grade列，y轴为tmp数据框中由genName变量指定的列，并且使用Grade列作为填充色
  ggplot(data = data,
         aes(x = Type, 
             y = expression, 
             fill = Type)) +
    # 设置填充颜色
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    scale_color_manual(values = mycolor)+
    scale_fill_manual(values = mycolor) + 
    # 绘制小提琴图，用于显示数据的分布，设置透明度、位置、线条宽度和颜色 
    geom_violin(alpha=0.4, position = position_dodge(width = .75), size=0.8, color="black") + 
    # 在小提琴图上添加箱线图，显示中位数、四分位数等，设置缺口、异常点大小、线条宽度、透明度和颜色  
    geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.5) +
    # 在箱线图的基础上添加点图，使用21号形状（带边框的圆形），设置大小、位置、颜色和透明度  
    geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha=0.3) +
    # 使用经典的主题样式  
    theme_classic() +
    # 设置y轴标签为genName变量的值  
    ylab(variable) +
    # 设置x轴标签为CancerName变量的值
    xlab("") +
    # 设置y轴的范围，基于tmp数据框的第1列的最小值和最大值的缩放  
    ylim(c(min(data[,variable])*0.9, max(data[,variable])*1.42)) +
    annotate(geom="text", 
             cex=6, 
             x=1.5, 
             y=max(data[,variable])*1.35,
             label=paste0("P ", ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3)))), color="black") + 
    # 设置图表的样式，包括面板边框、轴刻度线、轴标题和轴文本的样式 
    theme(
      # element_rect(): 用于定义矩形的元素。colour = "black": 设置面板边框的颜色为黑色。fill = NA: 设置面板的填充颜色为无即透明。size = 0.4: 设置边框线条的宽度为0.4个单位（这通常是磅或毫米，取决于具体的图形输出设备）。
      panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
      # axis.ticks = element_line(): 用于定义线条的元素。size = 0.2: 设置坐标轴刻度线的宽度为0.2个单位。color = "black": 设置坐标轴刻度线的颜色为黑色。
      axis.ticks = element_line(size=0.2, color="black"),
      # axis.ticks.length设置坐标轴刻度线的长度。这里使用了unit()函数来明确指定单位，这里是厘米（cm）。刻度线的长度被设置为0.2厘米。
      axis.ticks.length = unit(0.2,"cm"),
      # legend.position设置图例位置。"none": 不显示图例。在有些图中，尤其是当使用了多个geom_...()图层并且这些图层有不同的aes(fill = ...)或aes(color = ...)设置时，会自动生成图例。通过设置legend.position = "none可以隐藏这个图例。
      legend.position = "none",
      # axis.title设置坐标轴标签。element_text(): 用于定义文本的元素。size = 12: 设置坐标轴标签的字体大小为12
      axis.title = element_text(size = 12),
      #axis.text设置坐标轴刻度线标签。element_text(): 用于定义文本的元素。size = 10: 设置坐标轴刻度标签的字体大小为10
      axis.text = element_text(size = 10)) 
  ggsave(paste0(variable,"_",DatasetName,"_ggplot2_wilcox.pdf"), width = 4.5, height = 4)
  
}
mycolor = c("#BC3C29FF","#0072B5FF","#E18727FF",
            "#20854EFF","#7876B1FF","#6F99ADFF", 
            "#FFDC91FF","#EE4C97FF","#E64B35FF",
            "#4DBBD5FF","#00A087FF","#3C5488FF",
            "#F39B7FFF","#8491B4FF","#91D1C2FF",
            "#DC0000FF","#7E6148FF","#B09C85FF",
            "#3B4992FF","#EE0000FF","#008B45FF",
            "#631879FF","#008280FF","#BB0021FF",
            "#5F559BFF","#A20056FF","#808180FF",
            "#00468BFF","#ED0000FF","#42B540FF", 
            "#0099B4FF","#925E9FFF","#FDAF91FF",
            "#AD002AFF","#ADB6B6FF","#374E55FF",
            "#DF8F44FF","#00A1D5FF","#B24745FF",
            "#79AF97FF","#6A6599FF","#80796BFF",
            "#1f77b4",  "#ff7f0e",  "#279e68",
            "#d62728",  "#aa40fc",  "#8c564b",
            "#e377c2",  "#b5bd61",  "#17becf",
            "#aec7e8")
getkruskalggplot2 <- function(data,
                              Type = "Type",
                              variable,
                              levels = NULL,
                              mycolor = c("#0072B5FF","#BC3C29FF","#E18727FF",
                                          "#20854EFF","#7876B1FF","#6F99ADFF", 
                                          "#FFDC91FF","#EE4C97FF","#E64B35FF",
                                          "#4DBBD5FF","#00A087FF","#3C5488FF",
                                          "#F39B7FFF","#8491B4FF","#91D1C2FF",
                                          "#DC0000FF","#7E6148FF","#B09C85FF",
                                          "#3B4992FF","#EE0000FF","#008B45FF",
                                          "#631879FF","#008280FF","#BB0021FF",
                                          "#5F559BFF","#A20056FF","#808180FF",
                                          "#00468BFF","#ED0000FF","#42B540FF", 
                                          "#0099B4FF","#925E9FFF","#FDAF91FF",
                                          "#AD002AFF","#ADB6B6FF","#374E55FF",
                                          "#DF8F44FF","#00A1D5FF","#B24745FF",
                                          "#79AF97FF","#6A6599FF","#80796BFF",
                                          "#1f77b4",  "#ff7f0e",  "#279e68",
                                          "#d62728",  "#aa40fc",  "#8c564b",
                                          "#e377c2",  "#b5bd61",  "#17becf",
                                          "#aec7e8"),
                              DatasetName,
                              width = 4.5){
  suppressPackageStartupMessages(library(ggplot2,quietly = TRUE))
  suppressPackageStartupMessages(library(ggpubr,quietly = TRUE))
  suppressPackageStartupMessages(library(limma,quietly = TRUE))
  suppressPackageStartupMessages(library(ggExtra,quietly = TRUE))
  suppressPackageStartupMessages(library(reshape2,quietly = TRUE))
  suppressPackageStartupMessages(library(data.table,quietly = TRUE))
  suppressPackageStartupMessages(library(aplot,quietly = TRUE))
  suppressPackageStartupMessages(library(dplyr,quietly = TRUE))
  suppressPackageStartupMessages(library(hrbrthemes,quietly = TRUE))
  suppressPackageStartupMessages(library(ggtext,quietly = TRUE))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  p <- kruskal.test(expression ~ Type, data = data)$p.value
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  # 使用ggplot2函数进行绘图，数据来源于tmp数据框，x轴为Grade列，y轴为tmp数据框中由genName变量指定的列，并且使用Grade列作为填充色
  ggplot(data = data,
         aes(x = Type, 
             y = expression, 
             fill = Type)) +
    # 设置填充颜色
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    scale_color_manual(values = mycolor)+
    scale_fill_manual(values = mycolor) + 
    # 绘制小提琴图，用于显示数据的分布，设置透明度、位置、线条宽度和颜色 
    geom_violin(alpha=0.4, position = position_dodge(width = .75), size=0.8, color="black") + 
    # 在小提琴图上添加箱线图，显示中位数、四分位数等，设置缺口、异常点大小、线条宽度、透明度和颜色  
    geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.5) +
    # 在箱线图的基础上添加点图，使用21号形状（带边框的圆形），设置大小、位置、颜色和透明度  
    geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha=0.3) +
    # 使用经典的主题样式  
    theme_classic() +
    # 设置y轴标签为genName变量的值  
    ylab(variable) +
    # 设置x轴标签为CancerName变量的值
    xlab("") +
    # 设置y轴的范围，基于tmp数据框的第1列的最小值和最大值的缩放  
    ylim(c(min(data[,variable])*0.9, max(data[,variable])*1.42)) +
    annotate(geom="text", 
             cex=6, 
             x=(length(unique(data[,"Type"]))/2)+0.5, 
             y=max(data[,variable])*1.35,
             label=paste0("P ", ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3)))), color="black") + 
    # 设置图表的样式，包括面板边框、轴刻度线、轴标题和轴文本的样式 
    theme(
      # element_rect(): 用于定义矩形的元素。colour = "black": 设置面板边框的颜色为黑色。fill = NA: 设置面板的填充颜色为无即透明。size = 0.4: 设置边框线条的宽度为0.4个单位（这通常是磅或毫米，取决于具体的图形输出设备）。
      panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
      # axis.ticks = element_line(): 用于定义线条的元素。size = 0.2: 设置坐标轴刻度线的宽度为0.2个单位。color = "black": 设置坐标轴刻度线的颜色为黑色。
      axis.ticks = element_line(size=0.2, color="black"),
      # axis.ticks.length设置坐标轴刻度线的长度。这里使用了unit()函数来明确指定单位，这里是厘米（cm）。刻度线的长度被设置为0.2厘米。
      axis.ticks.length = unit(0.2,"cm"),
      # legend.position设置图例位置。"none": 不显示图例。在有些图中，尤其是当使用了多个geom_...()图层并且这些图层有不同的aes(fill = ...)或aes(color = ...)设置时，会自动生成图例。通过设置legend.position = "none可以隐藏这个图例。
      legend.position = "none",
      # axis.title设置坐标轴标签。element_text(): 用于定义文本的元素。size = 12: 设置坐标轴标签的字体大小为12
      axis.title = element_text(size = 12),
      #axis.text设置坐标轴刻度线标签。element_text(): 用于定义文本的元素。size = 10: 设置坐标轴刻度标签的字体大小为10
      axis.text = element_text(size = 10)) 
  ggsave(paste0(variable,"_",DatasetName,"_ggplot2_kruskal.pdf"), width = width, height = 4)
}

getkruskalggplot2_next <- function(data,
                                   Type = "Type",
                                   variable,
                                   levels = NULL,
                                   mycolor = c("#BC3C29FF","#0072B5FF","#E18727FF",
                                               "#20854EFF","#7876B1FF","#6F99ADFF", 
                                               "#FFDC91FF","#EE4C97FF","#E64B35FF",
                                               "#4DBBD5FF","#00A087FF","#3C5488FF",
                                               "#F39B7FFF","#8491B4FF","#91D1C2FF",
                                               "#DC0000FF","#7E6148FF","#B09C85FF",
                                               "#3B4992FF","#EE0000FF","#008B45FF",
                                               "#631879FF","#008280FF","#BB0021FF",
                                               "#5F559BFF","#A20056FF","#808180FF",
                                               "#00468BFF","#ED0000FF","#42B540FF", 
                                               "#0099B4FF","#925E9FFF","#FDAF91FF",
                                               "#AD002AFF","#ADB6B6FF","#374E55FF",
                                               "#DF8F44FF","#00A1D5FF","#B24745FF",
                                               "#79AF97FF","#6A6599FF","#80796BFF",
                                               "#1f77b4",  "#ff7f0e",  "#279e68",
                                               "#d62728",  "#aa40fc",  "#8c564b",
                                               "#e377c2",  "#b5bd61",  "#17becf",
                                               "#aec7e8"),
                                   width = 6,
                                   height = 4,
                                   DatasetName,
                                   title = ""){
  suppressPackageStartupMessages(library(ggplot2,quietly = TRUE))
  suppressPackageStartupMessages(library(ggpubr,quietly = TRUE))
  suppressPackageStartupMessages(library(limma,quietly = TRUE))
  suppressPackageStartupMessages(library(ggExtra,quietly = TRUE))
  suppressPackageStartupMessages(library(reshape2,quietly = TRUE))
  suppressPackageStartupMessages(library(data.table,quietly = TRUE))
  suppressPackageStartupMessages(library(aplot,quietly = TRUE))
  suppressPackageStartupMessages(library(dplyr,quietly = TRUE))
  suppressPackageStartupMessages(library(hrbrthemes,quietly = TRUE))
  suppressPackageStartupMessages(library(ggtext,quietly = TRUE))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  p <- kruskal.test(expression ~ Type, data = data)$p.value
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  pdf(paste0(variable,"_",DatasetName,"_Cloud_Rain_Kruskal.pdf"), width = width, height = height)
  plot01 = ggplot(data, aes(x = Type, 
                            y =data[,variable])) + 
    ggdist::stat_halfeye(aes(color=Type,
                             fill=Type),
                         adjust = .5, 
                         width = .7, 
                         .width = 0, 
                         justification = -.2, 
                         point_colour = NA) + 
    geom_boxplot(aes(color = Type),width = .2, outlier.shape = NA) + 
    geom_jitter(aes(color = Type),width = .05, alpha = .3) +
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    scale_color_manual(values = mycolor)+
    scale_fill_manual(values = mycolor) + 
    ylab(variable) +
    xlab("")+
    coord_flip()+
    labs(title = title,
         subtitle = paste0("Kruskal-Wallis Rank Sum Test ","P Value ",ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3))))
         #caption = "Visualization by <span style='color:#0057FF'>DataCharm</span>"
    ) +
    #hrbrthemes::theme_ipsum(base_family = "Arial Narrow") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
          plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                        size = 20, margin = margin(t = 1, b = 12)),
          plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
          plot.caption = element_markdown(face = 'bold',size = 12),
          legend.position = "none")
  print(plot01)
  dev.off()
}
getwilcoxonggplot2_next <- function(data,
                                    Type = "Type",
                                    variable,
                                    levels = NULL,
                                    mycolor = c("#BC3C29FF","#0072B5FF","#E18727FF",
                                                "#20854EFF","#7876B1FF","#6F99ADFF", 
                                                "#FFDC91FF","#EE4C97FF","#E64B35FF",
                                                "#4DBBD5FF","#00A087FF","#3C5488FF",
                                                "#F39B7FFF","#8491B4FF","#91D1C2FF",
                                                "#DC0000FF","#7E6148FF","#B09C85FF",
                                                "#3B4992FF","#EE0000FF","#008B45FF",
                                                "#631879FF","#008280FF","#BB0021FF",
                                                "#5F559BFF","#A20056FF","#808180FF",
                                                "#00468BFF","#ED0000FF","#42B540FF", 
                                                "#0099B4FF","#925E9FFF","#FDAF91FF",
                                                "#AD002AFF","#ADB6B6FF","#374E55FF",
                                                "#DF8F44FF","#00A1D5FF","#B24745FF",
                                                "#79AF97FF","#6A6599FF","#80796BFF",
                                                "#1f77b4",  "#ff7f0e",  "#279e68",
                                                "#d62728",  "#aa40fc",  "#8c564b",
                                                "#e377c2",  "#b5bd61",  "#17becf",
                                                "#aec7e8"),
                                    width = 6,
                                    height = 4,
                                    title = ""){
  suppressPackageStartupMessages(library(ggplot2,quietly = TRUE))
  suppressPackageStartupMessages(library(ggpubr,quietly = TRUE))
  suppressPackageStartupMessages(library(limma,quietly = TRUE))
  suppressPackageStartupMessages(library(ggExtra,quietly = TRUE))
  suppressPackageStartupMessages(library(reshape2,quietly = TRUE))
  suppressPackageStartupMessages(library(data.table,quietly = TRUE))
  suppressPackageStartupMessages(library(aplot,quietly = TRUE))
  suppressPackageStartupMessages(library(dplyr,quietly = TRUE))
  suppressPackageStartupMessages(library(hrbrthemes,quietly = TRUE))
  suppressPackageStartupMessages(library(ggtext,quietly = TRUE))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  p <- wilcox.test(expression ~ Type, data = data)$p.value
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  pdf(paste0(variable,"_","Cloud_Rain_Wilcoxon.pdf"), width = width, height = height)
  plot01 = ggplot(data, aes(x = Type, 
                            y =data[,variable])) + 
    ggdist::stat_halfeye(aes(color=Type,
                             fill=Type),
                         adjust = .5, 
                         width = .7, 
                         .width = 0, 
                         justification = -.2, 
                         point_colour = NA) + 
    geom_boxplot(aes(color = Type),width = .2, outlier.shape = NA) + 
    geom_jitter(aes(color = Type),width = .05, alpha = .3) +
    #ggsci::scale_color_nejm()+
    #ggsci::scale_fill_nejm() +
    scale_color_manual(values = mycolor)+
    scale_fill_manual(values = mycolor) + 
    ylab(variable) +
    xlab("")+
    coord_flip()+
    labs(title = title,
         subtitle = paste0("Wilcoxon Rank Sum Test ","P Value ",ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3))))
         #caption = "Visualization by <span style='color:#0057FF'>DataCharm</span>"
    ) +
    #hrbrthemes::theme_ipsum(base_family = "Arial Narrow") +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
          plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                        size = 20, margin = margin(t = 1, b = 12)),
          plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
          plot.caption = element_markdown(face = 'bold',size = 12),
          legend.position = "none")
  print(plot01)
  dev.off()
}

getwilcoxonggplot2_another <- function(data,
                                       Type = "Type",
                                       variable,
                                       levels = NULL,
                                       CancerName ,
                                       width = 6,
                                       height = 3){
  suppressPackageStartupMessages(library(ggplot2,quietly = TRUE))
  suppressPackageStartupMessages(library(ggpubr,quietly = TRUE))
  suppressPackageStartupMessages(library(limma,quietly = TRUE))
  suppressPackageStartupMessages(library(ggExtra,quietly = TRUE))
  suppressPackageStartupMessages(library(reshape2,quietly = TRUE))
  suppressPackageStartupMessages(library(data.table,quietly = TRUE))
  suppressPackageStartupMessages(library(aplot,quietly = TRUE))
  suppressPackageStartupMessages(library(dplyr,quietly = TRUE))
  data[,variable] = unlist(as.numeric(data[,variable]))
  data[,"expression"] = data[,variable]
  colnames(data)[colnames(data) == Type] = "Type"
  if(is.null(levels)){
    data[,"Type"] = factor(data[,"Type"], levels = unique(data[,"Type"])) 
  }else{
    data[,"Type"] = factor(data[,"Type"], levels = levels) 
  }
  p.val <- wilcox.test(expression ~ Type,data = data)
  p.lab <- paste0("P",ifelse(p.val$p.value < 0.001, " < 0.001",paste0(" = ",round(p.val$p.value, 3))))
  p.lab
  green <- alpha("#4DBBD5FF",0.5)
  cyan <- alpha("#E64B35FF",0.5)
  blue <- alpha("#4DBBD5FF",0.5)
  p_top <- ggplot(data, 
                  aes(x = expression, 
                      color = Type, 
                      fill = Type)
  ) +
    geom_density() +
    scale_color_manual(values = c(green,cyan,blue)) + 
    scale_fill_manual(values = c(green,cyan,blue)) +
    theme_classic() + 
    xlab(paste0("Estimated Expression of ", variable)) + 
    ylab(NULL) + 
    theme(legend.position = "none", 
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    geom_rug()
  p_top
  p_bot <- ggplot(data, aes(Type, expression, fill = Type)) + 
    geom_boxplot(aes(col = Type)) + 
    scale_fill_manual(values = c(green, cyan, blue)) + 
    scale_color_manual(values = c(green, cyan, blue)) + 
    xlab(NULL) + 
    ylab("Estimated Expression") + 
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 11,color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) + 
    annotate(geom="text",x = 0.2,
             hjust = 1,
             y = max(data[,"expression"])*1,
             size = 4, 
             angle = 270, 
             fontface = "bold",
             label = p.lab) +
    coord_flip()
  dat <- ggplot_build(p_bot)$data[[1]]
  p_bot <- p_bot + 
    geom_segment(data = dat, aes(x=xmin, 
                                 xend=xmax, 
                                 y=middle, 
                                 yend=middle), 
                 color="white", 
                 inherit.aes = F)
  p_bot
  p <- p_top %>% insert_bottom(p_bot, height = 0.4)
  pdf(paste0(CancerName,"_",variable,"_wilcoxon.pdf"), width = width,height = height)
  print(p)
  invisible(dev.off())
}

getgenepaired <- function(CancerName,GeneName){
  data = getRNA(CancerName)
  data = data[,which(substr(colnames(data),14,14)==0)]
  data = t(data)
  data_without_hubgene = data[, -grep(paste0("^",GeneName,"$"), colnames(data))]  
  n_cols = ncol(data_without_hubgene)  
  data_01 = matrix(0, nrow = nrow(data), ncol = n_cols)  
  for (i in seq_len(n_cols)) {  
    data_01[, i] = ifelse(data[, GeneName] > data_without_hubgene[, i], 1, 0)  
  }  
  colnames(data_01) = colnames(data_without_hubgene)  
  data_01 = as.data.frame(data_01)  
  rownames(data_01) = rownames(data_without_hubgene)  
  prop_ones = colMeans(data_01 == 1)  
  keep_cols = which(prop_ones <= 0.9)  
  data_01_filter = data_01[, keep_cols]  
  colnames(data_01_filter) = paste0(GeneName,"|",colnames(data_01_filter))
  return(data_01_filter)
}
getscRNAGeneExp <- function(scRNA,GeneName){
  expression_data <- GetAssayData(object = scRNA, assay = "RNA") %>% .[c(GeneName),] %>% as.matrix()
  if(length(GeneName)==1){
    colnames(expression_data) = GeneName
    expression_data = as.data.frame(expression_data)
  }else{
    expression_data = t(expression_data)
    expression_data = as.data.frame(expression_data)
  }
  return(expression_data)
}
getstRNAGeneExp <- function(scRNA,GeneName){
  expression_data <- GetAssayData(object = scRNA, assay = "Spatial") %>% .[c(GeneName),] %>% as.matrix()
  if(length(GeneName)==1){
    colnames(expression_data) = GeneName
    expression_data = as.data.frame(expression_data)
  }else{
    expression_data = t(expression_data)
    expression_data = as.data.frame(expression_data)
  }
  return(expression_data)
}
getscRNAVlnPlot <- function(scRNA,GeneName,width=4.5,height=4){
  library(Seurat)
  library(ggplot2)
  pdf(paste0(GeneName,"_VlnPlot.pdf"), width = width, height = height)
  plot = VlnPlot(scRNA, features = GeneName, pt.size = 0)+ 
    scale_color_manual(values = alpha(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF",
                                        "#1f77b4", "#ff7f0e", "#279e68", "#d62728","#aa40fc", "#8c564b", "#e377c2", "#b5bd61","#17becf","#aec7e8")
                                      ,0.5)) +
    scale_fill_manual(values = alpha(c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF",
                                       "#1f77b4", "#ff7f0e", "#279e68", "#d62728","#aa40fc", "#8c564b", "#e377c2", "#b5bd61","#17becf","#aec7e8")
                                     ,0.5))+
    scale_y_continuous(expand = c(0, 0))+
    xlab("Chemotherapy") +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=0.4), 
      axis.ticks = element_line(size=0.2, color="black"),
      axis.ticks.length = unit(0.2,"cm"),
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10))
  print(plot)
  dev.off()
}


getQ1Q2Q3Q4 <- function(data){
  getQuantile <- function(x){
    ifelse(x>quantile(x,.75),"Q1",
           ifelse(x>quantile(x,.5),"Q2",
                  ifelse(x>quantile(x,.25),"Q3",
                         "Q4")
           )
    )
  }
  new_quantile = apply(data[,1,drop=F], 2, getQuantile)
  colnames(new_quantile) = "Quantile"
  ExpressionType = cbind(data,new_quantile)
  idx = order(ExpressionType$Quantile)#获得排序后索引
  ExpressionType = ExpressionType[idx,]#按索引重新提取以达到排序效果
  return(ExpressionType)
}

getstrsplit <- function(data,whatstr = "_",whatcount = 1){
  varible = sapply(strsplit(data,whatstr,fixed = TRUE) ,'[',whatcount)
  return(varible)
}

getfilenames <- function(path){
  filenames = list.files(path = path, full.names = FALSE)  
  return(filenames)
}

getfilenames_next_step <- function(whatstr = "RDS$",filenames){
  files = grep(whatstr,filenames,value=T) 
  return(files)
}

getkmplot <- function(data,minprop=0.3,GeneName,CancerName){
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(dplyr))
  bestcut = surv_cutpoint(data,
                          time = colnames(data)[1],
                          event = colnames(data)[2],
                          variables = GeneName ,
                          minprop = minprop) 
  SurvivalType = ifelse(colnames(data)[2] == "OS","Overall Survival",
                        ifelse(colnames(data)[2] == "DSS","Disease-specific Survival",
                               ifelse(colnames(data)[2] == "PFS","Progression-free Survival",
                                      ifelse(colnames(data)[2] == "DFS","Disease-free Survival",
                                             ifelse(colnames(data)[2] == "RFS","Relapse-free Survival",
                                                    ifelse(colnames(data)[2] == "DFI","Disease-free Interval",
                                                           ifelse(colnames(data)[2] == "PFI","Progression-free Interval","Survival")
                                                    )
                                             )
                                      )
                               )
                        )
  )
  cutoff = bestcut$cutpoint[1,1]
  data$group = factor(ifelse(data[,GeneName] > cutoff, "High","Low"), levels = c("High","Low"))
  data$group2 = factor(ifelse(data[,GeneName] > median(data[,GeneName]),"High","Low"), levels = c("High","Low"))
  if(substr(rownames(data)[1],1,3)=="GSM"){
    rownames(data) = paste0(CancerName,"_",substr(gsub("GSM","",rownames(data)),nchar(gsub("GSM","",rownames(data)))-3,nchar(gsub("GSM","",rownames(data)))),"_Rename")
    data2 = data[order(data[,GeneName],decreasing = T),,drop = F]
    write.csv(data2,paste0(GeneName,"_",CancerName,"_",SurvivalType,"_KM.csv"))
  }else{
    write.csv(data,paste0(GeneName,"_",CancerName,"_",SurvivalType,"_KM.csv"))
  }
  diff = survdiff(Surv(data[,1], data[,2]) ~ group,data = data, na.action = na.exclude)
  pValue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  if(pValue < 0.001){
    pValue = "< 0.001"
  }else{
    pValue = paste0("= ",sprintf("%.03f",pValue))
  }
  fit = survfit(Surv(data[,1], data[,2]) ~ group, data = data)
  pdf(paste0(GeneName,"_",CancerName,"_",SurvivalType,"_The_Best_Cutoff.pdf"),onefile = FALSE,width = 5,height =5)
  plot(fit,
       lty = 2:3,
       lwd = 2,
       col = c("#E64B35FF", "#4DBBD5FF"),
       xlab = "Time (year)",
       ylab = paste0("Survival Rate (",SurvivalType,")"),
       main = paste("log-rank test p", pValue , "(",CancerName,")"),
       mark.time=T)
  legend("bottomright",
         lty = 2:3,
         c(paste0("High ", GeneName," (n= ",nrow(data[data$group=="High",]),")"), 
           paste0("Low  ", GeneName," (n= ",nrow(data[data$group=="Low",]), ")")),
         lwd=2,
         col=c("#E64B35FF", "#4DBBD5FF"))
  dev.off()
  
  diff2 = survdiff(Surv(data[,1], data[,2]) ~ group2, data = data, na.action = na.exclude)
  pValue2 = 1 - pchisq(diff2$chisq, length(diff2$n) - 1)
  if(pValue2<0.001){
    pValue2 = "< 0.001"
  }else{
    pValue2 = paste0("= ",sprintf("%.03f",pValue2))
  }
  fit2 <- survfit(Surv(data[,1], data[,2]) ~ group2, data = data)
  pdf(paste0(GeneName,"_",CancerName,"_",SurvivalType,"_The_median_cutoff.pdf"),onefile = FALSE,width = 5,height =5)
  plot(fit2,
       lty = 2:3,
       lwd = 2,
       col = c("#E64B35FF", "#4DBBD5FF"),
       xlab = "Time (year)",
       ylab = paste0("Survival Rate (",SurvivalType,")"),
       main = paste("log-rank test p", pValue2 , "(",CancerName,")"),
       mark.time=T)
  legend("bottomright",
         lty=2:3,
         c(paste0("High ", GeneName," (n= ",nrow(data[data$group2 == "High",]),")"), 
           paste0("Low  ", GeneName," (n= ",nrow(data[data$group2 == "Low",]), ")")),
         lwd=2,
         col=c("#E64B35FF", "#4DBBD5FF"))
  dev.off()
}

getkmplot_OS_Month <- function(data,group = "Subtype",variable){
  data = data[,c("OS.time","OS",group)]
  data$OS.time = data$OS.time/30.4
  data = as.data.frame(data)
  data = data[data$OS.time>0,]
  data = data[!(is.na(data$OS.time)),]
  data = data[!(is.na(data$OS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(OS.time, OS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(OS.time, OS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(OS.time, OS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,
                 legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Overall survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_OS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}

getkmplot_OS_Month_external <- function(data,group = "Subtype",variable){
  data = data[,c("OS.time","OS",group)]
  data$OS.time = data$OS.time*12
  data = as.data.frame(data)
  data = data[data$OS.time>0,]
  data = data[!(is.na(data$OS.time)),]
  data = data[!(is.na(data$OS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(OS.time, OS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(OS.time, OS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(OS.time, OS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,
                 legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Overall survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_OS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_DSS_Month <- function(data, group = "Subtype",variable){
  data = data[,c("DSS.time","DSS",group)]
  data$DSS.time = data$DSS.time/30.4
  data = as.data.frame(data)
  data = data[data$DSS.time>0,]
  data = data[!(is.na(data$DSS.time)),]
  data = data[!(is.na(data$DSS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(DSS.time, DSS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(DSS.time, DSS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DSS.time, DSS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,
                 legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Disease-specific Survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_DSS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_DSS_Month_external <- function(data, group = "Subtype",variable){
  data = data[,c("DSS.time","DSS",group)]
  data$DSS.time = data$DSS.time*12
  data = as.data.frame(data)
  data = data[data$DSS.time>0,]
  data = data[!(is.na(data$DSS.time)),]
  data = data[!(is.na(data$DSS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(DSS.time, DSS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(DSS.time, DSS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DSS.time, DSS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,
                 legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Disease-specific Survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_DSS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_PFI_Month <- function(data,group = "Subtype",variable){
  data = data[,c("PFI.time","PFI",group)]
  data$PFI.time = data$PFI.time/30.4
  data = as.data.frame(data)
  data = data[data$PFI.time>0,]
  data = data[!(is.na(data$PFI.time)),]
  data = data[!(is.na(data$PFI)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(PFI.time, PFI) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(PFI.time, PFI)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(PFI.time, PFI)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Progression-free Interval",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_PFI.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_PFS_Month_external <- function(data, group = "Subtype",variable){
  data = data[,c("PFS.time","PFS",group)]
  data$PFS.time = data$PFS.time*12
  data = as.data.frame(data)
  data = data[data$PFS.time>0,]
  data = data[!(is.na(data$PFS.time)),]
  data = data[!(is.na(data$PFS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(PFS.time, PFS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(PFS.time, PFS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(PFS.time, PFS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Progression-free Survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_PFS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_DFI_Month <- function(data,group = "Subtype",variable){
  data = data[,c("DFI.time","DFI",group)]
  data$DFI.time = data$DFI.time/30.4
  data = as.data.frame(data)
  data = data[data$DFI.time>0,]
  data = data[!(is.na(data$DFI.time)),]
  data = data[!(is.na(data$DFI)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(DFI.time, DFI) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(DFI.time, DFI)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DFI.time, DFI)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Disease-free Interval",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_DFI.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_DFS_Month_external <- function(data,group = "Subtype",variable){
  data = data[,c("DFS.time","DFS",group)]
  data$DFS.time = data$DFS.time*12
  data = as.data.frame(data)
  data = data[data$DFS.time>0,]
  data = data[!(is.na(data$DFS.time)),]
  data = data[!(is.na(data$DFS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(DFS.time, DFS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(DFS.time, DFS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DFS.time, DFS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Disease-free Survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_DFS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_RFS_Month_external <- function(data, group = "Subtype",variable){
  data = data[,c("RFS.time","RFS",group)]
  data$RFS.time = data$RFS.time*12
  data = as.data.frame(data)
  data = data[data$RFS.time>0,]
  data = data[!(is.na(data$RFS.time)),]
  data = data[!(is.na(data$RFS)),]
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  data$group = data[,group]
  fitd = survdiff(Surv(RFS.time, RFS) ~ group, 
                  data = data, 
                  na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(RFS.time, RFS)~ group,
                data = data,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(Surv(RFS.time, RFS)~ group,
                         data = data,
                         p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH
  mycol = brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) = gsub("group=", "", names(fit$strata))
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p = ggsurvplot(fit = fit,
                 conf.int = FALSE,
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
                 data  = data,
                 xlim = c(0,120),
                 size = 1,
                 break.time.by = 12,legend.title = "",
                 xlab = "Time (months)",
                 ylab  = "Relapse-free Survival",
                 risk.table.y.text = FALSE,
                 tables.height = 0.3)
  p
  p.lab = paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot = p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab)
  p
  if(length(unique(data$group)>2)){
    ## 添加配对表格
    addTab = as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
    addTab[is.na(addTab)] = "-"
    df = tibble(x = 120, y = 1, tb = list(addTab))
    p$plot = p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
    p
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(variable,"_km_pairwise_RFS.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}
getkmplot_quantile <- function(data,GeneName,DatasetName){
  SurvivalType = ifelse(colnames(data)[2] == "OS","Overall Survival",
                        ifelse(colnames(data)[2] == "DSS","Disease-specific Survival",
                               ifelse(colnames(data)[2] == "PFS","Progression-free Survival",
                                      ifelse(colnames(data)[2] == "DFS","Disease-free Survival",
                                             ifelse(colnames(data)[2] == "RFS","Relapse-free Survival",
                                                    ifelse(colnames(data)[2] == "DFI","Disease-free Interval",
                                                           ifelse(colnames(data)[2] == "PFI","Progression-free Interval","Survival")
                                                    )
                                             )
                                      )
                               )
                        )
  )
  Quantile_data = data[,GeneName,drop=F]
  Quantile = function(x){
    ifelse(x>quantile(x,.75),"Q1",
           ifelse(x>quantile(x,.5),"Q2",
                  ifelse(x>quantile(x,.25),"Q3",
                         "Q4")
           )
    )
  }
  Quantile_data_new = apply(Quantile_data, 2, Quantile)
  colnames(Quantile_data_new) = "Expression_Type"
  Quantile_data_merge = cbind(Quantile_data,Quantile_data_new)
  idx = order(Quantile_data_merge$Expression_Type)
  Quantile_data_merge_order = Quantile_data_merge[idx,]
  SameSamples = intersect(rownames(data),rownames(Quantile_data_merge_order))
  data = data[SameSamples,,drop=F]
  Quantile_data_merge_order = Quantile_data_merge_order[SameSamples,,drop=F]
  dat = cbind(data,Quantile_data_merge_order)
  dat = as.data.frame(dat)
  write.csv(dat,paste0(GeneName,"_",DatasetName,"_",colnames(data)[2],"_plotdata_All_Survival_Data_Q1Q2Q3Q4.csv"),row.names=T)
  dat$group = dat$Expression_Type
  dat = dat[,c(colnames(data)[1],colnames(data)[2],"group")]
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(survminer))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  formula_name <- as.formula(paste0("Surv(",colnames(dat)[1],",",colnames(dat)[2],")","~","group"))
  fitd = survdiff(Surv(dat[,1],dat[,2]) ~ group,data = dat, na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit = survfit(Surv(dat[,1],dat[,2]) ~ group,
                data = dat,
                type= "kaplan-meier",
                error = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
  ps = pairwise_survdiff(formula_name,
                         data = dat,
                         p.adjust.method = "none")
  mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]
  names(fit$strata) <- gsub("group=", "", names(fit$strata))
  p.lab <- paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p <- ggsurvplot(fit = fit,
                  conf.int = FALSE,
                  risk.table = TRUE,
                  risk.table.col = "strata",
                  palette = c("#E64B35FF", "#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF"), # KM曲线颜色
                  data  = dat,
                  xlim = c(0,120),
                  size = 1,
                  break.time.by = 12,
                  legend.title = "",
                  xlab = "Time (months)",
                  ylab  = SurvivalType,
                  risk.table.y.text = FALSE,
                  tables.height = 0.3)
  p
  
  p.lab <- paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3))))
  p$plot <- p$plot + annotate("text",x = 0, y = 0.1, hjust = 0,fontface = 4,label = p.lab);p
  ## 添加配对表格
  addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3))))
  addTab[is.na(addTab)] <- "-"
  df <- tibble(x = 120, y = 1, tb = list(addTab))
  p$plot <- p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
  p
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(paste0(GeneName,"_",DatasetName,"_",colnames(data)[2],"_km_pairwise_Q1Q2Q3Q4.pdf"), width = 6.5, height = 6)
  print(p)
  dev.off()
}
getlimmadegs_gene <- function(data,GeneName){
  library(limma)
  GeneData = data[GeneName,,drop=F] %>% t() %>% .[order(.[,1],decreasing = T),,drop = F]
  GeneData_high_samples = rownames(GeneData[1:(nrow(GeneData) * pct),,drop = F])
  GeneData_low_samples  = rownames(GeneData[nrow(GeneData):(nrow(GeneData) - nrow(GeneData) * pct + 1),,drop = F])
  data_high_samples = data[,GeneData_high_samples]
  case_number = ncol(data_high_samples)
  data_low_samples = data[,GeneData_low_samples]
  control_number = ncol(data_low_samples)
  data_high_low = cbind(data_high_samples,data_low_samples)
  Type =  factor(c(rep("case",case_number),rep("control",control_number)))
  design = model.matrix(~0 + factor(Type))
  colnames(design) = levels(Type)
  rownames(design) = colnames(data_high_low)
  contrast.matrix = makeContrasts(case-control, levels = design)
  fit = lmFit(data_high_low, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  limmma_results = topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  write.csv(limmma_results, paste0("easy_input_limma_",GeneName,"_",CancerName,"_",pctName,".csv"), quote = F)
  return(limmma_results)
}

getlimmadegs_gene_new <- function(data,GeneName,pct = pct,DatasetName){
  library(limma)
  library(tidyverse)
  GeneData = data[GeneName,,drop=F] %>% t() %>% .[order(.[,1],decreasing = T),,drop = F]
  GeneData_high_samples = rownames(GeneData[1:(nrow(GeneData) * pct),,drop = F])
  GeneData_low_samples  = rownames(GeneData[nrow(GeneData):(nrow(GeneData) - nrow(GeneData) * pct + 1),,drop = F])
  data_high_samples = data[,GeneData_high_samples]
  case_number = ncol(data_high_samples)
  data_low_samples = data[,GeneData_low_samples]
  control_number = ncol(data_low_samples)
  data_high_low = cbind(data_high_samples,data_low_samples)
  Type =  factor(c(rep("case",case_number),rep("control",control_number)))
  design = model.matrix(~0 + factor(Type))
  colnames(design) = levels(Type)
  rownames(design) = colnames(data_high_low)
  contrast.matrix = makeContrasts(case-control, levels = design)
  fit = lmFit(data_high_low, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  limmma_results = topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  write.csv(limmma_results, paste0("easy_input_limma_",GeneName,"_",CancerName,"_",DatasetName,"_",pct,".csv"), quote = F)
  return(limmma_results)
}

getlimmadegs_gene_fgsea <- function(data,GeneName,DatasetName,pct,gmt=NULL,gmtrds){
  GeneData = data[GeneName,,drop=F] %>% t() %>% .[order(.[,1],decreasing = T),,drop = F]
  GeneData_high_samples = rownames(GeneData[1:(nrow(GeneData) * pct),,drop = F])
  GeneData_low_samples  = rownames(GeneData[nrow(GeneData):(nrow(GeneData) - nrow(GeneData) * pct + 1),,drop = F])
  data_high_samples = data[,GeneData_high_samples]
  case_number = ncol(data_high_samples)
  data_low_samples = data[,GeneData_low_samples]
  control_number = ncol(data_low_samples)
  data_high_low = cbind(data_high_samples,data_low_samples)
  Type =  factor(c(rep("case",case_number),rep("control",control_number)))
  design = model.matrix(~0 + factor(Type))
  colnames(design) = levels(Type)
  rownames(design) = colnames(data_high_low)
  contrast.matrix = makeContrasts(case-control, levels = design)
  fit = lmFit(data_high_low, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  limmma_results = topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  write.csv(limmma_results, paste0("easy_input_limma_",GeneName,"_",DatasetName,"_",pct,".csv"), quote = F)
  limmma_results$logFC = as.numeric(as.character(limmma_results$logFC))
  limmma_results_order = limmma_results[order(limmma_results$logFC, decreasing = T), ]
  GSEA_input = limmma_results_order$logFC
  names(GSEA_input) = rownames(limmma_results_order)
  if(!is.null(gmt)){
    # 自定义函数读取GMT文件  
    readGMT <- function(gmtFilePath) {  
      # 初始化一个空列表  
      gmtList <- list()  
      
      # 读取GMT文件，跳过前两行（通常是注释）  
      gmtFile <- read.table(gmtFilePath, header = FALSE, sep = "\t", skip = 2, quote = "", stringsAsFactors = FALSE)  
      
      # 遍历GMT文件的每一行  
      for (i in 1:nrow(gmtFile)) {  
        # 获取基因集名称  
        geneSetName <- gmtFile[i, 1]  
        
        # 获取基因集成员，并去除可能的空格  
        geneSetMembers <- gmtFile[i, -c(1:2)]  # 假设第二列是描述信息，我们不需要它  
        geneSetMembers <- geneSetMembers[geneSetMembers != ""]  # 去除空字符串  
        
        # 将基因集及其成员添加到列表中  
        gmtList[[geneSetName]] <- geneSetMembers  
      }  
      
      return(gmtList)  
    }  
    gmtList = readGMT(gmt)
  }else{
    gmtList = readRDS(gmtrds)
  }
  fgseaRes = fgsea(pathways = gmtList, 
                   stats = GSEA_input,
                   minSize = 1,
                   maxSize = 500,
                   nperm = 10000)
  fgseaRes = fgseaRes[,1:7]
  return(fgseaRes)
}

getlimmadegs_gene_GSEA <- function(data,GeneName,DatasetName,pct,gmt=NULL,gmtrds){
  library(limma)
  library(tidyverse)
  GeneData = data[GeneName,,drop=F] %>% t() %>% .[order(.[,1],decreasing = T),,drop = F]
  GeneData_high_samples = rownames(GeneData[1:(nrow(GeneData) * pct),,drop = F])
  GeneData_low_samples  = rownames(GeneData[nrow(GeneData):(nrow(GeneData) - nrow(GeneData) * pct + 1),,drop = F])
  data_high_samples = data[,GeneData_high_samples]
  case_number = ncol(data_high_samples)
  data_low_samples = data[,GeneData_low_samples]
  control_number = ncol(data_low_samples)
  data_high_low = cbind(data_high_samples,data_low_samples)
  Type =  factor(c(rep("case",case_number),rep("control",control_number)))
  design = model.matrix(~0 + factor(Type))
  colnames(design) = levels(Type)
  rownames(design) = colnames(data_high_low)
  contrast.matrix = makeContrasts(case-control, levels = design)
  fit = lmFit(data_high_low, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  limmma_results = topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
  limmma_results = na.omit(limmma_results)
  write.csv(limmma_results, paste0("easy_input_limma_",GeneName,"_",DatasetName,"_",pct,".csv"), quote = F)
  limmma_results$logFC = as.numeric(as.character(limmma_results$logFC))
  limmma_results_order = limmma_results[order(limmma_results$logFC, decreasing = T), ]
  GSEA_input = limmma_results_order$logFC
  names(GSEA_input) = rownames(limmma_results_order)
  if(!is.null(gmt)){
    gmtList = read.gmt(gmt)
  }else{
    gmtList = readRDS(gmtrds)
  }
  suppressPackageStartupMessages(library(clusterProfiler))
  gseaRes = GSEA(geneList = GSEA_input,
                 minGSSize = 1,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 seed = TRUE,
                 TERM2GENE = gmtList)
  return(gseaRes)
}
getlimmadegs_gene_GSEA_customized <- function(data,GeneName,CancerName,DatasetName){
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(DOSE))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(ggrepel))
  limmma_results = getlimmadegs_gene_new(data = data, GeneName = GeneName, pct = pct,DatasetName = DatasetName)
  limmma_results = na.omit(limmma_results)
  limmma_results$logFC = as.numeric(as.character(limmma_results$logFC))
  limmma_results_order = limmma_results[order(limmma_results$logFC, decreasing = T), ]
  GSEA_input = limmma_results_order$logFC
  names(GSEA_input) = rownames(limmma_results_order)
  kk = GSEA(geneList = GSEA_input , 
            TERM2GENE = customized, 
            verbose = F,
            minGSSize = 2, 
            maxGSSize = 500, 
            nPerm = 10000, 
            pvalueCutoff = 1
  )
  write.csv(kk, paste0(GeneName,"_",CancerName,"_",DatasetName,"_GSEA.csv"))
  plot = gseaplot2(kk, geneSetID=c(1:dim(kk)[1]), 
                   color = c("#E64B35FF", "#4DBBD5FF","#3C5488FF","#00A087FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")[1:dim(kk)[1]])
  addTab = as.data.frame(as.matrix(ifelse(round(kk@result[["pvalue"]], 3) < 0.001, "<0.001", paste0(round(kk@result[["pvalue"]], 3)))),
                         as.matrix(paste0(round(kk@result$NES,3))))
  addTab_next = data.frame(row.names(addTab), addTab$V1)
  colnames(addTab_next) = c("NES","pvalue")
  rownames(addTab_next) = "GSEA"
  df = tibble(x = 1, y = 1, tb = list(addTab_next))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(ggpp))
  plot[[1]] = plot[[1]] + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)
  ggsave(paste0(GeneName,"_",CancerName,"_",DatasetName,"_GSEA.pdf"), plot = plot ,width=5, height=6)
}
getPanGeneZscore_TCGA_GTEx <- function(GeneName){
  data <- read.csv(paste0("/mnt/TCGA_GTEx/",GeneName,"_TCGA_GTEx.csv"))
  return(data)
}

getcortest_scatter <- function(data,GeneName1,GeneName2,DatasetName,
                               method = 'pearson'){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(ggExtra))
  suppressPackageStartupMessages(library(limma))
  data$shape = ifelse(data[,GeneName1] > median(data[,GeneName1]) & data[,GeneName2]>median(data[,GeneName2]),"Mixed",
                      ifelse(data[,GeneName1] > median(data[,GeneName1]) & data[,GeneName2]<=median(data[,GeneName2]),GeneName1,
                             ifelse(data[,GeneName1] <= median(data[,GeneName1]) & data[,GeneName2]>median(data[,GeneName2]),GeneName2,
                                    "Quiescent")
                      )
  )
  write.csv(data,paste0(GeneName1,"_",GeneName2,"_",DatasetName,"_plotdata_Cor2genes.csv"),row.names=T)
  data$shape = factor(data$shape, levels=levels(factor(data$shape)))
  shape = levels(factor(data$shape))
  length = length(levels(factor(data$shape)))
  bioCol = c("#BC3C29FF","#0072B5FF","#7876B1FF","#6F99ADFF")
  library(ggplot2)  
  p1 = ggplot(data, aes(data[,GeneName1], data[,GeneName2])) + 
    xlab(GeneName1) +
    ylab(GeneName2) +
    geom_point(aes(colour=shape,shape=shape))+
    scale_color_manual(values=bioCol[1:length])+ 
    scale_shape_manual(values=c(2,8,12,18))+
    geom_smooth(method="lm",formula = y ~ x) +
    theme_bw()+
    stat_cor(method = method, aes(x =data[,GeneName1], y =data[,GeneName2]))+
    theme(legend.position = "none")#bottom
  p2 = ggMarginal(p1, type="densigram", xparams=list(fill = alpha("#E18727FF",0.5)), yparams=list(fill = alpha("#20854EFF",0.5)))
  print(p2)
  p2
  ggsave(paste0(GeneName1,"_",GeneName2,"_",DatasetName,"_Scatter_diagram.pdf"), p2,width = 5.5, height = 5.5)
}

getcortest_fisher <- function(data,GeneName1,GeneName2,DatasetName,
                              method = 'pearson'){
  simdata = data
  r.cor = cor(simdata[,GeneName1], simdata[,GeneName2],method = method)
  r.cor
  p.cor = cor.test(simdata[,GeneName1], simdata[,GeneName2],method = method)$p.value
  p.cor
  simdata$GeneName1_classify = factor(cut(simdata[,GeneName1], quantile(simdata[,GeneName1]), 
                                          labels = c("Negative", "Weak", "Moderate", "Positive")), 
                                      levels = c("Positive", "Moderate", "Weak", "Negative")) # 根据例文设置的字符串排序
  simdata$GeneName2_classify = factor(cut(simdata[,GeneName2], quantile(simdata[,GeneName2]), 
                                          labels = c("Negative", "Weak", "Moderate", "Positive")), 
                                      levels = c("Negative", "Weak", "Moderate", "Positive")) # 根据例文设置的字符串排序
  tab_classify = as.data.frame.array(table(simdata$GeneName1_classify,simdata$GeneName2_classify)) # 将列联表转换为矩阵
  tab_classify
  # 分类较多时fisher'S exact test会需要更大的工作空间以及较长的时耗，一般设为1e9
  p.fisher <- fisher.test(tab_classify, workspace = 3e9,simulate.p.value=TRUE)$p.value 
  p.fisher
  blue   <- "#204F8D"
  lblue  <- "#498EB9"
  dwhite <- "#B6D1E8"
  white  <- "#E6EAF7"
  pdf(paste0(GeneName1,"_",GeneName2,"_",DatasetName,"_Fisher’s_exact_test.pdf"),width = 6,height = 6)
  par(bty="n", mgp = c(2,0.5,0), mar = c(5.1,6.1,4.1,2.1),tcl=-.25, font.main=3)
  par(xpd=NA)
  # 生成一个空白图
  plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)), # 产生左下和左上两个点，张开画布
       col = "white", # 点设置为白色
       xlab = "",xaxt = "n", # 不显示x坐标轴
       ylab = "",yaxt = "n") # 不显示y坐标轴)
  title(paste("Correlation between ",GeneName1,"and ",GeneName2,"\nrho = ", round(r.cor,2), "; ", "P.cor = ", format(p.cor, digits = 3, scientific = T),
              "\nP.fisher = ", format(p.fisher, digits = 3, scientific = T)), adj = 0, line = 0) # 科学计数法保留3位,左对齐title
  # 画纵坐标
  axis(2, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE) # 生成y坐标轴刻度，并位于每个间隔的中心位置
  text(y = 0.5:(ncol(tab_classify)-0.5), # 生成坐标轴标签
       par("usr")[1], 
       labels = rownames(tab_classify)[nrow(tab_classify):1], # 注意这里根据例文原图设置为倒置的坐标
       srt = 0, pos = 2, xpd = TRUE)
  mtext(GeneName2, side=2, line = 4.5) # 生成坐标轴名
  # 画横坐标
  axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE) # 生成x坐标轴刻度，并位于每个间隔的中心位置
  text(x = 0.5:(ncol(tab_classify)-0.5), # 生成坐标轴标签
       par("usr")[1] - 0.2, # 微调横坐标标签与坐标轴的相对位置，防止重合
       labels = colnames(tab_classify), 
       srt = 45, pos = 1, xpd = TRUE)
  mtext(GeneName1, side=1, line = 3.5) # 生成坐标轴名
  
  # 产生线性变化的颜色
  input_matrix <- as.matrix(tab_classify) 
  mat.max = max(input_matrix) # 输入矩阵的最大计数值
  unq.value <- unique(sort(as.vector(input_matrix))) # 得出所有独特的计数值
  rbPal <- colorRampPalette(c(white,dwhite,lblue,blue)) # 产生颜色区间函数
  col.vec <- rbPal(max(unq.value) + 1) # 产生计数值个数+1个颜色（+1目的是防止计数有0值出现）
  col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix)) # 生成空的颜色矩阵
  # 根据计数矩阵填充颜色矩阵
  for (i in 1:nrow(col.mat)) {
    for (j in 1:ncol(col.mat)) {
      col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
    }
  }
  # 通过矩形块产生热图
  x_size <- ncol(input_matrix)
  y_size <- nrow(input_matrix)
  my_xleft = rep(c(0:(x_size-1)),each = x_size) # 产生各个矩形的左x点
  my_xright = my_xleft + 1 # 产生各个矩形的右x点
  my_ybottom = rep(c((y_size-1):0),y_size) # 产生各个矩形的下y点
  my_ytop = my_ybottom + 1 # 产生各个矩形的上y点
  rect(xleft = my_xleft,
       ybottom = my_ybottom,
       xright = my_xright,
       ytop = my_ytop,
       col=col.mat, # 填充颜色
       border = F) # 取消矩形的边界
  text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1.3) # 填充计数值
  invisible(dev.off())
}

get2varibles_mean <- function(GeneName1,GeneName2,CancerName){
  GeneName1_data = getPanGeneTime(GeneName1)
  GeneName1_data = GeneName1_data[(GeneName1_data[,"CancerType"] == CancerName),]
  GeneName1_data = GeneName1_data[,GeneName1,drop=F]
  GeneName2_data = getPanGeneTime(GeneName2)
  GeneName2_data = GeneName2_data[(GeneName2_data[,"CancerType"] == CancerName),]
  GeneName2_data = GeneName2_data[,GeneName2,drop=F]
  GeneName1_GeneName2_merge_data = getmerge_row(GeneName1_data,GeneName2_data)
  GeneName1_GeneName2_merge_data_scale = scale(GeneName1_GeneName2_merge_data, center = TRUE, scale = TRUE)
  GeneName1_GeneName2_merge_data_scale =as.data.frame(GeneName1_GeneName2_merge_data_scale)
  GeneName1_GeneName2_merge_data_scale$Subtype <- 
    ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] >= 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] >= 0,paste0(GeneName1,"+ & ",GeneName2,"+"),
           ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] <= 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] <= 0,paste0(GeneName1,"- & ",GeneName2,"-"),
                  ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] > 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] < 0,paste0(GeneName1,"+ & ",GeneName2,"-"),
                         paste0(GeneName1,"- & ",GeneName2,"+"))))
  write.csv(GeneName1_GeneName2_merge_data_scale,paste0(GeneName1,"_",GeneName2,"_merge_data_scale.csv"))
  
  #-- -+ +- ++
  palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF")
  pdf(file = paste0(GeneName1,"_",GeneName2,"_",CancerName,"_Mean_classification.pdf"),width = 5, height = 5)
  par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,4.1,2.1,3.1)+.1, las=1, tcl=-.25)
  plot(NULL, NULL, # 绘制空白背景
       xlim =c(-4,4),
       ylim =c(-4,4),
       xlab = paste0(GeneName1," (z-score)"),
       ylab = paste0(gsub("_"," ",GeneName2)," (z-score)"))
  grid(col = "grey85", lty = 2, lwd = 1.5) # 添加网格线
  abline(v = 0, lty = 2, lwd = 2) # 添加水平0截断
  abline(h = 0, lty = 2, lwd = 2) # 添加垂直0阶段
  
  # 添加四个象限的散点
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#E64B35FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#F39B7FFF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#00A087FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#4DBBD5FF")
  
  # 绘制图例
  legend("topleft",
         legend = c("Subgroup",paste0(GeneName1,"- & ",GeneName2,"-"),paste0(GeneName1,"- & ",GeneName2,"+"),paste0(GeneName1,"+ & ",GeneName2,"-"),paste0(GeneName1,"+ & ",GeneName2,"+")),
         col = c(NA,"#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
         pch = 19,
         border = NA,
         bty = "n")
  dev.off()
  survival_time_data = getPanGeneTime(GeneName2)
  survival_time_data = survival_time_data[,c("OS.time","OS","DSS.time","DSS","DFI.time","DFI","PFI.time","PFI")]
  dat = getmerge_row(survival_time_data,GeneName1_GeneName2_merge_data_scale)
  return(dat)
}

get2varibles_mean_external <- function(GeneName1_GeneName2_merge_data,GeneName1,GeneName2,CancerName,DatasetName){
  GeneName1_GeneName2_merge_data_scale = GeneName1_GeneName2_merge_data
  GeneName1_GeneName2_merge_data_scale = as.data.frame(GeneName1_GeneName2_merge_data_scale)
  GeneName1_GeneName2_merge_data_scale$Subtype <- 
    ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] >= 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] >= 0,paste0(GeneName1,"+ & ",GeneName2,"+"),
           ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] <= 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] <= 0,paste0(GeneName1,"- & ",GeneName2,"-"),
                  ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] > 0 & GeneName1_GeneName2_merge_data_scale[,GeneName2] < 0,paste0(GeneName1,"+ & ",GeneName2,"-"),
                         paste0(GeneName1,"- & ",GeneName2,"+"))))
  write.csv(GeneName1_GeneName2_merge_data_scale,paste0(GeneName1,"_",GeneName2,"_",DatasetName,"_merge_data_scale.csv"))
  
  #-- -+ +- ++
  palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF")
  pdf(file = paste0(GeneName1,"_",GeneName2,"_",CancerName,"_Mean_classification.pdf"),width = 5, height = 5)
  par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,4.1,2.1,3.1)+.1, las=1, tcl=-.25)
  plot(NULL, NULL, # 绘制空白背景
       xlim =c(-4,4),
       ylim =c(-4,4),
       xlab = paste0(GeneName1," (z-score)"),
       ylab = paste0(gsub("_"," ",GeneName2)," (z-score)"))
  grid(col = "grey85", lty = 2, lwd = 1.5) # 添加网格线
  abline(v = 0, lty = 2, lwd = 2) # 添加水平0截断
  abline(h = 0, lty = 2, lwd = 2) # 添加垂直0阶段
  
  # 添加四个象限的散点
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#E64B35FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#F39B7FFF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#00A087FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#4DBBD5FF")
  
  # 绘制图例
  legend("topleft",
         legend = c("Subgroup",paste0(GeneName1,"- & ",GeneName2,"-"),paste0(GeneName1,"- & ",GeneName2,"+"),paste0(GeneName1,"+ & ",GeneName2,"-"),paste0(GeneName1,"+ & ",GeneName2,"+")),
         col = c(NA,"#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
         pch = 19,
         border = NA,
         bty = "n")
  dev.off()
  
  return(GeneName1_GeneName2_merge_data_scale)
}

get2varibles_median <- function(GeneName1,GeneName2,CancerName){
  GeneName1_data = getPanGeneTime(GeneName1)
  GeneName1_data = GeneName1_data[(GeneName1_data[,"CancerType"] == CancerName),]
  GeneName1_data = GeneName1_data[,GeneName1,drop=F]
  GeneName2_data = getPanGeneTime(GeneName2)
  GeneName2_data = GeneName2_data[(GeneName2_data[,"CancerType"] == CancerName),]
  GeneName2_data = GeneName2_data[,GeneName2,drop=F]
  GeneName1_GeneName2_merge_data = getmerge_row(GeneName1_data,GeneName2_data)
  GeneName1_GeneName2_merge_data_scale = scale(GeneName1_GeneName2_merge_data, center = TRUE, scale = TRUE)
  GeneName1_GeneName2_merge_data_scale =as.data.frame(GeneName1_GeneName2_merge_data_scale)
  GeneName1_GeneName2_merge_data_scale$Subtype <- 
    ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] >= median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] >= median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"+ & ",GeneName2,"+"),
           ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] <= median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] <= median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"- & ",GeneName2,"-"),
                  ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] > median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] < median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"+ & ",GeneName2,"-"),
                         paste0(GeneName1,"- & ",GeneName2,"+"))))
  write.csv(GeneName1_GeneName2_merge_data_scale,paste0(GeneName1,"_",GeneName2,"_merge_data_scale.csv"))
  
  #-- -+ +- ++
  palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF")
  pdf(file = paste0(GeneName1,"_",GeneName2,"_",CancerName,"_Median_classification.pdf"),width = 5, height = 5)
  par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,4.1,2.1,3.1)+.1, las=1, tcl=-.25)
  plot(NULL, NULL, # 绘制空白背景
       xlim =c(-max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName1])),max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName1]))),
       ylim =c(-max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName2])),max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName2]))),
       xlab = paste0(GeneName1," (z-score)"),
       ylab = paste0(gsub("_"," ",GeneName2)," (z-score)"))
  grid(col = "grey85", lty = 2, lwd = 1.5) # 添加网格线
  abline(v = median(GeneName1_GeneName2_merge_data_scale[,GeneName1]), lty = 2, lwd = 2) # 添加水平0截断
  abline(h = median(GeneName1_GeneName2_merge_data_scale[,GeneName2]), lty = 2, lwd = 2) # 添加垂直0阶段
  
  # 添加四个象限的散点
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#E64B35FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#F39B7FFF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#00A087FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#4DBBD5FF")
  
  # 绘制图例
  legend("topleft",
         legend = c("Subgroup",paste0(GeneName1,"- & ",GeneName2,"-"),paste0(GeneName1,"- & ",GeneName2,"+"),paste0(GeneName1,"+ & ",GeneName2,"-"),paste0(GeneName1,"+ & ",GeneName2,"+")),
         col = c(NA,"#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
         pch = 19,
         border = NA,
         bty = "n")
  dev.off()
  survival_time_data = getPanGeneTime(GeneName2)
  survival_time_data = survival_time_data[,c("OS.time","OS","DSS.time","DSS","DFI.time","DFI","PFI.time","PFI")]
  dat = getmerge_row(survival_time_data,GeneName1_GeneName2_merge_data_scale)
  return(dat)
}

get2varibles_median_external <- function(GeneName1_GeneName2_merge_data,GeneName1,GeneName2,CancerName,DatasetName){
  GeneName1_GeneName2_merge_data_scale = GeneName1_GeneName2_merge_data
  GeneName1_GeneName2_merge_data_scale = as.data.frame(GeneName1_GeneName2_merge_data_scale)
  GeneName1_GeneName2_merge_data_scale$Subtype <- 
    ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] >= median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] >= median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"+ & ",GeneName2,"+"),
           ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] <= median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] <= median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"- & ",GeneName2,"-"),
                  ifelse(GeneName1_GeneName2_merge_data_scale[,GeneName1] > median(GeneName1_GeneName2_merge_data_scale[,GeneName1]) & GeneName1_GeneName2_merge_data_scale[,GeneName2] < median(GeneName1_GeneName2_merge_data_scale[,GeneName2]),paste0(GeneName1,"+ & ",GeneName2,"-"),
                         paste0(GeneName1,"- & ",GeneName2,"+"))))
  write.csv(GeneName1_GeneName2_merge_data_scale,paste0(GeneName1,"_",GeneName2,"_",DatasetName,"_merge_data_scale.csv"))
  
  #-- -+ +- ++
  palette = c("#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF")
  pdf(file = paste0(GeneName1,"_",GeneName2,"_",CancerName,"_Median_classification.pdf"),width = 5, height = 5)
  par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,4.1,2.1,3.1)+.1, las=1, tcl=-.25)
  plot(NULL, NULL, # 绘制空白背景
       xlim =c(-max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName1])),max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName1]))),
       ylim =c(-max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName2])),max(abs(GeneName1_GeneName2_merge_data_scale[,GeneName2]))),
       xlab = paste0(GeneName1," (z-score)"),
       ylab = paste0(gsub("_"," ",GeneName2)," (z-score)"))
  grid(col = "grey85", lty = 2, lwd = 1.5) # 添加网格线
  abline(v = median(GeneName1_GeneName2_merge_data_scale[,GeneName1]), lty = 2, lwd = 2) # 添加水平0截断
  abline(h = median(GeneName1_GeneName2_merge_data_scale[,GeneName2]), lty = 2, lwd = 2) # 添加垂直0阶段
  
  # 添加四个象限的散点
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#E64B35FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"+ & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#F39B7FFF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"+")),GeneName2],
         pch = 19,
         col = "#00A087FF")
  
  points(GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName1],
         GeneName1_GeneName2_merge_data_scale[which(GeneName1_GeneName2_merge_data_scale$Subtype == paste0(GeneName1,"- & ",GeneName2,"-")),GeneName2],
         pch = 19,
         col = "#4DBBD5FF")
  
  # 绘制图例
  legend("topleft",
         legend = c("Subgroup",paste0(GeneName1,"- & ",GeneName2,"-"),paste0(GeneName1,"- & ",GeneName2,"+"),paste0(GeneName1,"+ & ",GeneName2,"-"),paste0(GeneName1,"+ & ",GeneName2,"+")),
         col = c(NA,"#4DBBD5FF","#00A087FF","#F39B7FFF","#E64B35FF"),
         pch = 19,
         border = NA,
         bty = "n")
  dev.off()
  return(GeneName1_GeneName2_merge_data_scale)
}

getChisqProportion_OS <- function(data_OS,DatasetName){
  data_OS$OS = as.vector(ifelse(data_OS$OS == 0,"Alive","Dead"))
  data_OS = data_OS[,c("Quantile","OS")]
  colnames(data_OS) = c("Group","Status")
  table_data_OS = table(data_OS)
  df_table_data_OS = as.data.frame(table_data_OS)
  chisq.test_result = chisq.test(table_data_OS)
  pvalue = chisq.test_result$p.value
  if(pvalue > 0.05){
  }else{
    if(pvalue < 0.001){
      pvalue = "p < 0.001"
    }else{
      pvalue = paste0("p = ",round(pvalue,3))
    }
    library(plyr)
    df = ddply(df_table_data_OS, .(Group), transform, percent = Freq/sum(Freq) * 100)
    df = ddply(df, .(Group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label = paste0(sprintf("%.0f", df$percent), "%")
    df$Status = factor(df$Status, levels = c("Alive","Dead"))
    library(ggplot2)
    p = ggplot(df, 
               aes(x = factor(Group), 
                   y = percent,  
                   fill = Status)
    ) +
      geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Dead" = "#E64B35FF", 
                                   "Alive" = "#4DBBD5FF")
      ) +
      ggtitle(paste0("Chi-Square test ",pvalue,"(",gsub("_","-",DatasetName),")")) +
      xlab("") + 
      ylab("Percent") +  
      guides(fill = guide_legend(title="Status")
      ) +
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      theme_bw() +
      theme(axis.text.x = ggtext::element_markdown(angle =90, 
                                                   size = 12, 
                                                   hjust = 1, 
                                                   vjust = 0.5),
            plot.margin = unit(c(1,1,1,1), "lines")
      )
    pdf(file=paste0(DatasetName,"_Survival_Chi-Square_test_OS.pdf"), width=5, height=5)
    print(p)
    dev.off()
  }
}

getChisqProportion_PFS <- function(data_PFS,DatasetName){
  data_PFS$PFS = as.vector(ifelse(data_PFS$PFS == 0,"Progression-free","Progression"))
  data_PFS = data_PFS[,c("Quantile","PFS")]
  colnames(data_PFS) = c("Group","Status")
  table_data_PFS = table(data_PFS)
  df_table_data_PFS = as.data.frame(table_data_PFS)
  chisq.test_result = chisq.test(table_data_PFS)
  pvalue = chisq.test_result$p.value
  if(pvalue > 0.05){
  }else{
    if(pvalue < 0.001){
      pvalue = "p < 0.001"
    }else{
      pvalue = paste0("p = ",round(pvalue,3))
    }
    library(plyr)
    df = ddply(df_table_data_PFS, .(Group), transform, percent = Freq/sum(Freq) * 100)
    df = ddply(df, .(Group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label = paste0(sprintf("%.0f", df$percent), "%")
    df$Status = factor(df$Status, levels = c("Progression-free","Progression"))
    library(ggplot2)
    p = ggplot(df, 
               aes(x = factor(Group), 
                   y = percent,  
                   fill = Status)
    ) +
      geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Progression" = "#E64B35FF", 
                                   "Progression-free" = "#4DBBD5FF")
      ) +
      ggtitle(paste0("Chi-Square test ",pvalue,"(",gsub("_","-",DatasetName),")")) +
      xlab("") + 
      ylab("Percent") +  
      guides(fill = guide_legend(title="Status")
      ) +
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      theme_bw() +
      theme(axis.text.x = ggtext::element_markdown(angle =90, 
                                                   size = 12, 
                                                   hjust = 1, 
                                                   vjust = 0.5),
            plot.margin = unit(c(1,1,1,1), "lines")
      )
    pdf(file=paste0(DatasetName,"_Survival_Chi-Square_test_PFS.pdf"), width=5, height=5)
    print(p)
    dev.off()
  }
}

getChisqProportion_RFS <- function(data_RFS,DatasetName){
  data_RFS$RFS = as.vector(ifelse(data_RFS$RFS == 0,"Relapse-free","Relapse"))
  data_RFS = data_RFS[,c("Quantile","RFS")]
  colnames(data_RFS) = c("Group","Status")
  table_data_RFS = table(data_RFS)
  df_table_data_RFS = as.data.frame(table_data_RFS)
  chisq.test_result = chisq.test(table_data_RFS)
  pvalue = chisq.test_result$p.value
  if(pvalue > 0.05){
  }else{
    if(pvalue < 0.001){
      pvalue = "p < 0.001"
    }else{
      pvalue = paste0("p = ",round(pvalue,3))
    }
    library(plyr)
    df = ddply(df_table_data_RFS, .(Group), transform, percent = Freq/sum(Freq) * 100)
    df = ddply(df, .(Group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label = paste0(sprintf("%.0f", df$percent), "%")
    df$Status = factor(df$Status, levels = c("Relapse-free","Relapse"))
    library(ggplot2)
    p = ggplot(df, 
               aes(x = factor(Group), 
                   y = percent,  
                   fill = Status)
    ) +
      geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Relapse" = "#E64B35FF", 
                                   "Relapse-free" = "#4DBBD5FF")
      ) +
      ggtitle(paste0("Chi-Square test ",pvalue,"(",gsub("_","-",DatasetName),")")) +
      xlab("") + 
      ylab("Percent") +  
      guides(fill = guide_legend(title="Status")
      ) +
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      theme_bw() +
      theme(axis.text.x = ggtext::element_markdown(angle =90, 
                                                   size = 12, 
                                                   hjust = 1, 
                                                   vjust = 0.5),
            plot.margin = unit(c(1,1,1,1), "lines")
      )
    pdf(file=paste0(DatasetName,"_Survival_Chi-Square_test_RFS.pdf"), width=5, height=5)
    print(p)
    dev.off()
  }
}

getChisqProportion_DFS <- function(data_DFS,DatasetName){
  data_DFS$DFS = as.vector(ifelse(data_DFS$DFS == 0,"Disease-free","Disease"))
  data_DFS = data_DFS[,c("Quantile","DFS")]
  colnames(data_DFS) = c("Group","Status")
  table_data_DFS = table(data_DFS)
  df_table_data_DFS = as.data.frame(table_data_DFS)
  chisq.test_result = chisq.test(table_data_DFS)
  pvalue = chisq.test_result$p.value
  if(pvalue > 0.05){
  }else{
    if(pvalue < 0.001){
      pvalue = "p < 0.001"
    }else{
      pvalue = paste0("p = ",round(pvalue,3))
    }
    library(plyr)
    df = ddply(df_table_data_DFS, .(Group), transform, percent = Freq/sum(Freq) * 100)
    df = ddply(df, .(Group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label = paste0(sprintf("%.0f", df$percent), "%")
    df$Status = factor(df$Status, levels = c("Disease-free","Disease"))
    library(ggplot2)
    p = ggplot(df, 
               aes(x = factor(Group), 
                   y = percent,  
                   fill = Status)
    ) +
      geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Disease" = "#E64B35FF", 
                                   "Disease-free" = "#4DBBD5FF")
      ) +
      ggtitle(paste0("Chi-Square test ",pvalue,"(",gsub("_","-",DatasetName),")")) +
      xlab("") + 
      ylab("Percent") +  
      guides(fill = guide_legend(title="Status")
      ) +
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      theme_bw() +
      theme(axis.text.x = ggtext::element_markdown(angle =90, 
                                                   size = 12, 
                                                   hjust = 1, 
                                                   vjust = 0.5),
            plot.margin = unit(c(1,1,1,1), "lines")
      )
    pdf(file=paste0(DatasetName,"_Survival_Chi-Square_test_DFS.pdf"), width=5, height=5)
    print(p)
    dev.off()
  }
}

getChisqProportion_DSS <- function(data_DSS,DatasetName){
  data_DSS$DSS = as.vector(ifelse(data_DSS$DSS == 0,"Alive","Dead"))
  data_DSS = data_DSS[,c("Quantile","DSS")]
  colnames(data_DSS) = c("Group","Status")
  table_data_DSS = table(data_DSS)
  df_table_data_DSS = as.data.frame(table_data_DSS)
  chisq.test_result = chisq.test(table_data_DSS)
  pvalue = chisq.test_result$p.value
  if(pvalue > 0.05){
  }else{
    if(pvalue < 0.001){
      pvalue = "p < 0.001"
    }else{
      pvalue = paste0("p = ",round(pvalue,3))
    }
    library(plyr)
    df = ddply(df_table_data_DSS, .(Group), transform, percent = Freq/sum(Freq) * 100)
    df = ddply(df, .(Group), transform, pos = (cumsum(Freq) - 0.5 * Freq))
    df$label = paste0(sprintf("%.0f", df$percent), "%")
    df$Status = factor(df$Status, levels = c("Alive","Dead"))
    library(ggplot2)
    p = ggplot(df, 
               aes(x = factor(Group), 
                   y = percent,  
                   fill = Status)
    ) +
      geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Dead" = "#E64B35FF", 
                                   "Alive" = "#4DBBD5FF")
      ) +
      ggtitle(paste0("Chi-Square test ",pvalue,"(",gsub("_","-",DatasetName),")")) +
      xlab("") + 
      ylab("Percent") +  
      guides(fill = guide_legend(title="Status")
      ) +
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      theme_bw() +
      theme(axis.text.x = ggtext::element_markdown(angle =90, 
                                                   size = 12, 
                                                   hjust = 1, 
                                                   vjust = 0.5),
            plot.margin = unit(c(1,1,1,1), "lines")
      )
    pdf(file=paste0(DatasetName,"_Survival_Chi-Square_test_DSS.pdf"), width=5, height=5)
    print(p)
    dev.off()
  }
}

getlist2dataframe <- function(original_list){
  gene_data_list = lapply(names(original_list), function(pathway) {  
    data.frame(  
      term = rep(pathway, length(original_list[[pathway]])),  
      gene = original_list[[pathway]],  
      stringsAsFactors = FALSE  
    )  
  })  
  transformed_df = do.call(rbind, gene_data_list) 
  return(transformed_df)
}

getgmt2list <- function(gmtFilePath){
  geneSets = GSEABase::getGmt(gmtFilePath)  
  return(geneSets)  
}  

getkegggseaplot = function(gsea_Result = gsea_Result,DatasetName){
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(magrittr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(ggplot2))
  gsea_Result_df = as.data.frame(gsea_Result)
  data_ann = read.table("/mnt/Merge/KEGG.ann.txt", sep="\t", header=TRUE, check.names=FALSE)
  KEGG_results = merge(gsea_Result_df, data_ann, by.x = "ID", by.y = "pathway")
  KEGG_results_plot = KEGG_results[,c("ID.y","Description.y","NES","pvalue","p.adjust")]
  colnames(KEGG_results_plot) = c("ID","Description","NES","pvalue","p.adjust")
  KEGG_results_plot = KEGG_results_plot[which(as.numeric(KEGG_results_plot$pvalue) < 0.05 & as.numeric(KEGG_results_plot$p.adjust) < 0.25),]
  KEGG_results_plot$ID = gsub("hsa","ko",KEGG_results_plot$ID)
  KEGG_results_plot$Direction = ifelse(KEGG_results_plot$NES > 0,"Activated", "Inhibited")
  KEGG_results_plot$NES = round(KEGG_results_plot$NES, 2)
  KEGG_results_plot$ID = ifelse(KEGG_results_plot$Description == "Aminoacyl-tRNA biosynthesis","ko00970",KEGG_results_plot$ID)
  KEGG_results_plot$ID = ifelse(KEGG_results_plot$Description == "Spliceosome","ko03040",KEGG_results_plot$ID)
  write.csv(KEGG_results_plot,paste0(DatasetName,"_KEGG_results_plot.csv"),row.names = F)
  KEGG_ann = read.csv(file = "/mnt/Merge/kegg_hierarchy.csv")
  levels1_2 = as.vector(data.frame(table(KEGG_ann$level1))[,1])
  levels2_2 = as.vector(data.frame(table(KEGG_ann$level2))[,1])
  ggData = KEGG_results_plot %>%
    mutate(level3 = KEGG_ann$level3[match(ID, KEGG_ann$id)],
           level2 = KEGG_ann$level2[match(ID, KEGG_ann$id)],
           level1 = KEGG_ann$level1[match(ID, KEGG_ann$id)]) %>%
    mutate(level1 = factor(level1,
                           c(levels1_2)),
           level2 = factor(level2,
                           c(levels2_2))
    ) %>%
    arrange(level1, level2, Direction, level3) %>%
    mutate(ID = factor(ID, rev(unique(ID))),
           level3_x = factor(level3, rev(unique(as.character(level3)))))
  head(ggData)
  ggData_l1 = data.frame(nrow(ggData) - (cumsum(table(ggData$level1)) - table(ggData$level1)/2))
  ggData_l1$start = ggData_l1$Freq - table(ggData$level1)/2
  ggData_l1$end = ggData_l1$Freq + table(ggData$level1)/2
  levels1_2 = data.frame(table(ggData$level1))
  levels1_2 = levels1_2[-which(levels1_2[,2] %in% c("0")),]
  levels1_2 = as.vector(levels1_2[,1])
  ggData_l1 = ggData_l1[ggData_l1$Var1 %in% levels1_2,]
  ggData$level2 = factor(ggData$level2,levels = unique(ggData$level2))
  ggData_l2 = data.frame(nrow(ggData) - (cumsum(table(ggData$level2)) - table(ggData$level2)/2))
  ggData_l2$start = ggData_l2$Freq - table(ggData$level2)/2
  ggData_l2$end = ggData_l2$Freq + table(ggData$level2)/2
  # bar 1
  ggplot(ggData) +
    geom_col(mapping = aes(ID, -log10(ggData$pvalue), fill = Direction),
             color = "white", 
             width = 0.75,
             show.legend = F) +
    geom_text(mapping = aes(ID, -log10(ggData$pvalue), label = round(-log10(ggData$pvalue),2)),
              #hjust = ifelse(ggData$NES > 0, 1, 0), 
              size = 2.5) +
    scale_y_continuous(limits = c(0,max(-log10(ggData$pvalue))*1.1),
                       expand = expansion()
    ) +
    scale_fill_manual(values = c(alpha("#E64B35FF",0.5), alpha("#4DBBD5FF",0.5))) +
    coord_flip() +
    theme_classic() +
    labs(x = NULL, y = NULL, title = "\n-log10(pvalue)") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p1
  # bar2
  ggplot(ggData) +
    geom_col(mapping = aes(ID, NES, fill = Direction),
             color = "white", width = 0.75, show.legend = F) +
    geom_text(mapping = aes(ID, NES, label = NES),
              hjust=ifelse(ggData$NES > 0, 1, 0), # 水平位置, 
              size = 2.5) +
    scale_y_continuous(limits = c(min(ggData$NES)*1.1,max(ggData$NES)*1.1),
                       expand = expansion()) +
    scale_fill_manual(values = c(alpha("#E64B35FF",0.5), alpha("#4DBBD5FF",0.5))) +
    coord_flip() +
    theme_classic() +
    labs(x = NULL, y = NULL, title = "\nNES") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p2
  # text ko
  ggplot(ggData) +
    geom_text(mapping = aes(ID, 0, label = ID, color = Direction),
              size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")) +
    scale_y_continuous(expand = expansion(), limits = c(0,1)) +
    coord_flip() +
    theme_void() +
    labs(x = NULL, y = NULL, title = "Pathway") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p0
  
  # text level3
  ggplot(ggData) +
    geom_text(mapping = aes(ID, 0, label = level3, color = Direction),
              size = 3, show.legend = F, hjust = 0) +
    scale_color_manual(values = c( "#E64B35FF", "#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")) +
    scale_y_continuous(expand = expansion(), limits = c(0,1)) +
    coord_flip() +
    theme_void() +
    labs(x = NULL, y = NULL, title = "Level 3 of KEGG\nfunctional Category") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p3
  
  # text level2
  ggplot(ggData_l2) +
    geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
    geom_text(mapping = aes(Freq, 0, label = Var1), #Freq
              size = 3, show.legend = F, hjust = 0) +
    scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
    scale_x_continuous(expand = expansion()) +
    coord_flip() +
    theme_void() +
    labs(x = NULL, y = NULL, title = "Level 2 of KEGG\nfunctional Category") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p4
  
  # text level1
  ggplot(ggData_l1) +
    geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
    geom_text(mapping = aes(Freq, 0, label = Var1),
              size = 3, show.legend = F, hjust = 0) +
    scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
    scale_x_continuous(expand = expansion()) +
    coord_flip() +
    theme_void() +
    labs(x = NULL, y = NULL, title = "Level 1 of KEGG\nfunctional Category") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8)) -> p5
  cowplot::plot_grid(p0,p1, p2,p3,p4,p5, align = "h", nrow = 1,
                     rel_widths = c(0.1,0.2,0.2,0.6, 0.5, 0.4))
  ggsave(paste0(DatasetName,"_KEGGhierarchy.pdf"), width = 12, height =4 + 0.125*nrow(KEGG_results_plot))
}

getspearman_cor1 = function(data,TME,nameit = paste0("_ST_cor.pdf")){
  #基因与免疫细胞相关性分析
  geneCor = data.frame()
  for(cell in colnames(TME)){
    for(gene in colnames(data)){
      x = as.numeric(TME[,cell])
      y = as.numeric(data[,gene])
      corT = cor.test(x, y, method="spearman")
      cor = corT$estimate
      pvalue = corT$p.value
      geneCor = rbind(geneCor, cbind(spec=gene, env=cell, r=cor, p=pvalue))
    }
  }
  geneCor$r = as.numeric(geneCor$r)
  geneCor$p = as.numeric(geneCor$p)
  geneCor$pd = ifelse(geneCor$p<0.05, ifelse(geneCor$r>0, "Postive", "Negative"), "Not")
  geneCor$r = abs(geneCor$r)
  geneCor = geneCor %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                                        labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6",">= 0.6")))
  
  #绘制图形
  qcorPlot = qcorrplot(correlate(TME, method="spearman"), type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(aes(colour = pd, size = rd), 
                data = geneCor, 
                curvature = nice_curvature()) +
    #设置图形的颜色和图例的名称
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
    scale_size_manual(values = c(0.5, 1, 2, 3)) +
    scale_colour_manual(values = c(Negative="#1B9E77", Not="#CCCCCC99", Postive="#D95F02","black")) +
    guides(size = guide_legend(title = "abs(Cor)",
                               override.aes = list(colour = "grey35"), 
                               order = 2),
           colour = guide_legend(title = "pvalue", 
                                 override.aes = list(size = 3), 
                                 order = 1),
           fill = guide_colorbar(title = "cell-cell cor", order = 3))
  
  #输出图形
  pdf(file = nameit, width=6, height=6)
  print(qcorPlot)
  dev.off()
}

