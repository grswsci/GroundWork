scRead10X <- function(dir,version="v5"){
  library(Seurat)
  samples_name=as.vector(list.files(dir))
  samples_name
  sampes_name_path = paste0(dir,"/",samples_name)
  sampes_name_path
  scRNAlist <- list()
  for(i in sampes_name_path){
    project <- gsub(paste0(dir,"/"),"",i)
    counts <- Read10X(data.dir = i)
    scRNAlist[[i]] <- CreateSeuratObject(counts, 
                                         project = project,
                                         min.cells = 3, 
                                         min.features = 200)
    scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project)
  }
  names(scRNAlist) <- samples_name
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
  scRNA <- JoinLayers(scRNA)
  if(version=="v5"){
    return(scRNA)
  }else{
    scRNA[["RNA"]] <- as(scRNA[["RNA"]], "Assay")
    return(scRNA)
  }
}

scQC <- function(scRNA,
                 pctMT=10,
                 minGene=500,
                 maxGene=4000,
                 maxUMI=15000){
  scRNA[["percent.MT"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA <- subset(scRNA, subset = percent.MT < pctMT &
                                  nCount_RNA < maxUMI &
                                  nFeature_RNA > minGene &
                                  nFeature_RNA < maxGene)
}

scPipeline <- function(scRNA,
                       method = "LogNormalize",
                       resolution = 0.1,
                       integration = NULL){
  if(method=="LogNormalize"){
  scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 2000)
  scRNA <- ScaleData(scRNA)
  scRNA <- RunPCA(scRNA, npcs=20,pc.genes=VariableFeatures(object = scRNA), verbose=FALSE)
  pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2
  pcSelect <- min(co1, co2)
  print(pcSelect)
  if(is.null(integration)){
  scRNA <- FindNeighbors(scRNA, reduction="pca",dims = 1:pcSelect)
  scRNA <- FindClusters(object = scRNA,reduction="pca",resolution = resolution)
  scRNA <- RunUMAP(object = scRNA, reduction="pca",dims = 1:pcSelect,check_duplicates = FALSE)
  return(scRNA)
  }else{
    library(harmony)
    scRNA <- RunHarmony(object=scRNA, group.by.vars="orig.ident",reduction.use = "pca",assay.use="RNA",max.iter.harmony = 20)
    scRNA <- FindNeighbors(scRNA, reduction="harmony",dims = 1:pcSelect)
    scRNA <- FindClusters(object = scRNA,reduction="harmony",resolution = resolution)
    scRNA <- RunUMAP(object = scRNA, reduction="harmony",dims = 1:pcSelect,check_duplicates = FALSE)
    return(scRNA)
  }
  }else if(method=="SCT"){
    scRNA <- SCTransform(scRNA)
    scRNA <- RunPCA(scRNA,assay = "SCT", npcs=20, pc.genes=VariableFeatures(object = scRNA), verbose=FALSE)
    pct <- scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcSelect <- min(co1, co2)
    print(pcSelect)
    if(is.null(integration)){
      scRNA <- FindNeighbors(scRNA, reduction="pca",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="pca",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="pca",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }else{
      library(harmony)
      scRNA <- RunHarmony(object=scRNA, group.by.vars="orig.ident",reduction.use = "pca",assay.use="RNA",max.iter.harmony = 20)
      scRNA <- FindNeighbors(scRNA, reduction="harmony",dims = 1:pcSelect)
      scRNA <- FindClusters(object = scRNA,reduction="harmony",resolution = resolution)
      scRNA <- RunUMAP(object = scRNA, reduction="harmony",dims = 1:pcSelect,check_duplicates = FALSE)
      return(scRNA)
    }
  }
}

getUMAP <- function(scRNA,Labels="seurat_clusters") {
  library(tidyverse)
  library(randomcoloR)
  library(ggrepel)
  RandomColor <-distinctColorPalette(1500)
  palette = RandomColor 
  umap = scRNA@reductions$umap@cell.embeddings %>%
    as.data.frame() %>%
    cbind(Cluster = scRNA@meta.data[,Labels])
  head(umap)
  colnames(umap)=c("UMAP_1","UMAP_2","Cluster")
  allcolour=c(RandomColor)
  p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = Cluster)) +  
    geom_point(size = 0.1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)
  p
  p2 <- p + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          axis.title = element_blank(),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill="white"))
  p2
  p3 <- p2 + 
    theme(legend.title = element_blank(),
          legend.key=element_rect(fill='white'), 
          legend.text = element_text(size=20), 
          legend.key.size=unit(1,'cm') ) + 
    guides(color = guide_legend(override.aes = list(size=5)))
  p3
  p4 <- p3 + geom_segment(aes(x = min(umap$UMAP_1) ,
                              y = min(umap$UMAP_2),
                              xend = min(umap$UMAP_1) +3, 
                              yend = min(umap$UMAP_2) ),
                          colour = "black", 
                          size=1,
                          arrow = arrow(length = unit(0.3,"cm")))+
    geom_segment(aes(x = min(umap$UMAP_1)  ,
                     y = min(umap$UMAP_2),
                     xend = min(umap$UMAP_1) , 
                     yend = min(umap$UMAP_2) + 3),
                 colour = "black", 
                 size=1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    annotate("text",
             x = min(umap$UMAP_1) +
               1.5, 
             y = min(umap$UMAP_2) -1, 
             label = "UMAP_1",
             color="black",
             size = 3, 
             fontface="bold" ) + 
    annotate("text",
             x = min(umap$UMAP_1) -1, 
             y = min(umap$UMAP_2) + 1.5, 
             label = "UMAP_2",
             color="black",
             size = 3, fontface="bold" ,angle=90)
  p4
  cell_type_med <- umap %>% group_by(Cluster) %>%
    summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  p6 <- p4 +geom_label_repel(aes(label=Cluster),
                             fontface="bold",
                             data = cell_type_med,
                             point.padding=unit(0.5, "lines")) +
    theme(legend.position = "none")
  p6
  return(p6)
}