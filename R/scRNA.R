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
    return(scRNA)
  }
}

