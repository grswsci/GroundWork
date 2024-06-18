getGPL <- function(GPLID){
  data_file <- system.file("data", paste0(GPLID,".RDS"), package = "GroundWork")
  GPL <- readRDS(data_file)
  GPL <- GPL[GPL$`Gene Symbol`!="",]
  return(GPL)
}

