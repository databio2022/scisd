library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
# library(future)
# plan("multiprocess", workers = 4)
options(future.globals.maxSize = 1e11)

# rds_list<-list.files('./',pattern='seuratObject.rds',recursive=T)
# rds_list<-rds_list[3:length(rds_list)]
# rds_list<-sub('/seuratObject.rds','',rds_list,fixed = TRUE)
# 
# obj_list<-list()
# count<-1
# 
# for(i in rds_list){
#   
#   dataset <- open_matrix_dir(dir = paste0('/data/leiwen/skin_disease_dataset/bpcells/',i))
#   metadata<-read.csv(paste0('/data/leiwen/skin_disease_dataset/bpcells/',i,'_metadata.csv'),row.names = 1)
#   metadata$sample<-paste(i,metadata$sample,sep='_')
#   skindisease <- CreateSeuratObject(counts = dataset,meta.data = metadata)
#   obj_list[[count]]<-skindisease
#   count<-count+1
#   
#   # gse<-sub('/seuratObject.rds','',i,fixed = TRUE)
#   # if (!file.exists(paste0('/data/leiwen/skin_disease_dataset/bpcells/',gse))){
#   #   obj<-readRDS(i)
#   #   obj <- RenameCells(obj, add.cell.id = gse)
#   #   
#   #   
#   #   
#   #   write_matrix_dir(mat = obj[["RNA"]]$counts, dir = paste0('/data/leiwen/skin_disease_dataset/bpcells/',gse))
#   #   metadata<-obj@meta.data
#   #   write.csv(metadata,paste0('/data/leiwen/skin_disease_dataset/bpcells/',gse,'_metadata.csv'),quote=FALSE)
#   # }
# }
# 
# obj<-merge(obj_list[[1]],obj_list[2:length(obj_list)],project='skin_disease')
# 
# obj[['RNA']]<-JoinLayers(obj[['RNA']])
# obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
# saveRDS(obj,'combined_stage1.rds')
obj<-readRDS('combined_stage1.rds')

obj[['RNA']]<-JoinLayers(obj[['RNA']])
split_sample<-stringr::str_split(obj$sample,'_')
split_sample_aa<-lapply(split_sample,FUN=function(x){x[1]})
obj$subject<-unlist(split_sample_aa)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$subject)
# bb<-names(table(obj$sample))[table(obj$sample)<100]
# obj<-obj[,!obj$sample%in%bb]
sketch_cell<-ceiling(300000/length(unique(obj$subject)))
print(sketch_cell)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = sketch_cell,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
saveRDS(obj,'combined_subject_stage2.rds')
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:15)
obj <- FindClusters(obj, resolution = 1, cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, dims = 1:15, return.model = T,reduction = "harmony")

obj<- ProjectIntegration(object=obj,sketched.assay='sketch',assay='RNA',reduction='harmony')

obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "harmony.full",
  sketched.assay = "sketch",
  sketched.reduction = "harmony",
  umap.model = "umap",
  dims = 1:15,
  refdata = list(cluster_full = "harmony_clusters")
)
obj <- RunUMAP(obj, reduction = "harmony.full", dims = 1:15, reduction.name = "umap.full",
                         reduction.key = "UMAP_full_")


saveRDS(obj,'combined_subject_stage3.rds')
obj[['sketch']]<-JoinLayers(obj[['sketch']])
Idents(obj)<-obj$harmony_clusters
markers <- FindAllMarkers(obj, only.pos = TRUE)
write.csv(markers,'subject_sketch_markers.csv')

