#/user/bin/Rscript
#integrate analysis for mulitple single cell sequencing samples, find cluster and markers
#any sample number
#self.R experifile sepcies

args=commandArgs(T)
species=args[2] # hsa or mmu

dir.create("outFigures")
dir.create("outData")

library(Seurat)
library(cowplot)
library(patchwork)
library(Matrix)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
#library(future)

#plan("multiprocess", workers = 8)

experi=read.table(args[1],header=TRUE,sep="\t",check.names=FALSE)#inputData/inputInfiles.txt
allsample=NULL

for(i in 1:length(experi$filePath)){
	d=read.table(experi$filePath[i],header=TRUE,sep="\t")
	d=d[!duplicated(d[,1]),]
	rownames(d)=d[,1]
	d=d[,-1]
	d=as.matrix(d)
	sparse <- Matrix(d, sparse = T )
	data=sparse
	# Initialize the Seurat object with the raw (non-normalized data).
	data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
	# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
	if(species == "hsa"){
		print("human")
		data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") #human
	}else if(species == "mmu"){
		print("mouse")
		data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-") #mouse
	}else if(species == "rno"){
		print("rat")
		data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^Mt-") #rat
	}
	
	for(j in 2:length(names(experi))){
		data[[names(experi)[j]]]=experi[i,j]
	}
	# Visualize QC metrics as a violin plot
	QC1=VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(paste0("outFigures/",experi$sample[i],".QC1.pdf"),QC1,width=16,height=9)
	# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
	# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
	plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	QC2=CombinePlots(plots = list(plot1, plot2))
	ggsave(paste0("outFigures/",experi$sample[i],".QC2.pdf"),QC2,width=16,height=9)
	#filter cells
	data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
	allsample=append(allsample,data)
}

if(length(allsample) == 1){
		alls=allsample[[1]]
}else{
		alls1=allsample
		alls1[[1]]=NULL
		alls=merge(x=allsample[[1]],y=alls1)
}

ifnb.list <- SplitObject(alls, split.by = "sample")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)

anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

combined <- IntegrateData(anchorset = anchors)


DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# UMAP and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

p1 <- DimPlot(combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file="outFigures/seurat_integrated_umap.pdf",width=16)
plot_grid(p1, p2)
dev.off()

pdf(file="outFigures/seurat_integrated_umap_splitgroup.pdf",width=16)
print(DimPlot(combined, reduction = "umap", split.by = "group"))
dev.off()


DefaultAssay(combined) <- "RNA"


marker_table=NULL

all_diff_gene=NULL

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR", latent.vars = "sample")
markers=markers[markers$p_val_adj <= 0.05,]
all_diff_gene=markers
marker_table=markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) # top 10

for(i in sort(unique(as.vector(combined$seurat_clusters)))){
	pdf(file=paste0(paste0("outFigures/cluster",i),"_top10_markers.pdf"),width=16,height=10.5)
	print(FeaturePlot(combined, features = marker_table[marker_table$cluster == i,]$gene, min.cutoff = "q10"))
	dev.off()
	p=VlnPlot(combined, features = marker_table[marker_table$cluster == i,]$gene)
	ggsave(paste0(paste0("outFigures/cluster",i),"_top10_markers_vlnplot.pdf"),p,width=20,height=12)

}

####heatmap####

###2fc###
combined <- ScaleData(combined, features = unique(markers[markers$avg_log2FC >= 1,]$gene)) #2fc
#pdf(file="outFigures/cluster_2fc_markers_heatmap.pdf",width=20,height=16)
#print(DoHeatmap(combined, features = markers[markers$avg_log2FC >= 1,]$gene,size = 3,draw.lines = FALSE)+theme(axis.text.y = element_blank()))
#dev.off()

####average heatmap####
cluster.averages <- AverageExpression(combined, return.seurat = TRUE)
p=DoHeatmap(object = cluster.averages,features=markers[markers$avg_log2FC >= 1,]$gene,size = 3,draw.lines = FALSE)+theme(axis.text.y = element_blank())

ggsave("outFigures/cluster_2fc_markers_heatmap_averageexp.pdf",p,width=8,height=12)

####top10
combined <- ScaleData(combined, features = unique(marker_table$gene)) #top 10
pdf(file="outFigures/cluster_top10_markers_heatmap.pdf",width=20,height=20)
print(DoHeatmap(combined, features = marker_table$gene,size = 3,draw.lines = FALSE))
dev.off()

####average heatmap####
cluster.averages <- AverageExpression(combined, return.seurat = TRUE)
p=DoHeatmap(object = cluster.averages,features=marker_table$gene,size = 3,draw.lines = FALSE)

ggsave("outFigures/cluster_top10_markers_heatmap_averageexp.pdf",p,width=10,height=20)

###cell number bar plot####
cell <- table(combined$sample) %>% as.data.frame()
colnames(cell) <- c("sample", "number")
cell$sample=as.factor(cell$sample)

p=ggplot(cell, aes(sample, number, fill=sample))+
  geom_bar(stat="identity")+
  #coord_flip()+
  #scale_fill_brewer(palette="Paired")+
  geom_text(aes(label = cell$number),vjust = -0.2, size = 6)+
  ggtitle("")+
  theme_base()+
  ylab("number of cells")+
  xlab("")+
  theme(axis.text = element_text(color = "black", face = "bold", size = 12),
        axis.title = element_text(color = "black", face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

ggsave("outFigures/cell_number_stat.pdf",p,height= 6,width = 6)

###ratio for each sample
xyz=table(Idents(combined), combined$sample)
write.table(xyz,"outData/cluster_countcells_stat_forsample.xls",col.names = NA,sep="\t",quote = FALSE)

cellration <- prop.table(xyz, margin = 2) %>% as.data.frame()

colnames(cellration)=c("cluster","sample","Freq")
cellration$cluster=as.factor(cellration$cluster)

p=ggplot(cellration, aes(sample, Freq, fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  #scale_fill_brewer(palette="Paired")+
  ggtitle("")+
  theme_base()+
  ylab("Fraction of cells")+
  xlab("")+
  theme(axis.text = element_text(color = "black", face = "bold", size = 12),
        axis.title = element_text(color = "black", face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

ggsave("outFigures/cluster_ratio_stat_forsample.pdf",p,height= 6,width = 6)

###ratio for each group
xyz=table(Idents(combined), combined$group)
write.table(xyz,"outData/cluster_countcells_stat_forgroup.xls",col.names = NA,sep="\t",quote = FALSE)

cellration <- prop.table(xyz, margin = 2) %>% as.data.frame()

colnames(cellration)=c("cluster","group","Freq")
cellration$cluster=as.factor(cellration$cluster)

p=ggplot(cellration, aes(group, Freq, fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  #scale_fill_brewer(palette="Paired")+
  ggtitle("")+
  theme_base()+
  ylab("Fraction of cells")+
  xlab("")+
  theme(axis.text = element_text(color = "black", face = "bold", size = 12),
        axis.title = element_text(color = "black", face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

ggsave("outFigures/cluster_ratio_stat_forgroup.pdf",p,height= 6,width = 6)

##find conserve markers across conditions
#for(i in sort(unique(as.vector(combined$seurat_clusters)))){
#	if(length(which(combined@meta.data$group=="treat" & combined@meta.data$seurat_clusters==i)) > 2 && length(which(combined@meta.data$group=="ctrl" & combined@meta.data$seurat_clusters==i)) > 2){
#		markers <- FindConservedMarkers(combined, ident.1 = i, grouping.var = "group", verbose = FALSE,only.pos=TRUE)
#		#markers=markers[markers$minimump_p_val < 0.05 & markers$STIM_avg_logFC > 1 & markers$STIM_pct.1 >= 0.25,]
#		markers=markers[markers$minimump_p_val <= 0.05,]
#		diff_gene=data.frame(markers,cluster=i)
#		all_diff_gene=rbind(all_diff_gene,diff_gene)
#		markers=markers[order(markers$treat_avg_log2FC,decreasing=T),][1:10,] #top 10 marker
#		markers=data.frame(markers,cluster=i)
#		print(rownames(markers))
#		markers.to.plot=append(markers.to.plot,as.character(rownames(markers))[1:3]) #top 3 markers
#		pdf(file=paste0(paste0("outFigures/cluster",i),"_top10_markers.pdf"),width=12,height=16)
#		print(FeaturePlot(combined, features = as.character(rownames(markers)), min.cutoff = "q10"))
#		dev.off()
#		marker_table=rbind(marker_table,markers)
#	}
#}
##find conserve markers across conditions

write.table(marker_table,"outData/top10marker.txt",quote = FALSE,sep="\t",col.names=NA)
write.table(all_diff_gene,"outData/all.diffgene_onlypos_padj005.txt",quote = FALSE,sep="\t",col.names=NA)


saveRDS(combined,"outData/seuratObject.rds")


#if(args[2] != "NA"){
#if(args[3] == "NA"){
##annotate clusters using singleR
#refdata=readRDS(args[2])
#test=GetAssayData(combined,slot = 'data')
#pred1 <- SingleR(test=test, ref=refdata, labels = refdata$label.main, 
#    clusters=cluster$seurat_clusters)
#pred2 <- SingleR(test=test, ref=refdata, labels = refdata$label.fine, 
#    clusters=cluster$seurat_clusters)
#
#write.table(pred1,"outData/main_celltype_from_builtIn_singleR.xls",col.names=NA,sep="\t",quote=FALSE)
#write.table(pred2,"outData/cellsubtype_from_builtIn_singleR.xls",col.names=NA,sep="\t",quote=FALSE)
#}else{
#
##custome refdata set
#	refdata=readRDS(args[2])
#	label=colData(refdata)[,as.numeric(args[3])]
#	pred1 <- SingleR(test=test, ref=refdata, labels = label, 
#	    clusters=cluster$seurat_clusters)
#
#	write.table(pred1,"outData/main_celltype_from_customRef_singleR.xls",col.names=NA,sep="\t",quote=FALSE)
#}
#}

