#!/user/bin/Rscript
#single cell cluster annotation,and differential analysis
#self.R cellann
#for no duplicate sample within each group

args=commandArgs(T)
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

combined=readRDS("seuratObject.rds")

cellanno=read.table(args[1],sep="\t",header=TRUE) #manual annotation results
cell_label=NULL
for(i in 1:length(levels(combined))){
	cells.use <- WhichCells(combined, idents = as.character(cellanno$clusters)[i])
	combined <- SetIdent(combined, cells = cells.use, value = as.character(cellanno$labels)[i])
}

p1 <- DimPlot(combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file="outFigures/seurat_integrated_cellann.pdf",width=16)
plot_grid(p1, p2)
dev.off()

pdf(file="outFigures/seurat_integrated_cellann_splitgroup.pdf",width=16)
print(DimPlot(combined, reduction = "umap", split.by = "group",label = TRUE, repel = TRUE))
dev.off()

#Identify differential expressed genes for sample cell type across conditions
combined$celltype.group <- paste(Idents(combined), combined$group, sep = "_")
combined$celltype <- Idents(combined)

Idents(combined) <- "celltype.group"

experi=read.table(args[2],header=TRUE,sep="\t",check.names=FALSE)
for(i in unique(cellanno$labels)){
	for(ca in experi$vs){
		if(ca != "no"){
			cas=strsplit(ca,"_vs_")     # 前面比后面
			if(length(which(combined@meta.data$group==cas[[1]][1] & combined@meta.data$celltype==i)) > 2 && length(which(combined@meta.data$group==cas[[1]][2] & combined@meta.data$celltype==i)) > 2){
				diff <- FindMarkers(combined, ident.1 = paste0(i,"_",cas[[1]][1]), ident.2 = paste0(i,"_",cas[[1]][2]), verbose = FALSE)
				write.table(diff,paste0("outData/",i,"_",ca,"_deanalysis.xls"),quote = FALSE,sep="\t",col.names=NA)
			}
		}
	}
}

####heatmap####

markers=read.table(args[3],header=TRUE,sep="\t",check.names=FALSE)  # only pos padj005

###2fc###
combined <- ScaleData(combined, features = unique(markers[markers$avg_log2FC >= 1,]$gene)) #2fc
#pdf(file="outFigures/cluster_2fc_markers_heatmap.pdf",width=20,height=16)
#print(DoHeatmap(combined, features = markers[markers$avg_log2FC >= 1,]$gene,size = 3,draw.lines = FALSE)+theme(axis.text.y = element_blank()))
#dev.off()
Idents(combined)=combined$celltype
levels(combined)=unique(cellanno$labels) #keep order
####average heatmap####
cluster.averages <- AverageExpression(combined, return.seurat = TRUE)
p=DoHeatmap(object = cluster.averages,features=markers[markers$avg_log2FC >= 1,]$gene,size = 3,draw.lines = FALSE)+theme(axis.text.y = element_blank())

ggsave("outFigures/celltype_2fc_markers_heatmap_averageexp.pdf",p,width=8,height=12)

####top10
marker_table=read.table(args[4],header=TRUE,sep="\t",check.names=FALSE)  # only top10

combined <- ScaleData(combined, features = unique(marker_table$gene)) #top 10
pdf(file="outFigures/celltype_top10_markers_heatmap.pdf",width=20,height=20)
print(DoHeatmap(combined, features = marker_table$gene,size = 3,draw.lines = FALSE))
dev.off()

####average heatmap####
cluster.averages <- AverageExpression(combined, return.seurat = TRUE)
p=DoHeatmap(object = cluster.averages,features=marker_table$gene,size = 3,draw.lines = FALSE)

ggsave("outFigures/celltype_top10_markers_heatmap_averageexp.pdf",p,width=10,height=20)


###ratio for each sample
xyz=table(Idents(combined), combined$sample)
write.table(xyz,"outData/celltype_countcells_stat_forsample.xls",col.names = NA,sep="\t",quote = FALSE)

cellration <- prop.table(xyz, margin = 2) %>% as.data.frame()

colnames(cellration)=c("celltype","sample","Freq")
cellration$celltype=as.factor(cellration$celltype)

p=ggplot(cellration, aes(sample, Freq, fill=celltype))+
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

ggsave("outFigures/celltype_ratio_stat_forsample.pdf",p,height= 6,width = 6)

###ratio for each group
xyz=table(Idents(combined), combined$group)
write.table(xyz,"outData/celltype_countcells_stat_forgroup.xls",col.names = NA,sep="\t",quote = FALSE)

cellration <- prop.table(xyz, margin = 2) %>% as.data.frame()

colnames(cellration)=c("celltype","group","Freq")
cellration$celltype=as.factor(cellration$celltype)

p=ggplot(cellration, aes(group, Freq, fill=celltype))+
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

ggsave("outFigures/celltype_ratio_stat_forgroup.pdf",p,height= 6,width = 6)

#export cell count matrix
counts_final=GetAssayData(combined,slot = 'counts')
tryCatch({
  write.table(counts_final,"outData/seurat_filter_counts_final.txt",col.names = NA,sep="\t",quote=FALSE)
},error=function(e){
  print(e)
  writeMM(obj = counts_final, file="seurat_filter_counts_final.mtx")
})
meta_data=combined@meta.data
write.table(meta_data,"outData/seurat_filter_counts_final_metadata.txt",col.names = NA,sep="\t",quote=FALSE)

saveRDS(combined,"outData/seuratObject.rds")

