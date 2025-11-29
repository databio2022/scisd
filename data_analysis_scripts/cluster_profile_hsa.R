#!/user/bin/Rscript
#2fc 0.05
#mouse
#GO and KEGG analyses
#self.R inputfile sym_col

args=commandArgs(T)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
gene<-read.table(args[1],header = T,sep="\t") #header

symbol=unique(as.character(gene[,as.numeric(args[2])])) #GENESYMBOL
eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = unique(as.character(eg[,2]))

egoMF <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pAdjustMethod = "fdr",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)

egoBP <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)

egoCC <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "fdr",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)


MFdot=dotplot(egoMF, showCategory=10,color="pvalue") + scale_color_gradient(high="#21908CFF",low="orange")
BPdot=dotplot(egoBP, showCategory=10,color="pvalue") + scale_color_gradient(high="#21908CFF",low="orange")
CCdot=dotplot(egoCC, showCategory=10,color="pvalue") + scale_color_gradient(high="#21908CFF",low="orange")

MFbar=barplot(egoMF, showCategory=10,color="pvalue")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")
BPbar=barplot(egoBP, showCategory=10,color="pvalue")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")
CCbar=barplot(egoCC, showCategory=10,color="pvalue")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")

#p.adjust
MFdot1=dotplot(egoMF, showCategory=10,color="p.adjust") + scale_color_gradient(high="#21908CFF",low="orange")
BPdot1=dotplot(egoBP, showCategory=10,color="p.adjust") + scale_color_gradient(high="#21908CFF",low="orange")
CCdot1=dotplot(egoCC, showCategory=10,color="p.adjust") + scale_color_gradient(high="#21908CFF",low="orange")

MFbar1=barplot(egoMF, showCategory=10,color="p.adjust")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")
BPbar1=barplot(egoBP, showCategory=10,color="p.adjust")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")
CCbar1=barplot(egoCC, showCategory=10,color="p.adjust")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")

kk <- enrichKEGG(gene = id,
                 organism = 'hsa',
                 pvalueCutoff = 1,
                 pAdjustMethod = "fdr",
                 qvalueCutoff = 1,
                 use_internal_data = TRUE)
kk=setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kkdot=dotplot(kk,showCategory=10,color="pvalue") + scale_color_gradient(high="#21908CFF",low="orange")
kkbar=barplot(kk, showCategory=10,color="pvalue")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")

#p.adjust
kkdot1=dotplot(kk,showCategory=10,color="p.adjust") + scale_color_gradient(high="#21908CFF",low="orange")
kkbar1=barplot(kk,showCategory=10,color="p.adjust")+scale_fill_gradient(high="#21908CFF",low="#FDE725FF")

MF=as.data.frame(egoMF)
BP=as.data.frame(egoBP)
CC=as.data.frame(egoCC)
KEGG=as.data.frame(kk)

outdir=paste0(strsplit(args[1],"\\.")[[1]][1],"_GOKEGG")
dir.create(outdir)
setwd(paste0("./",outdir))
i="Results_"
write.table(egoMF,paste0(i,"geneMF.xls"),quote = FALSE,sep="\t",col.names=NA)
write.table(egoBP,paste0(i,"geneBP.xls"),quote = FALSE,sep="\t",col.names=NA)
write.table(egoCC,paste0(i,"geneCC.xls"),quote = FALSE,sep="\t",col.names=NA)
write.table(KEGG,paste0(i,"geneKEGG.xls"),quote = FALSE,sep="\t",col.names=NA)

ggsave(paste0(i,"geneMFdot.pdf"),MFdot,width=9,height=9)
ggsave(paste0(i,"geneBPdot.pdf"),BPdot,width=9,height=9)
ggsave(paste0(i,"geneCCdot.pdf"),CCdot,width=9,height=9)
ggsave(paste0(i,"geneMFbar.pdf"),MFbar,width=12,height=9)
ggsave(paste0(i,"geneBPbar.pdf"),BPbar,width=12,height=9)
ggsave(paste0(i,"geneCCbar.pdf"),CCbar,width=12,height=9)
ggsave(paste0(i,"geneKEGGdot.pdf"),kkdot,width=9,height=9)
ggsave(paste0(i,"geneKEGGbar.pdf"),kkbar,width=12,height=9)

#p.adjust
ggsave(paste0(i,"geneMFdot-padj.pdf"),MFdot1,width=9,height=9)
ggsave(paste0(i,"geneBPdot-padj.pdf"),BPdot1,width=9,height=9)
ggsave(paste0(i,"geneCCdot-padj.pdf"),CCdot1,width=9,height=9)
ggsave(paste0(i,"geneMFbar-padj.pdf"),MFbar1,width=12,height=9)
ggsave(paste0(i,"geneBPbar-padj.pdf"),BPbar1,width=12,height=9)
ggsave(paste0(i,"geneCCbar-padj.pdf"),CCbar1,width=12,height=9)
ggsave(paste0(i,"geneKEGGdot-padj.pdf"),kkdot1,width=9,height=9)
ggsave(paste0(i,"geneKEGGbar-padj.pdf"),kkbar1,width=12,height=9)


warnings()

