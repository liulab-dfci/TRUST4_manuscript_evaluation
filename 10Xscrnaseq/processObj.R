library("Seurat")
library("ggplot2")
library("ggrepel")

seuratResult<-readRDS("./scRNA/10X_vdj_nextgem_hs_pbmc3/10X_vdj_nextgem_gex_scRNA_Object.rds")
seuratResult$RNA@meta.data

seuratResultCD8T<-seuratResult$RNA@meta.data[which(seuratResult$RNA@meta.data["assign.ident"] == "CD8Tcells"),]

head(seuratResultCD8T)

# Only compares cd8T cell population between gdT with remaining CD8 T cell
t4Barcodes<-read.table("./scRNA/10X_vdj_nextgem_hs_pbmc3/gdT_barcode.txt", header = FALSE, sep = "")
nonT4Barcodes<-read.table("./scRNA/10X_vdj_nextgem_hs_pbmc3/CD8_notGdT_barcode.txt", header=FALSE, 
                          sep="", dec=".")

markers<-FindMarkers(object = seuratResult$RNA, ident.1=as.vector(t4Barcodes$V1), 
                     ident.2=as.vector(nonT4Barcodes$V1))
head(markers)
write.csv(markers, "./scRNA/10X_vdj_nextgem_hs_pbmc3/gdT_vs_CD8nGdT_marker.csv")

# plot the results with scatter plot
data<-markers[ markers$p_val_adj < 1, ]
sp<-ggplot(data, aes(avg_logFC, -log10(p_val_adj), label=rownames(data)))+geom_point()+geom_text_repel()
sp<-sp+xlim(-2,2)+xlab("avg_logFC (gdT / otherCD8T)")
sp
