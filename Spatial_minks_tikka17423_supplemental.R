#I am analysing 10XVisium spatial transcriptomics data of mink incisors with Seurat, MILP, and other image analysis approaches. Pauli Tikka, 17.4.2023 


#' Based on partly:
#' https://elifesciences.org/articles/81525?utm_source=content_alert&utm_medium=email&utm_content=fulltext&utm_campaign=6-February-23-elife-alert
#' https://github.com/anacaetano/human-oral-mucosa-spatial/
#' #https://github.com/anacaetano/human-oral-mucosa-spatial/blob/main/OM-spatial.R
#' #http://bioconductor.riken.jp/packages/3.9/bioc/vignettes/SpatialCPie/inst/doc/SpatialCPie.html

require(dplyr);require(ggplot2);require(magick);require(cowplot);
require(ggpubr);library(caret);library(AppliedPredictiveModeling); library(BayesSpace)

ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial'; setwd(dir=ok)# 

save.image("minks_optia.RData")
load("minkse.RData"); load("newminks.RData"); load('minks_opti.RData'); 
load("minks_optia.RData") 

# load("mips.RData") 

#if data in csc, use paultik, pwd: Ze..., to https://sd-connect.csc.fi/browse/paulitik/90bd2275c4b4413b89deead4702615ac

#Seurat workflow #or: A1_spatial <- data_mink_E14.5_1...
#Let's mark 'N1_spatial's as 'cortexes'
A1_spatial <- cortex1; 
B1_spatial <- cortex2; 
C1_spatial <- cortex3; 
D1_spatial <- cortex4 #
DefaultAssay(A1_spatial)='Spatial';DefaultAssay(B1_spatial)='Spatial';DefaultAssay(C1_spatial)='Spatial';DefaultAssay(D1_spatial)='Spatial'

p1 <- DimPlot(A1_spatial, reduction = "umap", label = TRUE)#cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink")
p2 <- SpatialDimPlot(A1_spatial, label.size = 3, pt.size.factor = 3)#cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink")

jpeg("Mink1_basics.jpg", width = 10000, height = 10000, quality = 100,pointsize = 8, res=1200);p1 + p2;dev.off()
p1 <- DimPlot(B1_spatial, reduction = "umap", label = TRUE)#cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink")
p2 <- SpatialDimPlot(B1_spatial, label.size = 3, pt.size.factor = 3)#cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink")
jpeg("Mink2_basics.jpg", width = 10000, height = 10000, quality = 100,pointsize = 8, res=1200);p1 + p2;dev.off()
p1 <- DimPlot(C1_spatial, reduction = "umap", label = TRUE,pt.size = 8)
p2 <- SpatialDimPlot(C1_spatial, label.size = 8, pt.size.factor = 3)
jpeg("Mink3_basics.jpg", width = 10000, height = 10000, quality = 100,pointsize = 8, res=1200); p1+p2; dev.off()
p1 <- DimPlot(D1_spatial, reduction = "umap", label = TRUE,pt.size = 8)
p2 <- SpatialDimPlot(D1_spatial, label = TRUE, label.size = 3,pt.size.factor = 3)
jpeg("Mink4_basics.jpg", width = 10000, height = 10000, quality = 100,pointsize = 8, res=1200);p1 + p2;dev.off()


#Highlight regions, cortexes seem to be working but not the original big data...
# SpatialDimPlot(cortex1, cells.highlight = CellsByIdentities(object = cortex1, idents = 2), facet.highlight = TRUE, ncol = 1, pt.size.factor = 2.5,  cols.highlight = c("grey50","seagreen3"))

#Identification of spatial variable features
cortex1 <- FindSpatiallyVariableFeatures(cortex1, assay = "SCT", features = VariableFeatures(cortex1)[1:1000], selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(cortex1), selection.method = "markvariogram")
SpatialFeaturePlot(cortex1, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 3.5)
head(top.features, n=20)

#Find spatial variables by region (alternative example)
##Region 0
de_markers <- FindMarkers(A1_spatial, ident.1 = 2)
SpatialFeaturePlot(object = A1_spatial, features = rownames(de_markers)[1:6], alpha = c(0.1, 1), ncol = 3, pt.size.factor = 2.5)
head(de_markers, n=50)
write.table(de_markers, file = "A1-region_epit.csv", sep = ',', quote = F, col.names = NA)

#Integration analysis
#add single-cell rna seq data #The following dataset was used in the manuscript # (Caetano & Yianni et al. 2021) GSE152042
E11=readRDS(file = "D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Corrected/E11/Integration/tot/E11.5_grant int_ptt30323.rds") #this newest one is maybe the best
DimPlot( E11, reduction = "umap",pt.size = 1.5, label = TRUE, label.size =14 ,group.by = 'seurat_clusters') + ggtitle('') + guides(colour = guide_legend(override.aes = list(size=8)))

DefaultAssay(E11) <- "SCT"
E11 <- RunPCA(E11, npcs=50,verbose = FALSE);
ElbowPlot(object = E11, ndims = 80)# https://stackoverflow.com/questions/53843589/split-variable-on-every-other-row-to-form-two-new-columns-in-data-frame
pct <- E11[["pca"]]@stdev / sum(E11[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2);co1;co2;pcs
plot_df <- data.frame(pct = pct,cumu = cumu,rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + geom_text() + geom_vline(xintercept = 99, color = "grey") + geom_hline(yintercept = min(pct[pct > 1]), color = "grey") +theme_bw()
z=6 
E11 <- RunUMAP(E11, dims = 1:z,return.model=TRUE) #seed.use = saved.seed, huom. without seed.use estimate, basically no matter, return.model=TRUE
E11 <- FindNeighbors(E11, dims = 1:z) #  epim@meta.data= epim@meta.data[, c(1:15)]
E11 <- FindClusters(E11, resolution = 0.2) #0.7/0.8/1 is good for pseudo check the final resolution from above 'hei',
cal=0:(dim(table ( Idents( E11) ))-1); cyl=1:dim(table ( Idents( E11) ));current.cluster.ids <- cal; new.cluster.ids <- cyl; #or e.g.:
E11@active.ident <- plyr::mapvalues(x =  E11@active.ident , from = current.cluster.ids, to = new.cluster.ids);
E11$seurat_clusters <- plyr::mapvalues(x =  E11$seurat_clusters , from = current.cluster.ids, to = new.cluster.ids);
DimPlot( E11, reduction = "umap",pt.size = 3.5, label = TRUE, label.size =8 ,group.by = 'seurat_clusters') +ggtitle('') + guides(colour = guide_legend(override.aes = list(size=8)))
name='data_e11_Xclust_tikka5423.rds'; saveRDS(E11, file = name)

#A1
DefaultAssay(E11)='SCT'
anchors_1 <- FindTransferAnchors(reference = E11, query = A1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors_1, refdata = E11$seurat_clusters, prediction.assay = TRUE, weight.reduction = A1_spatial[["pca"]], dims = 1:30, k.weight = 32)
A1_spatial[["predictions"]] <- predictions.assay
DefaultAssay(cortex1) <- "predictions"
SpatialFeaturePlot(cortex1, features = "1", ncol = 1, crop = TRUE)
cortex1 <- FindSpatiallyVariableFeatures(cortex1, assay = "predictions", selection.method = "markvariogram", features = rownames(cortex1), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex1), 4)
SpatialPlot(object = cortex1, features = top.clusters, ncol = 2,pt.size.factor = 5) #These 1s and 2s does not say anything, what genes they are?

#B1
anchors <- FindTransferAnchors(reference = E11, query = B1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = E11$seurat_clusters, prediction.assay = TRUE, weight.reduction = B1_spatial[["pca"]], dims = 1:30, k.weight = 27) #changed: k.weight = 32"
B1_spatial[["predictions"]] <- predictions.assay
DefaultAssay(B1_spatial) <- "predictions"
SpatialFeaturePlot(B1_spatial, features = "1", ncol = 1, crop = TRUE, pt.size.factor = 2)
B1_spatial <- FindSpatiallyVariableFeatures(B1_spatial, assay = "predictions", selection.method = "markvariogram", features = rownames(B1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(B1_spatial), 4)
SpatialPlot(object = B1_spatial, features = top.clusters, ncol = 2, pt.size.factor = 3)

#C1
anchors <- FindTransferAnchors(reference = E11, query = C1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = E11$seurat_clusters, prediction.assay = TRUE, weight.reduction = C1_spatial[["pca"]], dims = 1:30, k.weight = 28) #"was 32"
C1_spatial[["predictions"]] <- predictions.assay
DefaultAssay(C1_spatial) <- "predictions"
SpatialFeaturePlot(C1_spatial, features = "1", ncol = 1, crop = TRUE, pt.size.factor = 3)
C1_spatial <- FindSpatiallyVariableFeatures(C1_spatial, assay = "predictions", selection.method = "markvariogram",features = rownames(C1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(C1_spatial), 4)
SpatialPlot(object = C1_spatial, features = top.clusters, ncol = 2, pt.size.factor = 3)

#D1
D1_spatial <- SCTransform(D1_spatial, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = E11, query = D1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = E11$seurat_clusters, prediction.assay = TRUE, weight.reduction = D1_spatial[["pca"]], dims = 1:30, k.weight = 23) #32
D1_spatial[["predictions"]] <- predictions.assay
DefaultAssay(D1_spatial) <- "predictions"
SpatialFeaturePlot(D1_spatial, features = "1", ncol = 1, crop = TRUE)
D1_spatial <- FindSpatiallyVariableFeatures(D1_spatial, assay = "predictions", selection.method = "markvariogram", features = rownames(D1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(D1_spatial), 4)
SpatialPlot(object = D1_spatial, features = top.clusters, ncol = 2,pt.size.factor = 5)

#BayesSpace
#Preprocessing data
set.seed(102)
DefaultAssay(cortex1)='Spatial'
sce_A1=as.SingleCellExperiment(cortex1)
colnames(sce_A1)=paste("Spot",1:length(colnames(sce_A1)),sep="_")
reducedDims(sce_A1)$PCA=reducedDims(sce_A1)$PCA[1:length(colnames(sce_A1)),]
sce_A1@colData$row=-cortex1@images[["slice1"]]@coordinates[["row"]]
sce_A1@colData$col=-cortex1@images[["slice1"]]@coordinates[["col"]]
sce_A1@colData$imagerow=-cortex1@images[["slice1"]]@coordinates[["imagerow"]]
sce_A1@colData$imagecol=-cortex1@images[["slice1"]]@coordinates[["imagecol"]]

DefaultAssay(cortex2)='Spatial'
sce_B1=as.SingleCellExperiment(cortex2)
colnames(sce_B1)=paste("Spot",1:length(colnames(sce_B1)),sep="_")
sce_B1@colData$row=-cortex2@images[["slice1"]]@coordinates[["row"]]
sce_B1@colData$col=-cortex2@images[["slice1"]]@coordinates[["col"]]
sce_B1@colData$imagerow=-cortex2@images[["slice1"]]@coordinates[["imagerow"]]
sce_B1@colData$imagecol=-cortex2@images[["slice1"]]@coordinates[["imagecol"]]


DefaultAssay(C1_spatial)='Spatial' #or cortex3
sce_C1=as.SingleCellExperiment(C1_spatial)
colnames(sce_C1)=paste("Spot",1:length(colnames(C1_spatial)),sep="_")
sce_C1@colData$row=-C1_spatial@images[["slice1"]]@coordinates[["row"]]
sce_C1@colData$col=-C1_spatial@images[["slice1"]]@coordinates[["col"]]
sce_C1@colData$imagerow=-C1_spatial@images[["slice1"]]@coordinates[["imagerow"]]
sce_C1@colData$imagecol=-C1_spatial@images[["slice1"]]@coordinates[["imagecol"]]
DefaultAssay(C1_spatial)='Spatial'

sce_D1=as.SingleCellExperiment(cortex4)
colnames(sce_D1)=paste("Spot",1:length(colnames(sce_D1)),sep="_")
sce_D1@colData$row=-cortex4@images[["slice1"]]@coordinates[["row"]]
sce_D1@colData$col=-cortex4@images[["slice1"]]@coordinates[["col"]]
sce_D1@colData$imagerow=-cortex4@images[["slice1"]]@coordinates[["imagerow"]]
sce_D1@colData$imagecol=-cortex4@images[["slice1"]]@coordinates[["imagecol"]]

# euna=make.names(colnames(sce_A1), unique=TRUE); n1=euna;n11=substr(n1, 1, 16); # colnames(sce_A1)=euna ;; # colnames(sce_A1) <- paste(colnames(sce_A1),"new",sep="_")  
# sce@colData #pitäis olla solujen nimien sijaan spot_N ja spatial cluster ja..:row       col  imagerow  imagecol sizeFactor spatial.cluster
sce_A1 <- spatialPreprocess(sce_A1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)
sce_B1 <- spatialPreprocess(sce_B1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)
sce_C1 <- spatialPreprocess(sce_C1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)
sce_D1 <- spatialPreprocess(sce_D1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)


#Clustering
## understand number of clusters, https://edward130603.github.io/BayesSpace/articles/BayesSpace.html
set.seed(149)
sce_A1 <- qTune(sce_A1, qs=seq(2, 15),platform="Visium",d=6);qPlot(sce_A1) #test, but not ok: sce_A1@colData$spatial.cluster=sce_A1@colData$ident but the rows and columns did,jei! :)
sce_B1 <- qTune(sce_B1, qs=seq(2, 15), platform="Visium", d=6);qPlot(sce_B1)
sce_C1 <- qTune(sce_C1, qs=seq(2, 15), platform="Visium", d=6); qPlot(sce_C1)
sce_D1 <- qTune(sce_D1, qs=seq(2, 15), platform="Visium", d=6); qPlot(sce_D1)

## clustering
sce_A1 <- spatialCluster(sce_A1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE); head(colData(sce_A1))
sce_B1 <- spatialCluster(sce_B1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE);head(colData(sce_B1))
sce_C1 <- spatialCluster(sce_C1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE);head(colData(sce_C1))
sce_D1 <- spatialCluster(sce_D1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE);head(colData(sce_D1))

#visualising spatial clusters

ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial/devio'; setwd(dir=ok)#

clusterPlot(sce_A1); 
clusterPlot(sce_B1)
clusterPlot(sce_C1); 
clusterPlot(sce_D1)

jpeg('sce_A1_c.jpg', width =10000 , height = 10000, quality = 100,pointsize = 12, res=1000);clusterPlot(sce_A1);dev.off()
jpeg('sce_B1_c.jpg', width =10000 , height = 10000, quality = 100,pointsize = 12, res=1000);clusterPlot(sce_B1);dev.off()
jpeg('sce_C1_mod.jpg', width =10000 , height = 10000, quality = 100,pointsize = 12, res=1000);clusterPlot(sce_C1);dev.off()
jpeg('sce_D1_c.jpg', width =10000 , height = 10000, quality = 100,pointsize = 12, res=1000);clusterPlot(sce_D1);dev.off()


# https://lmweber.org/OSTA-book/image-segmentation-visium.html#identify-number-of-cells-per-spot
# https://lmweber.org/OSTA-book/space-ranger-visium.html#space-ranger-visium
#http://research.libd.org/VistoSeg/
#http://research.libd.org/VistoSeg/step-4-gui-to-count-nuclei-in-a-visium-spot.html

## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- CreateSeuratObject(counts=logcounts(sce_C1),assay='Spatial',meta.data=as.data.frame(colData(sce_C1)))
sobj <- SetIdent(sobj, value = "spatial.cluster")
sobj@assays$Spatial@scale.data <- sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
## Select top n markers from each cluster (by log fold change)
top_markers <- FindAllMarkers(sobj, assay='Spatial', slot='data',group.by='spatial.cluster',only.pos=FALSE,test.use = "wilcox",logfc.threshold = 0.01,min.pct = 0.01) 
top_markers %>% dplyr::group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC) -> top7#added the class, https://stackoverflow.com/questions/65008837/group-by-using-the-function-for-didn%C2%B4t-work-r
jpeg('sce_C1_heatmap.jpg', width =10000 , height = 10000, quality = 100,pointsize = 12, res=1000);
DoHeatmap(sobj, features = top7$gene, slot='scale.data', group.by = "spatial.cluster",
          angle=0, size=6, label = TRUE, raster=FALSE,group.bar.height = 0.035)+ guides(scale = "none") + 
  theme(text = element_text(size = 15));dev.off()#+ guides(shape = guide_legend(override.aes = list(size=20)))
#https://github.com/satijalab/seurat/issues/1437 
top_markers['Sox2',]; top_markers['Pitx2',]

sobj <- CreateSeuratObject(counts=logcounts(sce_D1),assay='Spatial',meta.data=as.data.frame(colData(sce_D1)))
sobj <- SetIdent(sobj, value = "spatial.cluster")
sobj@assays$Spatial@scale.data <- sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
## Select top n markers from each cluster (by log fold change)
top_markers <- FindAllMarkers(sobj, assay='Spatial', slot='data',group.by='spatial.cluster',only.pos=FALSE,test.use = "wilcox",logfc.threshold = 0.01,min.pct = 0.01) 
top_markers %>% dplyr::group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) -> top7#added the class, https://stackoverflow.com/questions/65008837/group-by-using-the-function-for-didn%C2%B4t-work-r
jpeg('sce_C1_heatmap_30.jpg', width =10000 , height = 15000, quality = 100,pointsize = 12, res=1000);
DoHeatmap(sobj, features = top7$gene, slot='scale.data', group.by = "spatial.cluster",angle=0, size=6, label = TRUE, raster=FALSE,group.bar.height = 0.035)+ guides(scale = "none") + 
  theme(text = element_text(size = 15));dev.off()#+ guides(shape = guide_legend(override.aes = list(size=20)))
#https://github.com/satijalab/seurat/issues/1437 
top_markers['Sox2',]
top_markers['Pitx2',]

sce_A1 <- spatialCluster(sce_A1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE); head(colData(sce_A1))
sce_B1 <- spatialCluster(sce_B1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE); head(colData(sce_B1))
sce_C1 <- spatialCluster(sce_C1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE); head(colData(sce_C1))
sce_D1 <- spatialCluster(sce_D1, q=4, platform = "Visium", d=6, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE); head(colData(sce_D1))

markersu=c('Pitx2','Sox2','Acta2', 'Bmi1', 'Igfbp5',  
           'Lef1', 'Wnt10a', 'Wnt10b', 'Wnt4',  'Cdh5',  
           'Notch1', 'Notch2', 'Gli1', 'Ptch1',  'Sostdc1', 
           'Bmp4',   'Epcam', 'Msx1', 'Msx2','Vim')#'Krt14','Shh','Fgf4','Fgf8','Lgr5','Cdh1',

#Mink1 specs:
# draw_image(img,scale = 0.895, x = -0.02,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800)
# ggdraw() +draw_image(img,scale = 0.88, x = -0.033,y=0.015) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800)
# ggsave(paste("Mink_S1_background_Pitx2_all_topmark_wil",name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")
#Mink2 specs:
# ggdraw() +draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800)
# ggsave(paste("Mink_S2_background_Pitx2_wil",name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")}
# Mink3 specs:
# ggdraw() +draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) # best so far, except resolution+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
# name=featuresh[i]; ggdraw() +draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) # best so far, 13, scale = 0.7)
# ggsave(paste("Mink_S3_background_Pitx2_all_topmark_auc",name,".png"),width = 8000,height = 5000,units = "px",dpi = 1000, bg = "white")}
# Mink4 specs:
# ggdraw() +draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800)
ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial'; setwd(dir=ok)
img <- image_read('Mink3_mand.jpg')

ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial/Mink3/h_norm'; setwd(dir=ok)
#The best so far:# ggdraw() +draw_image(img2,x = 1,scale = 1.65) +featurePlot(sce_C1,'Sox2',alpha=0.5,low='#C1CDCD',high='red')+theme(aspect.ratio = 800 / 800)
p=BayesSpace::featurePlot(sce_D1,'Pitx2',alpha=0.4,low='#C1CDCD',high='red') #+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
ggdraw() +
# draw_image(img,scale = 0.88, x = -0.029,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
# draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
# draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) # 4, best so far, except resolution+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
ggsave("Mink_S3_background_Pitx2_oke.png",width = 6000,height = 6000,units = "px",dpi = 800,bg = "white")

sce_test=sce_A1
# sum(sce_test['Pitx2',]@assays@data[['logcounts']]>2)
# ho=which(sce_test['Pitx2',]@assays@data[['logcounts']]>2.0)
# sce_test_6=sce_test[,ho]

#Most epithelial:
library(pathviewr)
#Mesenchymal regression:
sobje <- CreateSeuratObject(counts=logcounts(sce_test),assay='Spatial',meta.data=as.data.frame(colData(sce_test)))
sobje <- SetIdent(sobje, value = "spatial.cluster")
sobje@assays$Spatial@scale.data <- sobje@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
# DefaultAssay( sobje) <- "SCT" #this is important! #epix<- subset(x = E11, idents = c(4,5,6,7,8,9,11,12)) #3,4,5,6,7,8
GAD1=GetAssayData(object = sobje) #dim(GAD1) # > dim(GAD1) # [1] 16154    6970
GAD11=rowMeans(sobje) #GAD1[rownames(GAD1)=='Sox2',]  #length(GAD11) is 3800, which(rownames(GAD1)=='Sox2'), 3773:which(names(GAD11[rev(order(GAD11))])=='Sox2') -> ok to find out sox2
fet=names(GAD11[(order(GAD11))][length(GAD11):(length(GAD11)-55)])
jpeg("Top55_S3_Heatmap_high.jpg", width = 10000, height = 15000, quality = 100,pointsize = 16, res=1200);
DoHeatmap(sobje, features = fet, size = 3, slot='data') +labs(color = "Cluster");dev.off() 

r1=c('Xist','Tshz2','Crabp1','Nusap1','Nap1l1','Dcn','Hmcn1','Serbp1'); #'Krt18'
r2 <- rownames(sce_test)[rownames(sce_test) %in% rownames(sce_test)[grep("*Rik",rownames(sce_test))]];
r3 <- rownames(sce_test)[rownames(sce_test) %in% rownames(sce_test)[grep("^Mir",rownames(sce_test))]];
r4=rownames(sce_test)[rownames(sce_test) %in% rownames(sce_test)[grep("^Col",rownames(sce_test))]];
r5=c('Vim','Wnt5a','Acta1','Prdx1') #index out of bounds: Phox1 (i.e. Prxx1) Gapdh Hrrt Hbb
lim1=ceiling(quantile(GAD11[(order(GAD11))][length(GAD11):(length(GAD11)-300)],0.85)); e=GAD11[(order(GAD11))][length(GAD11):(length(GAD11)-300)];
r6=names(e[e>lim1])
# mesg=unique(c(r1,r2,r3,r4,r5,r6))
mesg=unique(c(r2,r3,r4,r5))
C<-GetAssayData(object = sobje, slot = "counts"); mesg=rownames(C[rownames(C) %in% mesg,]); mesga <- colSums(C[mesg,])/Matrix::colSums(C)*100
sobje <- AddMetaData(sobje, mesga, col.name = "percent_mes") 
hist(sobje@meta.data$percent_mes,breaks=10)

sce_test
epig=c('Pitx2','Sox2','Epcam','Krt15')
D_ep=colMeans(sce_test[epig,]@assays@data[['logcounts']])/colMeans(sce_test[mesg,]@assays@data[['logcounts']])
normalized = (D_ep-min(D_ep))/(max(D_ep)-min(D_ep))
quantile(normalized,0.95); hist(as.vector(normalized),breaks=50)
df=data.frame(normalized, stringsAsFactors = TRUE)
# row.names(df)='N(1-(Vim-Pitx2)/Vim)'
df2a=df %>% filter(if_all(everything(), ~ . >0.2)) #this did not work:df[which(normalized>0.25),], but https://stackoverflow.com/questions/45427444/how-to-preserve-base-data-frame-rownames-upon-filtering-in-dplyr-chain
sce_test_ep=sce_test[,rownames(df2a)]
exprs_mat = SummarizedExperiment::assay(sce_test_ep, 'logcounts')
set.seed(19)
result2 = scSEGIndex(exprs_mat = exprs_mat)
#lower the better, let's take 200 and compare that to tgg
result2
df <- data.frame(x = seq(1:1000),y = result2[with(result2,order(-segIdx)),][1:1000,'segIdx'])
plot(df)
ae=find_curve_elbow(df, plot_curve = TRUE);ae #:) [1] 130
comp2=result2[with(result2,order(-segIdx)),][1:ae,]
cp2=rownames(comp2)

D_mes=colMeans(sce_test[mesg,]@assays@data[['logcounts']])/(colMeans(sce_test[epig,]@assays@data[['logcounts']])+quantile(colMeans(sce_test[epig,]@assays@data[['logcounts']]),0.25)/2+0.01)
normalized = (D_mes-min(D_mes))/(max(D_mes)-min(D_mes))
quantile(normalized,0.95); hist(as.vector(normalized),breaks=50)
df=data.frame(normalized, stringsAsFactors = TRUE)
# row.names(df)='N(1-(Vim-Pitx2)/Vim)'
df2b=df %>% filter(if_all(everything(), ~ . >0.5)) #this did not work:df[which(normalized>0.25),], but https://stackoverflow.com/questions/45427444/how-to-preserve-base-data-frame-rownames-upon-filtering-in-dplyr-chain
sce_test_mes=sce_test[,rownames(df2b)]
exprs_mat = SummarizedExperiment::assay(sce_test_mes, 'logcounts');set.seed(19)
result3 = scSEGIndex(exprs_mat = exprs_mat)
#lower the better, let's take 200 and compare that to tgg
result3
df <- data.frame(x = seq(1:500),y = result3[with(result3,order(-segIdx)),][1:500,'segIdx'])
plot(df)
ae=find_curve_elbow(df, plot_curve = TRUE) #:) [1] 130
comp3=result3[with(result3,order(-segIdx)),][1:ae,]
cp3=rownames(comp3)
# sum(sce_test['Pitx2',]@assays@data[['logcounts']]>2)
# ho=which(sce_test['Pitx2',]@assays@data[['logcounts']]>2.0)
# sce_test_ep=sce_test[,ho]
sce_test_ep[cp2,]
hist(colMeans(sce_test_ep[cp2,]@assays@data[['logcounts']]))
sce_test_mes[cp3,]

# ho=which(colnames(sce_test) %in% colnames(df2))
p=BayesSpace::featurePlot(sce_test_6,'Pitx2',alpha=0.4,low='#C1CDCD',high='red') #+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
ggdraw() +
  # draw_image(img,scale = 0.92, x = -0.02,y=-0.065) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
  # draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
  draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
  # draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) # 4, best so far, except resolution+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
ggsave("Mink_S3_Epithelial Spots as Pitx2&Vim Calc.png",width = 6000,height = 6000,units = "px",dpi = 800,bg = "white")


#Identification of spatial variable features
sobje <- CreateSeuratObject(counts=logcounts(sce_test_6),assay='Spatial',meta.data=as.data.frame(colData(sce_test_6)))
sobje <- SetIdent(sobje, value = "spatial.cluster")
sobje@assays$Spatial@scale.data <- sobje@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
# DefaultAssay( sobje) <- "SCT" #this is important! #epix<- subset(x = E11, idents = c(4,5,6,7,8,9,11,12)) #3,4,5,6,7,8
GAD1=GetAssayData(object = sobje) #dim(GAD1) # > dim(GAD1) # [1] 16154    6970
GAD11=rowMeans(sobje) #GAD1[rownames(GAD1)=='Sox2',]  #length(GAD11) is 3800, which(rownames(GAD1)=='Sox2'), 3773:which(names(GAD11[rev(order(GAD11))])=='Sox2') -> ok to find out sox2
featuresh=names(GAD11[(order(GAD11))][length(GAD11):(length(GAD11)-16)])
jpeg("Top16_S3_Heatmap_Pitx2_high.jpg", width = 12000, height = 9000, quality = 100,pointsize = 16, res=1200);
DoHeatmap(sobje, features = featuresh, size = 3, slot='data') +labs(color = "Cluster");dev.off() 

sobjee <- CreateSeuratObject(counts=logcounts(sce_test),assay='Spatial',meta.data=as.data.frame(colData(sce_test)))
sobjee <- SetIdent(sobjee, value = "spatial.cluster")
sobjee@assays$Spatial@scale.data <- sobjee@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
sobjee@meta.data$ident=1
sobjee@meta.data[colnames(df2),'ident']=2 #or ho
top_markers <- FindMarkers(sobjee, assay='Spatial', slot='data',group.by='ident',ident.1 = 2,only.pos=FALSE,test.use = "wilcox",logfc.threshold = 0.01,min.pct = 0.01) 
top_markers2 <- FindMarkers(sobjee, assay='Spatial', slot='data',group.by='ident',ident.1 = 2,only.pos=FALSE,test.use = "roc",logfc.threshold = 0.01,min.pct = 0.01) 
featureshe1=rownames(top_markers[top_markers[,'p_val_adj']<0.05,]);featureshe1 #hope this is more than 10
featureshe2=rownames(top_markers2[top_markers2[,'myAUC']>0.87,]);featureshe2 #85% or 70% etc. so as to get about 20

jpeg("Top_DE_wilcox3_S3_Heatmap_Pitx2_high.jpg", width = 12000, height = 9000, quality = 100,pointsize = 16, res=1200);
DoHeatmap(sobje, features = featureshe1, size = 3, slot='data') +labs(color = "Cluster");dev.off() 

# sce <- exampleSCE()
# BayesSpace::featurePlot(sce_test, "Pitx2")#Why I now need to add BayesSpace? 15323

for (i in 1:length(featuresh)) {
  p=BayesSpace::featurePlot(sce_test,featuresh[i],alpha=0.4,low='#C1CDCD',high='red')#+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
  name=featuresh[i]; ggdraw() +
    # draw_image(img,scale = 0.88, x = -0.029,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
  draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
  # draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
  # draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) # best so far, 13, scale = 0.7)
  ggsave(paste("Mink_S3_background_Pitx2_high",name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")}

for (i in 1:length(featureshe1)) {
    p=BayesSpace::featurePlot(sce_test,featureshe1[i],alpha=0.4,low='#C1CDCD',high='red')#+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
    name=featureshe1[i]; ggdraw() +
      # draw_image(img,scale = 0.88, x = -0.029,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
    draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
    # draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
    # draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800)  # best so far, 13, scale = 0.7)
    ggsave(paste("Mink_S3_background_Pitx2_wil",name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")}
  
for (i in 1:length(featureshe2)) {
    p=BayesSpace::featurePlot(sce_test,featureshe2[i],alpha=0.4,low='#C1CDCD',high='red')#+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
    name=featureshe2[i]; ggdraw() +
      # draw_image(img,scale = 0.88, x = -0.029,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
    draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
    # draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
    # draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800)# best so far, 13, scale = 0.7)
    ggsave(paste("Mink_S3_background_Pitx2_auc",name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")}

ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial'; setwd(dir=ok); img <- image_read('Mink4_mand.jpg')
ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial/Mink4/h3/basic'; 
setwd(dir=ok); sce_test=sce_D1; namea="Mink_S3_background_"
markersu=c('Pitx2','Sox2','Acta2', 'Bmi1', 'Igfbp5',  'Lef1', 'Wnt10a', 'Wnt10b', 'Wnt4',  'Cdh5',  'Fgf8',
'Notch1', 'Notch2', 'Gli1', 'Ptch1',  'Sostdc1', 'Bmp4',   'Epcam', 'Msx1', 'Msx2','Vim')#not in slice4:'Krt14','Shh','Fgf4','Lgr5','Cdh1', #Krt14 also not...

for (i in 9:length(markersu)) {
  p=featurePlot(sce_test,markersu[i],alpha=0.4,low='#C1CDCD',high='red')#+geom_point(alpha = 0.6)# draw_plot(p, x = -0.1, y = -0.13, scale = 0.7)
  name=markersu[i]; ggdraw() +  
  # draw_image(img,scale = 0.88, x = -0.029,y=0.012) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #1
  # draw_image(img,scale = 1.01, x = -0.003,y=0.003) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #
  # draw_image(img,scale = 1, x = -0.022,y=0.001) +draw_plot(p, scale = 1)+theme(aspect.ratio = 800 / 800) #3
  draw_image(img,scale = 0.99, x = -0.02,y=0.013) +draw_plot(p, scale = 1.0,x = 0.008)+theme(aspect.ratio = 800 / 800) #4
  ggsave(paste(namea,name,".png"),width = 6000,height = 6000,units = "px",dpi = 800, bg = "white")}

# f1=featuresh; f2=featuresh; f3=  featuresh; f4=  featuresh
length(unique(c(f1,f2,f3,f4)))
length(Reduce(intersect, list(f1,f2,f3,f4))) #60% same

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/#change-columnrow-span-of-a-plot
 #https://satijalab.org/seurat/articles/spatial_vignette.html

#Enhanced gene expression
# markers <- list(); markers[["Basal keratinocyte"]] <- c("KRT5", "KRT14", "COL7A1", "CXCL14"); markers[["Suprabasal keratinocyte"]] <- c("SBSN", "SPRR1B")
# sce_A1.enhanced <- enhanceFeatures(sce_A1.enhanced, sce_A1,model="xgboost",feature_names=purrr::reduce(markers, c),nrounds=0, eta=1, max.depth=2, colsample_bytree=0.7, colsample_bylevel=0.8)
# sum_counts <- function(sce_A1, features) {if (length(features) > 1) {colSums(logcounts(sce_A1)[features, ]) } else {logcounts(sce_A1)[features, ]}}

#Giotto ligand-receptor analysis
# remotes::install_github("RubD/Giotto") #https://rubd.github.io/Giotto_site/
results_folder = 'D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial/Mink1/res'
python_path = 'C:/Hyapp/Anaconda' #You may need to test this with laptop, since here it is not working
## provide path to visium folder
data_path = 'G:/scRNAseq_minks/mink_E14-5_1/outs/'

## create instructions
instrs = createGiottoInstructions(save_dir = results_folder,save_plot = TRUE,show_plot = TRUE,python_path = python_path)
## Create Giotto object & Process Data
## directly from visium folder
 sc <- spark_connect(master = "local") #spark_disconnect(sc) https://spark.rstudio.com/get-started/
A1 = createGiottoVisiumObject(visium_dir = data_path, expr_data = 'raw',png_name = 'tissue_lowres_image.png',gene_column_index = 2, instructions = instrs) # smfishHmrf trendsceek SPARK multinet RTriangle FactoMiner

## update and align background image
# problem: image is not perfectly aligned
spatPlot(gobject = A1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,save_param = list(save_name = '2_a_spatplot_image')) #nyt alko toimia :) lienee toi python miniconda installointi: https://giottosuite.readthedocs.io/en/master/gettingstarted.html#part2-python-giotto-requirements

showGiottoImageNames(A1)
A1 = updateGiottoImage(A1, image_name = 'image',xmax_adj = 800, xmin_adj =600, ymax_adj = 600, ymin_adj = 600)
spatPlot(gobject = A1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,save_param = list(save_name = '2_b_spatplot_image_adjusted'))

## check metadata
pDataDT(A1)
## compare in tissue with provided jpg
spatPlot(gobject = A1, cell_color = 'in_tissue', point_size = 2,cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),save_param = list(save_name = '2_c_in_tissue'))

## subset on spots that were covered by tissue
metadata = pDataDT(A1)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
A1 = subsetGiotto(A1, cell_ids = in_tissue_barcodes)

## filter
A1 <- filterGiotto(gobject = A1,expression_threshold = 0.5,gene_det_in_min_cells = 50,
                   min_det_genes_per_cell = 10,
                   expression_values = c('raw'),
                   verbose = T)

## normalize
A1 <- normalizeGiotto(gobject = A1, scalefactor = 6000, verbose = T)

## add gene & cell statistics
A1 <- addStatistics(gobject = A1)

## visualize
spatPlot2D(gobject = A1, show_image = T, point_alpha = 0.7,save_param = list(save_name = '2_d_spatial_locations'))
spatPlot2D(gobject = A1, show_image = T, point_alpha = 0.7,cell_color = 'nr_genes', color_as_factor = F,save_param = list(save_name = '2_e_nr_genes'))
spatDimPlot(A1, cell_color = 'leiden_clus',plot_alignment = 'horizontal', spat_point_size = 2)

## Dimension Reduction
## highly variable genes (HVG)
A1 <- calculateHVG(gobject = A1,save_param = list(save_name = '3_a_HVGplot')) 


## run PCA on expression values (default)
A1 <- runPCA(gobject = A1, center = TRUE, scale_unit = T)
screePlot(A1, ncp = 30, save_param = list(save_name = '3_b_screeplot'))
plotPCA(gobject =A1,save_param = list(save_name = '3_c_PCA_reduction'))

#Datan kaivuu:
pca_estudio=A1@dimension_reduction$cells$pca[[1]][3][[1]][,1:2]
plot(pca_estudio)

## run UMAP and tSNE on PCA space (default)
##UMAP
A1 <- runUMAP(A1, dimensions_to_use = 1:8,return.model=TRUE)
plotUMAP(gobject = A1,save_param = list(save_name = '3_d_UMAP_reduction'))
##tSNE
A1 <- runtSNE(A1, dimensions_to_use = 1:8)
plotTSNE(gobject = A1,save_param = list(save_name = '3_e_tSNE_reduction'))

##Clustering
## sNN network (default)
A1 <- createNearestNetwork(gobject = A1, dimensions_to_use = 1:8, k = 15)
## Leiden clustering
A1 <- doLeidenCluster(gobject = A1, resolution = 0.35)
plotUMAP(gobject = A1,cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,save_param = list(save_name = '4_a_UMAP_leiden'))
A1 <- doLouvainCluster(gobject = A1, resolution = 3)
plotUMAP(gobject = A1,cell_color = 'louvain_clus', show_NN_network = T, point_size = 2.5,save_param = list(save_name = '4_a_UMAP_louvain'))

##Spatial Grid
A1 <- createSpatialGrid(gobject = A1, sdimx_stepsize = 400,sdimy_stepsize = 400,minimum_padding = 0)
spatPlot(A1,  show_grid = T,grid_color = 'red', spatial_grid_name = 'spatial_grid', save_param = c(save_name = '8_grid')) #cell_color = 'cell_types',

##Spatial Network
## Delaunay network: stats + creation
plotStatDelaunayNetwork(gobject = A1, maximum_distance = 400, save_param = c(save_name = '9_a_delaunay_network'))
A1 = createSpatialNetwork(gobject = A1, minimum_k = 0)
showNetworks(A1); spatPlot(gobject = A1, show_network = T,network_color = 'blue', spatial_network_name = 'Delaunay_network',save_param = c(save_name = '9_b_delaunay_network'))

##Spatial Genes
## kmeans binarization
kmtest = binSpect(A1)
spatGenePlot(A1, expression_values = 'scaled',genes = kmtest$genes[1:20], cow_n_col = 4, point_size = 1.5,save_param = c(save_name = 'spatial_genes_km'))
spatGenePlot(A1, expression_values = 'scaled',genes = kmtest$genes[21:40], cow_n_col = 4, point_size = 1.5,save_param = c(save_name = 'spatial_genes_km2'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = kmtest$genes[41:60], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_km3'))

##Rank Binarisation
ranktest = binSpect(A1, bin_method = 'rank')
spatGenePlot(A1, expression_values = 'scaled',genes = ranktest$genes[1:20], cow_n_col = 4, point_size = 1.5,save_param = c(save_name = 'spatial_genes_rank'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = ranktest$genes[21:40], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_rank2'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = ranktest$genes[41:60], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_rank3'))

# https://rubd.github.io/Giotto_site/articles/subset_giotto.html

## Cell Neighborhood analysis
cell_proximities = cellProximityEnrichment(gobject = A1,cluster_column = 'leiden_clus',spatial_network_name = 'Delaunay_network',adjust_method = 'fdr',number_of_simulations = 1000)  #'cell_types' as leiden_clu

## barplot
cellProximityBarplot(gobject = A1,CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5, save_param = c(save_name = '12_a_barplot_cell_cell_enrichment'))

## heatmap
cellProximityHeatmap(gobject = A1, CPscore = cell_proximities, order_cell_types = T, scale = T, color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))

## network
cellProximityNetwork(gobject = A1, CPscore = cell_proximities, remove_self_edges = T,only_show_enrichment_edges = F,save_param = c(save_name = '12_c_network_cell_cell_enrichment'))

## network with self-edges
cellProximityNetwork(gobject = A1, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5),
                     save_param = c(save_name = '12_d_network_cell_cell_enrichment_self',
                                    base_height = 5, base_width = 5, save_format = 'pdf'))

## select top 25th highest expressing genes
gene_metadata = fDataDT(A1)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr_det)

quantile(gene_metadata$mean_expr_det)
high_expressed_genes = gene_metadata[mean_expr_det > 1.31]$gene_ID

# https://rubd.github.io/Giotto_site/articles/mouse_visium_kidney_200916.html
scran_markers_subclusters = findMarkers_one_vs_all(gobject = A1,method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')

topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes


# violinplot
violinPlot(A1, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = c(save_name = '6_d_violinplot_scran', base_width = 5))

A1 <- createSpatialGrid(gobject = A1,
                                   sdimx_stepsize = 400,
                                   sdimy_stepsize = 400,
                                   minimum_padding = 0)

spatPlot(A1, cell_color = 'leiden_clus', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid', 
         save_param = c(save_name = '8_grid'))

spatPlot(gobject = A1, show_network = T,
         network_color = 'blue', spatial_network_name = 'distance_network',
         point_size = 2.5, cell_color = 'leiden_clus',
         save_param = list(save_name = 'distance_network_spatPlot'))

# LR expression
C<-GetAssayData(object = epim, slot = "counts"); helloa='Sox2'; Sox2_score <- colSums(C[helloa,])/Matrix::colSums(C)*100
epim <- AddMetaData(epim, Sox2_score, col.name = "Sox2") 
epim[["Sox2"]] <- PercentageFeatureSet(epim, pattern = "^Sox2")
nozeromydata <- epim[["Sox2"]][ epim[["Sox2"]] != 0 ]
m=median(nozeromydata )
epim[["Sox2_ok"]]=epim[["Sox2"]][,1]>m

DimPlot(epim, group.by = "Sox2") 

epim@meta.data$celltype_aggregate = paste(epim@meta.data$seurat_clusters, epim@meta.data$Sox2,sep = "_") # user adaptation required on own dataset

DimPlot(epim, group.by = "Sox2_ok",pt.size = 5)

table(epim@meta.data$seurat_clusters, epim@meta.data$Sox2_ok)

epim@meta.data$celltype_aggregate = paste(epim@meta.data$seurat_clusters, epim@meta.data$Sox2_ok,sep = "_") # user adaptation required on own dataset
DimPlot(epim, group.by = "celltype_aggregate")

epim@meta.data$celltype_aggregate %>% table() %>% sort(decreasing = TRUE)
celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
epim = SetIdent(epim, value = epim[[celltype_id]])
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
organism = "mouse" # user adaptation required on own dataset
if(organism == "mouse"){
  lr_network = lr_network %>% mutate(ligand = convert_human_to_mouse_symbols(ligand), receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()
  
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}

niches = list(
  "Sox2_High_niche" = list(
    "sender" = c('3_TRUE','10_TRUE', '11_TRUE'),
    "receiver" = c('7_TRUE')),
  "Sox2_Low_niche" = list(
    "sender" = c('1_FALSE', '2_FALSE', '9_FALSE'),
    "receiver" = c('5_FALSE'))) 

assay_oi = "SCT" # other possibilities: RNA,...
DefaultAssay(epim)='SCT'
DE_sender = calculate_niche_de(seurat_obj  = epim %>% subset(features = lr_network$ligand %>% unique()), 
                               niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
DE_receiver = calculate_niche_de(seurat_obj = epim %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

expression_pct = 0.10
DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)

# https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet_pEMT.md
#Jääty:
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
specificity_score_spatial = "lfc"


# LR activity changes 
#https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md
#https://github.com/saeyslab/nichenetr/issues/118
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2#Sec10

nozeromydata2 <- DE_sender_receiver$ligand_score[ DE_sender_receiver$ligand_score != 0 ]
me=median(nozeromydata2 )

median(DE_sender_receiver$ligand)

# LR_data = data.table::fread('D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial/Mink1/res')
# LR_data
# 
LR_data[, ligand_det := ifelse(ligand %in% A1@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(receptor %in% A1@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]

lig=quantile(DE_sender_receiver$scaled_ligand_score)[4]
rcp=quantile(DE_sender_receiver$scaled_receptor_score)[4]

lig_rcp=quantile(DE_sender_receiver$scaled_avg_score_ligand_receptor)[4]

as.list(DataFrame(DE_sender_receiver[DE_sender_receiver[,'scaled_ligand_score']>lig,'ligand']))
DE_sender_receiver[DE_sender_receiver[,'scaled_receptor_score']>rcp,'receptor ']

select_ligands=names(as.list(DE_sender_receiver[DE_sender_receiver[,'scaled_avg_score_ligand_receptor']>lig_rcp,'ligand'])$ligand)

# select_ligands = names(as.list(DE_sender_receiver[DE_sender_receiver[,'scaled_ligand_score']>lig,'ligand'])$ligand)
  # LR_data_det$ligand
select_receptors = names(as.list(DE_sender_receiver[DE_sender_receiver[,'scaled_avg_score_ligand_receptor']>lig_rcp,'receptor'])$receptor)
  # names(as.list(DataFrame(DE_sender_receiver[DE_sender_receiver[,'scaled_receptor_score']>rcp,'receptor'])$receptor))
  # LR_data_det$receptor

tv=A1@norm_scaled_expr
iog=DataFrame(rownames(tv))
ig2=as.list(iog)$rownames.tv.
# tve=tv[ig2 %in% select_ligands,]
v1e=unique(select_ligands)
vie1=v1e[!v1e %in% ig2]
v1e=v1e[!v1e %in% vie1]

v2e=unique(select_receptors)
vie2=v2e[!v2e %in% ig2]
v2e=v2e[!v2e %in% vie2]

# v2e= unique(select_receptors)
# sum(ig2 %in% v1e)

limita=DE_sender_receiver[DE_sender_receiver[,'scaled_avg_score_ligand_receptor']>lig_rcp,c('ligand','receptor')]
cond=names(as.list(limita[,'ligand'])$ligand) %in% v1e
cond2=names(as.list(limita[,'receptor'])$receptor) %in% v2e

c3= cond & cond2 #the right condition is to have both with &

select_ligands=names(as.list(limita[c3,'ligand'])$ligand)
select_receptors=names(as.list(limita[c3,'receptor'])$receptor)

## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = A1,
                                   cluster_column = 'leiden_clus',#'cell_types', 
                                   random_iter = 500,
                                   gene_set_1 = select_ligands,
                                   gene_set_2 = select_receptors, 
                                   verbose = FALSE) #no no, nyt toimi

## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(A1,
                                     spatial_network_name = 'Delaunay_network',
                                     cluster_column = 'leiden_clus', 
                                     random_iter = 500,
                                     gene_set_1 = select_ligands,
                                     gene_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'none')

## select top LR ##
selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)
top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]
top_LR_cell_ints <-  top_LR_cell_ints[1:25]

# anyNaN(A1)
# [1] FALSE
# > sum(is.na(A1@raw_exprs))# [1] 0
# > sum(is.na(A1@norm_expr))# [1] 0
# > sum(is.na(A1@norm_scaled_expr))# [1] 0...

plotCCcomDotplot(gobject = A1,
                 comScores = spatial_all_scores, #spatialsocres pitäis olla nan vapaa: sum(is.nan(DataFrame(spatial_all_scores)[,18])) or sum(is.na(spatial_all_scores)), and: sum(spatial_all_scores == '')
                 selected_LR = top_LR_ints, #no nans
                 selected_cell_LR = top_LR_cell_ints, #no nans
                 cluster_on = 'PI',aggl_method='median',
                 save_param = c(save_name = '14_a_communication_dotplot', save_format = 'pdf'))

# Footer
# © 2023 GitHub, Inc.
# human-oral-mucosa-spatial/OM-spatial.R at main · anacaetano/human-oral-mucosa-spatial · GitHub


#Deconvolution, https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html

#%%Solving the number of cells in the spots
#https://www.r-orms.org/mixed-integer-linear-programming/practicals/problem-tsp/
# https://towardsdatascience.com/linear-programming-in-r-444e9c199280
#Starting with epit spots:

LinkedDimPlot(data1a);

data=data1a
ok='D:/Data for Seurat Analysis/filtered_gene_bc_matrices/Data Filtered and Processed in Seurat/Spatial'; setwd(dir=ok)# 
i2=row.names(data@images[["slice1"]]@coordinates) #which(ii==grep('CACGTG', ii, value=TRUE)) check this
#fig 2. tot=c('GTCCGAGA','GGACAAGTT','ATGGCAGCC','ATACGGTGAA','ATCGGTTACC') #needs to have 6-7 to get unique
tot=c('CTGGGTAGG','ACTCTCT','GTTTGGCCC')
jep=c(); for(i in 1:length(tot)) {jep=append(jep,which(i2==grep(tot[i], i2, value=TRUE)))}
data_es=data[,jep]

#M spots
#fig. 2: tot=c('ACATCAGCT','AACTTTACGG','CACGGGATT','ATGTTCGTC','AGTGGTTGCG') #needs to have 6-7 to get unique
tot=c('CGCCTCCC','GTGAGGACA','CTGCGTG') #needs to have 6-7 to get unique

jep=c(); for(i in 1:length(tot)) {jep=append(jep,which(i2==grep(tot[i], i2, value=TRUE)))}
data_ms=data[,jep]

#Border (M and E) spot:
tot=c('GAGATCTTCC','CCCAGTTAAGG')
jep=c(); for(i in 1:length(tot)) {jep=append(jep,which(i2==grep(tot[i], i2, value=TRUE)))}
data_brd=data[,jep]

set.seed(102); 
test=data_es
DefaultAssay(test)='Spatial'
sce_A1n=as.SingleCellExperiment(test)
colnames(sce_A1n)=paste("Spot",1:length(colnames(sce_A1n)),sep="_")
reducedDims(sce_A1n)$PCA=reducedDims(sce_A1n)$PCA[1:length(colnames(sce_A1n)),]
sce_A1n@colData$row=-test@images[["slice1"]]@coordinates[["row"]]
sce_A1n@colData$col=-test@images[["slice1"]]@coordinates[["col"]]
sce_A1n@colData$imagerow=-test@images[["slice1"]]@coordinates[["imagerow"]]
sce_A1n@colData$imagecol=-test@images[["slice1"]]@coordinates[["imagecol"]]

# mes=sce_A1n
# brd=sce_A1n
# epit=sce_A1n

sce_A1n=brd
exprs_mat = SummarizedExperiment::assay(sce_A1n, 'logcounts'); set.seed(19)
result2 = scSEGIndex(exprs_mat = exprs_mat);#lower the better, let's take 200 and compare that to tgg
result2
df <- data.frame(x = seq(1:dim(result2)[1]),y = result2[with(result2,order(-segIdx)),][1:dim(result2)[1],'segIdx'])
plot(df)
ae=find_curve_elbow(df, plot_curve = TRUE);ae #:) [1] 130
comp2=result2[with(result2,order(-segIdx)),][1:150,]
cp2b=rownames(comp2)

# total=rowSums(sce_A1n[cp2e,]@assays@data[['counts']])
# total=unname(total)

#https://www.r-orms.org/mixed-integer-linear-programming/practicals/problem-course-assignment/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763026/
# https://www.r-orms.org/mixed-integer-linear-programming/packages/modelling-milp/
sce_A1n=brd
sgen=c(cp2b)#,cp2,cp2b
spot=2
# t2=sce_A1n@assays@data[['counts']][order(-sce_A1n[,]@assays@data[['counts']]),][sce_A1n@assays@data[['counts']][order(-sce_A1n[,]@assays@data[['counts']]),]>5]
# hist(sce_A1n@assays@data[['counts']][order(-sce_A1n[,]@assays@data[['counts']]),][sce_A1n@assays@data[['counts']][order(-sce_A1n[,]@assays@data[['counts']]),]>2],breaks=100,xlim=c(0,20))
total=rowSums(sce_A1n[sgen,]@assays@data[['counts']])
# total=t2
b=unname(rowSums(sce_A1n[sgen,spot]@assays@data[['counts']])) #rowSums(sce_A1n[sgen,1]@assays@data[['counts']]) sum(b>0)
ok=as.vector(sce_A1n[sgen,spot]@assays@data[['counts']]>0)# b=unname(rowSums(sce_A1n[sgen,1]@assays@data[['counts']]))
b=sce_A1n[sgen,spot]@assays@data[['counts']][ok,1]# b=t2
b1=unname(b)
n <- 30
m=length(names(b))
v <- unname(matrix(1, n, 1)[,1])
avgf=dim(sce_A1n)[2] #the number of spots in your evaluation
w=matrix(total[names(b)]/(avgf*n),m,n)

model <- MIPModel() %>%
  # 1 if cell i is assigned to spot
  add_variable(x[i], i = 1:n, type = "binary") %>%
  # maximize the cell preferences
  set_objective(sum_expr(v[i] * x[i], i = 1:n),"max") %>%
  # we cannot exceed the capacity 
  add_constraint(sum_expr(w[j, i] * x[i], i = 1:n) <= b1[j], j = 1:m) #%>% 
model
result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
sum(result$solution)

#border 1: 12, 2:15

#first epithelial spot: 10 cells, 12 for next cells.. 21, 11, 10
#Mesenchyme: 15, 15, 16, 15, 30

#first dataset image:
#epit: 15,,11
x1=c(12,15)
#Number of mes cells in the 'border': 3, and epit: 30 (?)
#Tää formulaatio olis lähinnä: https://web.mit.edu/15.053/www/AMP-Chapter-09.pdf
#or:https://ds-pl-r-book.netlify.app/optimization-in-r.html
# mes=sce_A1n
# brd=sce_A1n
# epit=sce_A1n

xe=colSums(epit[cp2e,]@assays@data[['counts']]) #
#first epithelial spot: 10 cells, 12 for next cells.. 21, 11, 10
# ne=c(10,12,21,11,10)

ne =c(15,11,11)

# Emulti=median(ceiling(xe/ne))
xm=colSums(mes[cp2,]@assays@data[['counts']])
# nm=c(15,15,16,15,30)
nm=c(12,15,9)
# Mmulti=median(ceiling(xm/nm))
#Mesenchyme: 15, 15, 16, 15, 30
xbe=colSums(brd[cp2e,]@assays@data[['counts']]) #61
xbm=colSums(brd[cp2,]@assays@data[['counts']]) #101

e=round(mean(xe),0);ee=round(mean(ne),0);xbee=unname(plyr::round_any(xbe, 5, f = floor)) #sum(c(71,85,80))/3;ee=sum(c(10,11,10))/3 #https://datacornering.com/round-roundup-rounddown-trunc-in-r/
ecf=round(xbee/e*ee,0)
m=round(mean(xm),0);mm=round(mean(nm),0);xbmm=unname(plyr::round_any(xbm, 5, f = floor)) #sum(c(29,35,47))/3 sum(c(15,16,15))/3
mcf=round(xbmm/m*mm,0)
Ncepit=round(ecf/(ecf+mcf)*x1,0)
Ncmes=round(mcf/(ecf+mcf)*x1,0) #seems to be working !

#8 cells are epit and the rest (22) are mesenchymal like
# xe/xm
# ne/nm# xbe/xbm*0.33333



