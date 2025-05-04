setwd("./")

options(stringsAsFactors = FALSE)
#### Seurat version 4.2 include sctransform v2
#### use sctransform v2 pipeline
library(future)
library(doFuture)
library(doParallel)
library(doRNG)
registerDoFuture()
plan("multisession", workers = 10)
#plan("sequential")
options(future.globals.maxSize = 400*1024^3, future.seed=TRUE)
set.seed(42)
#library("Seurat")
library(ggplot2)
library(patchwork)
library(dplyr)
#library("hdf5r")
library(doMC)
library(openxlsx)
library(gridExtra)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(VennDiagram)
library(clustree)
library(NMF)
library(harmony)
library(reshape2)
source("../../scripts/helper.r")
old_path <- Sys.getenv("PATH")

####################################################################################################################################
### filter cells
###### QC & loading data
datapaths <- c(list.files("./LBD-1", full.names = TRUE),
               list.files("./LBD-3", full.names = TRUE),
               list.files("./LBD-4", full.names = TRUE),
               list.files("./LBD-5", full.names = TRUE),
               list.files("./LBD-6", full.names = TRUE),
               list.files("./LBD-7", full.names = TRUE),
               list.files("./LBD-8", full.names = TRUE),
               list.files("./LBD-9", full.names = TRUE),
               list.files("./LBD-11", full.names = TRUE),
               list.files("./LBD-12", full.names = TRUE))
#set up sample names from paths
samplenames <- sapply(strsplit(gsub("/outs", "", basename(dirname(datapaths))), "-"), function(i) if(length(i) ==2) paste0(i[1],'_',i[2]) else i[1])

#### load sampleinfo file 
sampleinfo <- readxl::read_xlsx("Zhao_SpatialTranscriptomics_metadata.xlsx")

###### correct sample IDs
sampleinfo$datapath <- datapaths[match(sampleinfo$`Sample ID`,samplenames)]
sampleinfo <- sampleinfo[sampleinfo$`Sample ID` %in% samplenames,]
colnames(sampleinfo) <- c("sampleid", "patientid", "age", "sex", "tissue", "anatomy_region", "disease_state", "tissue_type", "panel", "slide_serial","capture_area", "image_id","apoe_geno","snrnaseq","datapath")
data("cc.genes")
###### do not use targeted data for now
#sampleinfo_sel <- sampleinfo[!sampleinfo$sampleid %in% c("KH3a", "KH3b"),]
sampleinfo_sel <- sampleinfo

#### change to v5 if you want to work with Seurat5
options(Seurat.object.assay.version = "v3")
dataall <- foreach(i = 1:nrow(sampleinfo_sel)) %dorng% {
  print(paste0("loading data ", i))
  tmp <- Load10X_Spatial_fixed(data.dir = sampleinfo_sel$datapath[i])
  metadata <- as.data.frame(sampleinfo_sel[rep(i,ncol(tmp)),1:14])
  rownames(metadata) <- colnames(tmp)
  tmp <- AddMetaData(tmp, metadata)
  tmp <- filtercells(tmp, nCount_cutoff = 500, saveqc = TRUE, qcoutput = sampleinfo_sel$sampleid[i])
  message(paste0('predict cell cycle phase for sample ', sampleinfo_sel$sampleid[i]))
  ### use normalize data before cellcyclescoring
  Idents(tmp) <- "orig.ident"
  tmp <- NormalizeData(tmp,assay = "Spatial")
  tmp <- CellCycleScoring(tmp, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  tmp <- SCTransform(tmp, assay = "Spatial", vst.flavor = "v2", vars.to.regress = c("percent.mt","Phase"))
  #name <- sampleinfo_sel$sampleid[i]
  tmp
}

####################################################################################################################################
#### try integration all samples
object.features <- SelectIntegrationFeatures(object.list = dataall, nfeatures = 3000)
object.features <- setdiff(object.features,unique(unlist(cc.genes)))
message("Generating merged data for harmony")
data2 <- merge(dataall[[1]], y = c(dataall[[2]],dataall[[3]],dataall[[4]],dataall[[5]],dataall[[6]],dataall[[7]],dataall[[8]],
                                   dataall[[9]],dataall[[10]]),  project = "seurat_harmony", merge.data = TRUE)
names(data2@images) <- sampleinfo_sel$sampleid
VariableFeatures(data2) <- object.features
table(Idents(data2))

saveRDS(data2, "data2.merged.rds") 
#data2 <- readRDS("data2.merged.rds")

rm(dataall)
rm(sampleinfo)
rm(sampleinfo_sel)

data2 <- RunPCA(object = data2, assay = "SCT", verbose = FALSE, npcs = 30)
message('Running harmony')
data2 <- RunHarmony(
  object = data2,
  assay.use = "SCT",
  reduction = "pca",
  dims.use = 1:30,
  group.by.vars = "slide_serial",
  plot_convergence = TRUE, kmeans_init_nstart=20, kmeans_init_iter_max=100
  #max.iter.cluster = 100, 
  #max.iter.harmony = 25,
)

data2 <- RunTSNE(data2, assay="SCT", reduction = "harmony", dims = 1:30)
data2 <- RunUMAP(object = data2, assay = "SCT", reduction = "harmony", dims = 1:30)
data2 <- RunPCA(object = data2, assay = "SCT", dims = 1:30)

saveRDS(data2, "integrated_nCount500_v3.rds")

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))


data2 <- FindNeighbors(object = data2, assay = "SCT", reduction = "harmony", dims = 1:30)
tmp <- FindClusters(data2, resolution = seq(0.1,1.5,0.1))
#clustree(tmp, prefix = "SCT_snn_res.")
p <- clustree(tmp, prefix = "SCT_snn_res.")
ggsave("./qc/LBD_integration_clustree_081823.png", p, width = 12, height = 12, dpi = 300, units = "in")
### use 0.4 resolution 
data2 <- FindClusters(object = data2, resolution = 0.2, cluster.name = "SCT_LowRes_0.2")
data2$SCT_LowRes_0.2 <- factor(data2$SCT_LowRes_0.2, levels = sort(as.numeric(levels(data2$SCT_LowRes_0.2))))

data2 <- FindClusters(object = data2, resolution = 0.4, cluster.name = "SCT_OptRes_0.4")
data2$SCT_OptRes_0.4 <- factor(data2$SCT_OptRes_0.4, levels = sort(as.numeric(levels(data2$SCT_OptRes_0.4))))

data2 <- FindClusters(object = data2, resolution = 0.8, cluster.name = "SCT_HighRes_0.8")
data2$SCT_HighRes_0.8 <- factor(data2$SCT_HighRes_0.8, levels = sort(as.numeric(levels(data2$SCT_HighRes_0.8))))

data2 <- FindClusters(object = data2, resolution = 0.6, cluster.name = "SCT_LayerRes_0.6")
data2$SCT_LayerRes_0.6 <- factor(data2$SCT_LayerRes_0.6, levels = sort(as.numeric(levels(data2$SCT_LayerRes_0.6))))


saveRDS(object = data2, file = "tmp.data2.rds")


# basic dimplot
p <- DimPlot(data2, reduction = "umap", group.by = "SCT_HighRes_0.8", label = F, combine = T, ncol = 1, raster = F, shuffle = T)
ggsave("./figures/LBD_integration_HighRes_umap.png", p, width = 12, height = 6*2/2, dpi = 300, units = "in")
p <- DimPlot(data2, reduction = "umap", group.by = "SCT_LowRes_0.2", label = F, combine = T, ncol = 1, raster = F, shuffle = T)
ggsave("./figures/LBD_integration_LowRes_umap.png", p, width = 12, height = 6*2/2, dpi = 300, units = "in")
p <- DimPlot(data2, reduction = "umap", group.by = "SCT_OptRes_0.4", label = F, combine = T, ncol = 1, raster = F, shuffle = T)
ggsave("./figures/LBD_integration_OptRes_umap.png", p, width = 12, height = 6*2/2, dpi = 300, units = "in")

Idents(data2) <- "SCT_OptRes_0.4"
p <- DimPlot(data2, reduction = "umap", split.by  = c("disease_state"), label = TRUE, combine = T, ncol = 2, raster = F, shuffle = T)
ggsave("./figures/LBD_integration_disease_umap.png", p, width = 12, height = 6*2/2, dpi = 300, units = "in")
### group by sampleid
p <- DimPlot(data2, reduction = "umap", group.by = "sampleid", label = F, combine = T, ncol = 1, raster = F, shuffle = T)
ggsave("./figures/LBD_integration_sampleid_umap.png", p, width = 12, height = 6*2/2, dpi = 300, units = "in")

### get scalesizes (function in helper script)
imagescales <- getImagePointSizes(data2)

### set up plotting to adjust for point size again probably a function
dimColors <- colorRampPalette(colors = brewer.pal(n = 9, name = "Blues"))

splots <- c()
aplots <- c()
#gplots <- c()
for (image in names(imagescales)) {
p <- SpatialFeaturePlot(data2, features = "SNCA", images = image, pt.size.factor = imagescales[[image]]*2,
                        alpha = c(0.1,1), image.alpha = 0.3)
p1 <- SpatialFeaturePlot(data2, features = "APOE", images = image, pt.size.factor = imagescales[[image]]*2,
                         alpha = c(0.1,1), image.alpha = 0.3)
#p2 <- SpatialFeaturePlot(data2, features = "GBA", images = image, pt.size.factor = imagescales[[image]]*2,
#                        alpha = c(0.1,1), image.alpha = 0.3)
splots[[image]] <- p + labs(title = image) + CenterTitle()
aplots[[image]] <- p1 + labs(title = image) + CenterTitle()
#gplots[[image]] <- p2 + labs(title = image) + CenterTitle()
}
### remove all objects that start with p and either a single character or nothing
rm(list = ls(pattern = "^p.*$"))

### create the final images
splot<- wrap_plots(splots, guides = 'collect', ncol = 2)
aplot <- wrap_plots(aplots, guides = 'collect', ncol = 2)
#gplot <- wrap_plots(gplots, guides = 'collect', ncol = 2)

## save the images to file 
ggsave('./figures/LBD_SNCA_spatialplot.png', splot, width = 8, height = 15, dpi = 600, units = "in")
ggsave('./figures/LBD_APOE_spatialplot.png', aplot, width = 12, height = 15, dpi = 600, units = "in")
ggsave('./figures/LBD_GBA_spatialplot.png', gplot, width = 12, height = 15, dpi = 600, units = "in")

### remove all objects that have plot in the name
rm(list = ls(pattern = "*plot"))


#### saved at end of 8/18/23
#saveRDS(object = data2, file = "tmp.data2.rds")
data2 <- readRDS("tmp.data2.rds")

### clustering results do not look good for integration all samples
## Determine DE genes per cluster per sample
#workflowpath <- "/research/bsi/tools/pipelines/scrna_10x/tertiary_pipeline/seurat/v3/"
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
#make it easier to find the markers
data2 <- PrepSCTFindMarkers(data2)
saveRDS(object = data2, file = "data2.prepped.rds")


data2 <- readRDS("data2.prepped.rds")
data2 <- readRDS("bothCohorts_spatial_withAdjExpAssays_090123.rds")
#choose clusters to work with
selectedclusters <- as.numeric(which(!apply(table(data2$disease_state, data2$SCT_OptRes_0.4),2,function(i) any(i==0)))-1)
Idents(data2) <- "SCT_OptRes_0.4"
### dropping cluster 11 too few cells to analyze
#selectedclusters <- as.numeric(seq(0,10,1))
cat("selected clusters: ", selectedclusters)
#Start actually calling DEs
# cluster_list <- foreach(i = 0:(length(unique(data2$SCT_OptRes_0.4))-1),.inorder=TRUE, .errorhandling='stop') %dopar% {
#   FindMarkers(data2, ident.1 = i, 
#                        #grouping.var = "disease_state", 
#                        verbose = TRUE,min.cells.group=1,
#                        min.pct=0,min.cells.feature=0,logfc.threshold=0,assay = 'SCT',  max.cells.per.ident = 1000,
#                        slot = 'data') 
# }

### running as background job
cluster_list_all <- FindAllMarkers(data2,
                      #grouping.var = "disease_state", 
                      verbose = TRUE,min.cells.group=1,
                      min.pct=0,min.cells.feature=0,logfc.threshold=0,assay = 'SCT', 
                      #max.cells.per.ident = 1000,
                      slot = 'data', ) 


saveRDS(cluster_list, "cluster_list_090123.rds")

#cluster_list <- readRDS("cluster_list_090123.rds")
cluster_list <- split(cluster_list, cluster_list$cluster)
names(cluster_list) <- paste0("cluster", names(cluster_list))

### and BH adjusted p-val
cluster_list <- lapply(1:length(unique(data2$SCT_OptRes_0.4)),function(i) {
  tmp <- cluster_list[[i]]
  if( (i-1) %in% selectedclusters) {
      #tmp[,2] <- tmp[,2]*log2(2) ### this doesn't do anything( just multiplies by 1...)
      colnames(tmp)[2] <- 'avg_log2FC'
      tmp[['p_val_adj_BH']] <- p.adjust(tmp[,1],method = 'BH')
      tmp[order(tmp[,2],decreasing=T),]}
})
names(cluster_list) <- paste0("cluster",unlist(lapply(1:length(cluster_list), function(x) { unlist(cluster_list[[x]][1,"cluster"])})))


cluster_list_filter <- lapply(1:length(unique(data2$SCT_OptRes_0.4)), function(i) {
  tmp <- cluster_list[[i]]
  if( (i-1) %in% selectedclusters)
        tmp[(abs(tmp[,2])>log2(1.2) & tmp[,8]<0.05),]
})

conserved_cluster_list <- readRDS("conserved_cluster_list2023-09-01.rds")
conserved_cluster_list <- lapply(names(conserved_cluster_list), function(x){
  conserved_cluster_list[[x]]$cluster <- gsub(pattern = "cluster",replacement = "",x)
  conserved_cluster_list[[x]]
})

conserved_cluster_list <- lapply(1:length(unique(data2$SCT_OptRes_0.4)),function(i) {
  tmp <- conserved_cluster_list[[i]]
  if( (i-1) %in% selectedclusters) {
    tmp[,2] <- tmp[,2]*log2(2)
    tmp[,7] <- tmp[,7]*log2(2)
    tmp[,12] <- tmp[,12]*log2(2)
    #colnames(tmp)[2] <- gsub('logFC','log2FC',colnames(tmp)[2])
    #colnames(tmp)[7] <- gsub('logFC','log2FC',colnames(tmp)[7])
    #colnames(tmp)[12] <- gsub('logFC','log2FC',colnames(tmp)[12])
    tmp[[gsub('_p_val','_p_val_adj_BH',colnames(tmp)[1])]] <- p.adjust(tmp[,1],method = 'BH')
    tmp[[gsub('_p_val','_p_val_adj_BH',colnames(tmp)[6])]] <- p.adjust(tmp[,6],method = 'BH')
    tmp[[gsub('_p_val','_p_val_adj_BH',colnames(tmp)[11])]] <- p.adjust(tmp[,11],method = 'BH')
    tmp <- tmp[,c(1:5,19,6:10,20,11:15,21,16,17,18)]
    tmp <- tmp[which(tmp[,2]*tmp[,8]*tmp[,14] >0),]
    tmp[order(tmp[,2],decreasing=T),]} else{
      tmp[,2] <- tmp[,2]*log2(2)
      colnames(tmp)[2] <- 'avg_log2FC'
      tmp[['p_val_adj_BH']] <- p.adjust(tmp[,1],method = 'BH')
      tmp[order(tmp[,2],decreasing=T),]}
})
names(conserved_cluster_list) <- paste0("cluster",unlist(lapply(1:length(conserved_cluster_list), function(x) { unlist(conserved_cluster_list[[x]][1,"cluster"])})))


conserved_cluster_list_filter <- lapply(1:length(unique(data2$SCT_OptRes_0.4)), function(i) {
  tmp <- conserved_cluster_list[[i]]
  if( (i-1) %in% selectedclusters)
    ## only select positive conserved genes
    tmp[(tmp[,2]>log2(1.2) & tmp[,6]<0.05) &
          (tmp[,8]>log2(1.2) & tmp[,12]<0.05) &
          (tmp[,14]>log2(1.2) & tmp[,18]<0.05),] else
      tmp[(tmp[,2]>log2(1.2) & tmp[,6]<0.05),]
})
names(conserved_cluster_list_filter) <- paste0("cluster",unlist(lapply(1:length(conserved_cluster_list_filter), function(x) { unlist(conserved_cluster_list_filter[[x]][1,"cluster"])})))

conserved_cluster_list_filter <- lapply(conserved_cluster_list_filter, function(i) {
  i$gene <- rownames(i)
  d1 <- cbind(gene_symbol=i$gene,human_anno[match(i$gene,human_anno$hgnc_symbol),-1],i,stringsAsFactors=F) 
  rownames(d1) <- d1[,1]
  d1 <- d1[,colnames(d1) != "gene"]
  d1
})

saveRDS(conserved_cluster_list_filter,"conserved_cluster_list_filtered_2023-09-01.rds")

###add annotations to full list and filtered
human_anno <- read.delim('./seurat/v4.1/docs/human_anno.txt',
                         header=T,stringsAsFactors=F,check.names=F)
cluster_list <- lapply(cluster_list, function(i) {
  d1 <- cbind(gene_symbol=i$gene,human_anno[match(i$gene,human_anno$hgnc_symbol),-1],i,stringsAsFactors=F)
  rownames(d1) <- d1[,1]
  d1 <- d1[,colnames(d1) != "gene"]
  d1
})
names(cluster_list) <- paste0("cluster",unlist(lapply(1:length(cluster_list), function(x) { unlist(cluster_list[[x]][1,"cluster"])})))

cluster_list_filter <- lapply(cluster_list_filter, function(i) {
  d1 <- cbind(gene_symbol=i$gene,human_anno[match(i$gene,human_anno$hgnc_symbol),-1],i,stringsAsFactors=F) 
  rownames(d1) <- d1[,1]
  d1 <- d1[,colnames(d1) != "gene"]
  d1
})
names(cluster_list_filter) <- paste0("cluster",unlist(lapply(1:length(cluster_list_filter), function(x) { unlist(cluster_list_filter[[x]][1,"cluster"])})))

saveRDS(cluster_list_filter, "cluster_list_filter_090123.rds")
saveRDS(cluster_list, "cluster_list_annotated_090123.rds")

cluster_list_filter_up <- lapply(names(cluster_list_filter), function(i) {
  tmp <- cluster_list_filter[[i]]
  tmp[tmp[,"avg_log2FC"]>0,]
})
names(cluster_list_filter_up) <- paste0("cluster",unlist(lapply(1:length(cluster_list_filter_up), function(x) { unlist(cluster_list_filter_up[[x]][1,"cluster"])})))

saveRDS(cluster_list_filter_up,"cluster_list_filter_PosOnly_090123.rds")

### write data to excel sheets
output <- paste0(getwd(),"/DEGs")
dir.create(output)
prefix <- "LBD"
wb <- createWorkbook()
for(i in 0:(length(unique(data2$SCT_OptRes_0.4))-1)) {
  addWorksheet(wb, paste0('cluster',i))
  writeData(wb, sheet = i+1,cluster_list[[i+1]],rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, paste0(output,'/',prefix,'_marker_genes.xlsx'), overwrite = T)

wb <- createWorkbook()
for(i in 0:(length(unique(data2$SCT_OptRes_0.4))-1)) {
  addWorksheet(wb, paste0('cluster',i))
  writeData(wb, sheet = i+1,cluster_list_filter[[i+1]],rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, paste0(output,'/',prefix,'_marker_genes_filter.xlsx'), overwrite = T)

wb <- createWorkbook()
for(i in 0:(length(unique(data2$SCT_OptRes_0.4))-1)) {
  addWorksheet(wb, paste0('cluster',i))
  writeData(wb, sheet = i+1,cluster_list_filter_up[[i+1]],rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, paste0(output,'/',prefix,'_marker_genes_filter_PosOnly.xlsx'), overwrite = T)

saveRDS(cluster_list,"annotated_conservedClusterList.rds")

#### fix factors/levels
data2$sampleid <- factor(data2$sampleid, levels = paste0("LBD_", c(1,3:9,11,12)))
#data2$sampleid <- factor(data2$sampleid, levels = paste0("LBD_", c(1,3:9,11,12)))
data2$SCT_OptRes_0.4 <- factor(data2$SCT_OptRes_0.4, levels=sort(as.numeric(levels(data2$SCT_OptRes_0.4))))

data2$seurat_clusters <- data2$SCT_OptRes_0.4

saveRDS(data2, "bothCohorts_spatial_withAdjExpAssays_090123.rds")

### select out genes specific to clusters
genelist_conserv <- sapply(0:(length(unique(data2$seurat_clusters))-1), function(i) {rownames(cluster_list[[i+1]])[1]})
labels <- paste0('Cluster', 0:(length(genelist_conserv)-1),':',genelist_conserv)
labels <- labels[!is.na(genelist_conserv)]
labels <- rep(labels,each=2)
genelist_conserv <- genelist_conserv[!is.na(genelist_conserv)]
p <- lapply(genelist_conserv,function(i) {
  limits <- c(0,max(data2@assays$SCT@data[rownames(data2@assays$SCT@data) == i,]))
  p1<- SpatialFeaturePlot(data2, features = i, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.2,0.8),
                          images = imageNo3, 
                          #stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = i, 
                           pt.size.factor = 1, 
                           #alpha = c(0.2,0.8),
                           images = "LBD_3", 
                           #stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1, p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
})
p <- unlist(p,recursive = F)
p <- lapply(1:length(p),function(i) p[[i]] + labs(title = labels[i]) + theme(plot.title = element_text(hjust = 0.5)))
m <- grid.arrange(grobs=p,ncol=length(data2@images),nrow=length(p)/length(data2@images))
ggsave(paste0(output,'/',prefix,'_conserve_gene.png'),m,width=6*length(data2@images),height=6*length(p)/length(data2@images),limitsize = F,dpi=300)

### barchart of proportion of spots in each cluster
temp <- melt(table(data2$SCT_OptRes_0.4,data2$sampleid))
temp <- temp %>% group_by(Var2) %>% mutate(pct=value/sum(value))
gp <- ggplot() + geom_col(aes(x=Var1, y=pct,fill=Var2),position = position_dodge(), data=temp) +
  theme_classic()+ scale_x_continuous(breaks=0:(length(unique(temp[[1]]))-1)) +
  theme(legend.title = element_blank(), legend.text = element_text(size=15), axis.text = element_text(size=15), axis.title = element_text(size=18,face="bold")) +
  labs(x='cluster number', y='percentage of spots')
ggsave('./figures/LBD_pctofspotsbarplot.png',gp,width=12,height=12)



#find spatially variable features
assayList <- c("Spatial","SCT","adjExp_Ast","adjExp_Endo","adjExp_Ex","adjExp_In","adjExp_Mic","adjExp_OPC","adjExp_Olig","adjExp_Peri")
top.features <- c()
foreach(assay %in% assayList) %dorng% {
  SeuratObj <- FindSpatiallyVariableFeatures(SeuratObj,
                                       assay = assay, selection.method = "markvariogram")
  top.features[[assay]] <- SpatiallyVariableFeatures(data2, selection.method = "markvariogram", assay = assay)
}
saveRDS(SeuratObj, "spatial_withAdjExpAssays_wSpatVarFeat_071223.rds")
saveRDS(top.features, "spatiallyVariable_features.rds")

output <- "./figures"
prefix <- "LBD"
topSpatial <- head(top.features,6)
labels <- paste0('Cluster', 0:(length(topSpatial)-1),':',topSpatial)
labels <- labels[!is.na(topSpatial)]
labels <- rep(labels,each=2)
p <- lapply(topSpatial,function(i) {
  limits <- c(0,max(test@assays$SCT@data[rownames(test@assays$SCT@data) == i,]))
  p1 <- SpatialFeaturePlot(test, features = i,combine = F)
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
})
p <- unlist(p,recursive = F)
p <- lapply(1:length(p),function(i) p[[i]] + labs(title = labels[i]))
m <- grid.arrange(grobs=p,ncol=length(test@images),nrow=length(p)/length(test@images))
ggsave(paste0(output,'/',prefix,'_spatialFeatures.png'),m,width=6*length(test@images),height=6*length(p)/length(test@images),limitsize = F,dpi=300)

top.features <- human_anno[match(top.features,human_anno$hgnc_symbol),]
wb <- createWorkbook()
addWorksheet(wb, "spatialFeatures")
writeData(wb, sheet = "spatialFeatures",top.features,rowNames=FALSE,colNames=TRUE)
saveWorkbook(wb, "./LBD_spatialFeatures.xlsx", overwrite = TRUE)


####################################################################################################################################
############################### determine the best strategy for integration
######## merge data for Ascending region
#dataall <- readRDS("dataall.rds")
#sampleinfo <- read.delim("sampleinfo.txt", header = T)
#sampleinfo_sel <- sampleinfo[!sampleinfo$sampleid %in% c("KH3a", "KH3b"),]
index <- which(sampleinfo_sel$anatomy_region == "Frontal Cortex")
data3 <- merge(dataall[[index[1]]], dataall[index[2:length(index)]])
names(data3@images) <- sampleinfo_sel$sampleid[index]

DefaultAssay(data3) <- "SCT"
VariableFeatures(data3) <- setdiff(unique(c(sapply(index, function(i) VariableFeatures(dataall[[i]])))),unique(unlist(cc.genes)))
#### Identify optimal # of clusters
data3 <- RunPCA(data3, verbose = FALSE)
data3 <- FindNeighbors(data3, dims = 1:30)
tmp <- FindClusters(data3, resolution = seq(0.1,1.5,0.2))
clustree(tmp, prefix = "SCT_snn_res.")
p <- clustree(tmp, prefix = "SCT_snn_res.")
ggsave("./figures/LBD_merge_clustree.png", p, width = 12, height = 12, dpi = 300, units = "in")
#### choose value 0.4 as resolution number
data3 <- FindClusters(data3, verbose = FALSE,resolution=0.4)
data3 <- RunUMAP(data3, dims = 1:30)
p <- DimPlot(data3, reduction = "umap", split.by  = c("sampleid"), label = TRUE)
ggsave("./figures/LBD_merge_umap.png", p, width = 12, height = 6, dpi = 300, units = "in")
p <- SpatialDimPlot(data3, label = TRUE)
ggsave('./figures/LBD_merge_spatialplot.png', p, width = 12, height = 6, dpi = 300, units = "in")

#############################################
### requested by Zonghua
#############################################
analysis_genes <- readxl::read_xlsx("../metadata/CellMarkers.xlsx")
colnames(analysis_genes) <- gsub(" ", "_", colnames(analysis_genes))
colnames(analysis_genes) <- gsub("/", "_", colnames(analysis_genes))
output="../processing/figures"

### basic visualization of each gene
for( column in 1:length(analysis_genes)){ 
  prefix=colnames(analysis_genes)[column]
  print(prefix)
  genelist <- analysis_genes[,column][!is.na(analysis_genes[,column])]
  droplist <- NULL
  for(idx in 1:length(genelist)){
    tmp <- NA
    try(tmp <- FetchData(data2,genelist[idx]))
    if(!is.data.frame(tmp)){
      droplist <- c(droplist,idx)
    }
  }
  rm(tmp)
  print(paste(prefix,"dropping: "))
  print(droplist)
  if (!is.null(droplist)){
    genelist <- genelist[-droplist]  
  }
  labels <- paste0(prefix,':',genelist)
  labels <- rep(labels,each=6)
  p <- lapply(genelist,function(i) {
    limits <- c(0,max(data2@assays$SCT@data[rownames(data2@assays$SCT@data) == i,]))
    p1 <- SpatialFeaturePlot(data2, features = i,combine = F, alpha = c(0.1,.99))
    p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  })
  p <- unlist(p,recursive = F)
  p <- lapply(1:length(p),function(i) p[[i]] + labs(title = labels[i]) + theme(plot.title = element_text(hjust = 0.5)))
  m <- grid.arrange(grobs=p,ncol=length(data2@images),nrow=length(p)/length(data2@images))
  ggsave(paste0(output,'/',prefix,'_analysisFeatures.png'),m,width=6*length(data2@images),height=6*length(p)/length(data2@images),limitsize = F,dpi=300)
}

### add annotation modules
prefix="./figures/"
for(column in 1:length(analysis_genes)) {
  genelist <- analysis_genes[,column][!is.na(analysis_genes[,column])]
  annoName <- colnames(analysis_genes)[column]
  scoreName <- paste0(annoName,"1")
  data2 <- AddModuleScore(data2, features =list(genelist), assay = "SCT", name = annoName)
  limits <- c(min(data2[[scoreName]]), max(data2[[scoreName]]))
  p<- SpatialFeaturePlot(data2, features = scoreName, 
                         pt.size.factor = 3, alpha = c(0.1,.99)) &
    scale_fill_gradientn(limits=limits,colors = SpatialColors(n = 100), name=annoName)
  ggsave(paste0(prefix,annoName,"_AnnotationPredictions.png"), p, limitsize = F, dpi=300)
  # p <- lapply(scoreName,function(i) {
  #   limits <- c(min(test$scoreName),max(test$scoreName))
  #   p1 <- SpatialFeaturePlot(test, features = i,combine = F)
  #   p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  # })
  #p <- unlist(p,recursive = F)
  #p <- lapply(1:length(p),function(i) p[[i]] + labs(title = labels[i]))
  #m <- grid.arrange(grobs=p,ncol=length(data2@images),nrow=length(p)/length(data2@images))
}

### cluster DEGs
degData2 <- FindAllMarkers(data2, assay = 'SCT', slot = 'data', min.cells.group=1,
                           min.pct=0,min.cells.feature=0,logfc.threshold=0,
                           verbose = TRUE)
saveRDS(degData2, file = "degData2.rds")
wb <- createWorkbook()
for(clustNum in unique(degData2$cluster)){
  tmp <- degData2[degData2$cluster == clustNum,]
  addWorksheet(wb = wb, sheetName = paste0("degFeatures_cluster",clustNum))
  writeData(wb = wb, paste0("degFeatures_cluster",clustNum), tmp,rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, "./DEGs/LBD_degFeaturesPerCluster.xlsx", overwrite = TRUE)

wb <- createWorkbook()
for(clustNum in unique(degData2$cluster)){
  tmp <- degData2[((degData2$cluster == clustNum) & (degData2$p_val_adj <= 0.05) & 
                     (abs(degData2$avg_log2FC) >= 0.1)),]
  addWorksheet(wb = wb, sheetName = paste0("degFeatures_cluster",clustNum))
  writeData(wb = wb, paste0("degFeatures_cluster",clustNum), tmp,rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, "./DEGs/LBD_degFeaturesPerCluster_filtered.xlsx", overwrite = TRUE)

wb <- createWorkbook()
for(clustNum in unique(degData2$cluster)){
  tmp <- degData2[((degData2$cluster == clustNum) & (degData2$p_val_adj <= 0.05) & 
                     (abs(degData2$avg_log2FC) >= 0.1)),]
  
  tmp2 <- degData2[((degData2$cluster != clustNum) & (degData2$p_val_adj <= 0.05) & 
                      (abs(degData2$avg_log2FC) >= 0.1)),]
  tmp <- tmp[!tmp$gene %in% tmp2$gene, ]
  addWorksheet(wb = wb, sheetName = paste0("unique-degFeatures_cluster",clustNum))
  writeData(wb = wb, paste0("unique-degFeatures_cluster",clustNum), tmp,rowNames=FALSE,colNames=TRUE)
}
saveWorkbook(wb, "./DEGs/LBD_unique-degFeaturesPerCluster_filtered.xlsx", overwrite = TRUE)

### plots of DEG marker genes
genelist <- sapply(0:(length(unique(data2$seurat_clusters))-1), function(i) {filteredDEG_list[[i+1]]$gene[1]})
labels <- paste0('Cluster', 0:(length(genelist)-1),':',genelist)
labels <- labels[!is.na(genelist)]
labels <- rep(labels,each=6)

genelist <- genelist[!is.na(genelist)]
p <- lapply(genelist,function(i) {
  label <- grep(i, labels)
  out <- gsub(":.*", "", labels[[label]])
  limits <- c(0,max(data2@assays$SCT@data[rownames(data2@assays$SCT@data) == i,]))
  p1<- SpatialFeaturePlot(data2, features = i, 
                          pt.size.factor = 2.4, 
                          #alpha = c(0.2,0.9),
                          images = imageNo3, 
                          #stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = i, 
                           pt.size.factor = 1.5, 
                           #alpha = c(0.2,0.9),
                           images = "LBD_3", 
                           #stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1, p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = paste0(imageOrder[[i]],"_",out)) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- wrap_plots(p1,ncol=3, guides = 'collect')
  ggsave(paste0("./DEGs/LBD_",out,"_topclustermarker.png"), m, limitsize = F, dpi=300, width = 8, height = 4)
})
p <- unlist(p,recursive = F)
p <- lapply(1:length(p),function(i) p[[i]] + labs(title = labels[i]) + theme(plot.title = element_text(hjust = 0.5)))
m <- grid.arrange(grobs=p,ncol=length(data2@images),nrow=length(p)/length(data2@images))
output <- "./DEGs"
prefix <- "LBD"
ggsave(paste0(output,'/',prefix,'_topclustermarker_gene.png'),m,width=6*length(data2@images),height=6*length(p)/length(data2@images),limitsize = F,dpi=300)

##### nFeature_Spatial, nCount-Spatial, percent.mt , percent.rb, and percent.hb on umap
prefix="./qc/"
for(feature in c("nFeature_Spatial", "nCount_Spatial", "percent.mt" , "percent.rb", "percent.hb")){
  try({p <- FeaturePlot(data2, split.by = "sampleid", features = feature, 
              keep.scale = "feature" , order = TRUE, combine = F) #, min.cutoff = "q1", max.cutoff = "q99"))
  m <- wrap_plots(p, guides = 'collect') * ylab(label = NULL) * xlab(label = NULL) #* NoAxes(keep.text = F, keep.ticks = T)
  ggsave(paste0(prefix,feature,"_UMAP.png"), m, limitsize = F, dpi=300)})
}

### Layer vis:
layer23 <- "CUX2"
layer4 <- "RORB"
layer5 <- c("PCP4","HTR2C")
layer6 <- c("OPRK1","NR4A2")
simpleLayer <- list(layer23,layer4,layer5,layer6)
data2 <- AddModuleScore(data2, features = simpleLayer, name = "SimpleLayer")
names(data2@meta.data)[grep("Simple",names(data2@meta.data))] <- c("SimpleLayer23",
                                                                   "SimpleLayer4",
                                                                   "SimpleLayer5",
                                                                   "SimpleLayer6")

df <- data2@meta.data[c("SimpleLayer23",
                        "SimpleLayer4",
                        "SimpleLayer5",
                        "SimpleLayer6")]
df$SimplePrediction <- ifelse(apply(df, 1, max, na.rm=TRUE) < 0,
                        "Unknown",
                        names(df)[max.col(df)])

simpLayerCols <- c(#"Layers23" = "#5E4FA2",
               "SimpleLayer23" = "#5E4FA2",
               #"cLayers3" = "#BEE4A0",
               "SimpleLayer4" = "#88CFA4",
               "SimpleLayer5" = "#FFFFBF",
               "SimpleLayer6" = "#F88D52",
               "Unknown" = "#9E0142")

data2@meta.data$SimplePrediction <- df$SimplePrediction


p1<- SpatialDimPlot(data2, group.by = "SimplePrediction", 
                    cols = simpLayerCols,
                          pt.size.factor = 2.5, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
names(p1) <- imageNo3
p2 <- SpatialDimPlot(data2, group.by = "SimplePrediction", 
                     cols = simpLayerCols, 
                           pt.size.factor = 1.5625, 
                           #alpha = c(0.1,0.9),
                           images = "LBD_3", 
                           stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
names(p2) <- "LBD_3"
p1 <- c(p1,p2)
p1 <- p1[imageOrder]
#p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
m <- wrap_plots(p1, guides = 'collect', nrow = 1)
ggsave(paste0("./figures/SimplePrediction_spatialDistribution_051523.png"), m, limitsize = F, dpi=300, width = 14, height = 4)



genes <- c(layer23,layer4,layer5,layer6)

DefaultAssay(data2) <- "SCT"
lapply(genes, function(gene) {
  limits <- c(0,max(data2@assays$SCT@data[rownames(data2@assays$SCT@data) == gene,]))
  p1<- SpatialFeaturePlot(data2, features = gene, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = gene, 
                           pt.size.factor = 1, 
                           #alpha = c(0.1,0.9),
                           images = "LBD_3", 
                           stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1,p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- grid.arrange(grobs=p1,ncol=length(data2@images),nrow=length(p1)/length(data2@images))
  ggsave(paste0("./figures/",gene,"_spatialDistribution_051523.png"), m, limitsize = F, dpi=300, width = 12, height = 8)
})

controlGenes <- c("GAPDH", "HPRT1","B2M","ACTB","HMBS","PPIA")
gliosis <- c("AKT3","LRPPRC","ADAR","TFG","FARS2","LRRK2","ADH1C",
             "COX15","CSF1R","CYP27A1","ARX","DTYMK","DYRK1A","EIF2B1",
             "EIF4G1","C9orf72","ERCC6","ETFA","ETFB","ETFDH","BRAT1",
             "PMPCA","DNAJC13","TARDBP","NUP62","MTOR","FUS","FBXO7",
             "NDUFAF3","GIGYF2","GBA1","GLUD2","TSEN54","GRN","TBK1",
             "UBQLN2","HTT","HSD17B4","LMNB1","CHCHD10","MAN2B1","MAPT",
             "ATXN3","MOCS1","MOCS2","MT-ATP6","MT-TT","NDUFB8","NDUFS2",
             "NDUFS8","NR4A2","PRKN","PAX2","PDHA1","KDM3B","ATP6V1A",
             "SERPINI1","PIGA","PIK3CA","PLP1","POLG","GDAP2","AVP","VPS35",
             "NAXD","NGLY1","PRNP","PSEN1","THOC2","KCNT1","KMT2C","RANBP2",
             "BCS1L","ATXN2","ATXN8OS","ZNF335","SNCA","SURF1","CNTN2",
             "TBCD","TBP","TLR3","TYROBP","VCP","VRK1","YY1","NARS2","EHMT1",
             "L2HGDH","PLA2G6","SQSTM1","EIF2B4","EIF2B3","EIF2B2","EIF2B5",
             "SLC25A46","AP4M1","LONP1","SNCAIP","SCO2")
layerGenes <- read.csv("layerMarkerGenes_LIBD_HumanPilot.csv", header = T)
layerGeneList <- split(layerGenes, ~models)
layerGeneList <- lapply(layerGeneList, function(x) {
  x$gene_name
})
layerGeneList$Layer3
SpatialFeaturePlot(data2, features = layerGeneList$Layer3) & NoLegend()

DefaultAssay(data2) <- "SCT"
lapply(controlGenes, function(gene) {
  limits <- c(0,max(data2@assays$SCT@data[rownames(data2@assays$SCT@data) == gene,]))
  p1<- SpatialFeaturePlot(data2, features = gene, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = gene, 
                           pt.size.factor = 1, 
                           #alpha = c(0.1,0.9),
                           images = "LBD_3", 
                           stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1,p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- grid.arrange(grobs=p1,ncol=length(data2@images),nrow=length(p1)/length(data2@images))
  ggsave(paste0("./figures/",gene,"_ControlGenes_spatialDistribution_051523.png"), m, limitsize = F, dpi=300, width = 12, height = 8)
})

DefaultAssay(data2) <- "adjExp_Ex"
lapply(names(layerGeneList), function(layerNum) {
  limits <- range(data2@meta.data[[layerNum]])
  p1<- SpatialFeaturePlot(data2, features = layerNum, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = layerNum, 
                           pt.size.factor = 1, 
                           #alpha = c(0.1,0.9),
                           images = "LBD_3", 
                           stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1,p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- grid.arrange(grobs=p1,ncol=length(data2@images),nrow=length(p1)/length(data2@images))
  ggsave(paste0("./figures/Ex-expression_",layerNum,"_LayerPredictLIBD_spatialDistribution_051523.png"), m, limitsize = F, dpi=300, width = 12, height = 8)
  dev.off()
})


layers <- paste0("Layer",seq(1,6)) 
combinedLayers <- lapply(layers, function(lay){
  tmpList <- names(layerGeneList)[grep(lay, names(layerGeneList))]
  out <- c()
  for (tmp in tmpList) {
    out <- c(out,layerGeneList[[tmp]])
  }
  out
  }) 
names(combinedLayers) <- layers

data2 <- AddModuleScore(data2, features = combinedLayers, name="cLayers_Ex", assay = "adjExp_Ex", )
data2 <- AddModuleScore(data2, features = list(combinedLayers[[1]]), name="cLayers_Ast", assay = "adjExp_Ast", )


lapply(paste0("cLayers_Ex",seq(1,6)), function(layerNum) {
  limits <- c(0,max(data2@meta.data[[layerNum]]))
  p1<- SpatialFeaturePlot(data2, features = layerNum, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(data2, features = layerNum, 
                           pt.size.factor = 1, 
                           #alpha = c(0.1,0.9),
                           images = "LBD_3", 
                           stroke = 0.01,
                           #min.cutoff = "q1", 
                           #max.cutoff = "q99", 
                           combine=F
  )
  names(p2) <- "LBD_3"
  p1 <- c(p1,p2)
  p1 <- p1[imageOrder]
  p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- grid.arrange(grobs=p1,ncol=length(data2@images),nrow=length(p1)/length(data2@images))
  ggsave(paste0("./figures/Ast-expression",layerNum,"_combinedLayerPredictLIBD_spatialDistribution_052523.png"), m, limitsize = F, dpi=300, width = 12, height = 8)
  dev.off()
})

layerCols <- c("cLayers1" = "#5E4FA2",
               "cLayers2" = "#54AEAD",
               "cLayers3" = "#BEE4A0",
               "cLayers4" = "#FFFFBF",
               "cLayers5" = "#FDBE6F",
               "cLayers6" = "#E95D46",
               "Unknown" = "#9E0142")

head(max.col(data2@meta.data[c(paste0("cLayers",seq(1,6)))]))
df <- data2@meta.data[c(paste0("cLayers",seq(1,6)))]
df$prediction <- ifelse(apply(df, 1, max, na.rm=TRUE) < 0,
                                          "Unknown",
                                          names(df)[max.col(df)])
data2@meta.data$predictionLayer <- df$prediction

p1<- SpatialDimPlot(data2, group.by = "predictionLayer", 
                    cols = layerCols,
                    pt.size.factor = 2.5, 
                    #alpha = c(0.1,0.9),
                    images = imageNo3, 
                    stroke = 0.01,
                    #min.cutoff = "q1", 
                    #max.cutoff = "q99", 
                    combine=F
)
names(p1) <- imageNo3
p2 <- SpatialDimPlot(data2, group.by = "predictionLayer", 
                     cols = layerCols, 
                     pt.size.factor = 1.5625, 
                     #alpha = c(0.1,0.9),
                     images = "LBD_3", 
                     stroke = 0.01,
                     #min.cutoff = "q1", 
                     #max.cutoff = "q99", 
                     combine=F
)
names(p2) <- "LBD_3"
p1 <- c(p1,p2)
p1 <- p1[imageOrder]
#p1 <- lapply(p1,function(j) j + scale_fill_gradientn(limits=limits,colours = SpatialColors(n = 100)))
p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = labels[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
#p1 <- unlist(p1,recursive = F)
m <- wrap_plots(p1, guides = 'collect', nrow = 1)
ggsave(paste0("./figures/CombinedLayerPrediction_LIBD_spatialDistribution_051523.png"), m, limitsize = F, dpi=300, width = 14, height = 4)

