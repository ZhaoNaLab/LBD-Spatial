setwd("./2ndCohort/spacexr")

#### convert seurat object to anndata for neigborhood enrichment analysis
#### use R 4.2.2
options(stringsAsFactors = FALSE)
.libPaths(c("./R-4.2.2-latest/", "./Rlib4.2.2/"))
library("Seurat", lib.loc = "./R-4.2.2-latest//")
library('SeuratData')
library(SeuratDisk)
library(ggplot2)
library("hdf5r", lib.loc = "./R-4.2.2-latest/")
library(future)
library(doFuture)
registerDoFuture()
plan("multisession", workers = 10)
options(future.globals.maxSize = 400*1024^3)
#options(timeout = 600000000) ### set this to avoid timeout error
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE, lib="./Rlib4.2.2/")
library(spacexr)
library(Matrix)
library(doParallel)
#install.packages("celda")
#library(celda)
library(RColorBrewer)
library(gridExtra)
library(Matrix)
library(ggpubr)
library(patchwork)
library(dplyr)
source("../scripts/spacexr_functions.r")

### to make pipeline-able
# scRNAseq rds location: scrna = "./Seurat.YJ.final.rds"
# annotation column name: ann = "Subcluster"
# spatial assay name: spAssay = "Spatial"
# spatial rds location: spatial = "tmp.data2.rds"
# doublet mode: doublet = "full"
# output new seurat rds location: out = "deconvolutionOutput.rds"


### prep scRNAseq as reference for deconvolution
# need 3 things:
# 1. untransformed count-level data counts.csv
# 2. cell_types: A named (by cell barcode) factor of cell type for each cell. 
#     The ‘levels’ of the factor would be the possible cell type identities
# 3. optional nUMI:Optional, a named (by cell barcode) list of total counts or UMI’s appearing at each pixel. 
#    If not provided, nUMI will be assumed to be the total counts appearing on each pixel.
snRNAseq_reference <- readRDS("./Seurat.YJ.final.rds")
snRNAseq_reference[["roundedDecontX"]] <- CreateAssayObject(counts=decontXcounts(snRNAseq_reference@assays$decontX@counts),
                                                            check.matrix = T)

Idents(snRNAseq_reference) <- "MajorCluster"
#counts <- GetAssayData(snRNAseq_reference, assay = "RNA", slot = "counts")
counts <- snRNAseq_reference@assays$roundedDecontX@counts
cluster <- as.factor(snRNAseq_reference$MajorCluster)
#names(cluster) <- colnames(snRNAseq_reference)
#nUMI <- snRNAseq_reference$nCount_RNA
#names(nUMI) <- colnames(snRNAseq_reference)
#nUMI <- colSums(counts)



#reference <- Reference(counts, cluster, nUMI)
reference <- Reference(counts, cluster)
saveRDS(reference,"./spacexr_reference.rds")
reference <- readRDS("../../spacexr_reference.rds")

rm(nUMI, cluster,counts, snRNAseq_reference)

#dataMat_spatial <- readRDS("tmp.data2.rds")


#rm(data2)
#spatial.obj <- data2
spatial.obj <- readRDS("../data2.prepped.rds")

### use SCT counts
bulk_spatial_list2 <- lapply(names(spatial.obj@images), function(image) {
  #imageOut <- paste0(image,".RCTD")
  imageQuery <- paste0(image,".query.SCT")
  query.counts <- LayerData(spatial.obj, assay = "SCT", layer = "counts")[, Cells(spatial.obj[[image]])]
  coords <- GetTissueCoordinates(spatial.obj, which = "centroids", image = image)
  query <- SpatialRNA(coords, query.counts)
  #rm(coords, query.counts)
  #RCTD <- create.RCTD(query, reference, max_cores = 20)
  #puck <- run.RCTD.local(RCTD, doublet_mode = "full")
  #doubPuck <- run.RCTD.local(RCTD, doublet_mode = "doublet")
  #saveRDS(puck, paste0(image,".RCTD.rds"))
  #saveRDS(doubPuck, paste0(image,".RCTD.doublet.rds"))
  saveRDS(query, file = paste0(imageQuery,"_SCT.rds"))
  assign(imageQuery, query)
  query
})

### rerun as replicates using SCT counts
replicate_names <- names(spatial.obj@images)
group_ids <- c(1,3,3,1,2,1,2,1,2,2)
RCTD.reps2 <- create.RCTD.replicates(spatialRNA.replicates = bulk_spatial_list2, 
                                    reference = reference, 
                                    replicate_names = replicate_names, 
                                    group_ids = group_ids, 
                                    max_cores = 20,
                                    CELL_MIN_INSTANCE=20)


### re-run with SCT
t1 <- Sys.time()
RCTD.reps2 <- run.RCTD.replicates.local(RCTD.reps2, doublet_mode = 'full')
t2 <- Sys.time()
t2-t1
saveRDS(RCTD.reps2, file = "RCTD.reps_SCT.rds")

RCTD.reps <- readRDS("RCTD.reps.rds")
#testing on LBD_6 move into loop
#up number of cores to submit and speed it up, failing due to build location of spacexr
#RCTD <- create.RCTD(query, reference, max_cores = 1)
#update lib paths
#doublet.test <- create.RCTD(query, reference, max_cores = 1)

### run the deconvolution
#doublet_mode: how many cell types predicted per spot
#             "doublet" 1-2 cell types
#             "full"  no restrictions on number of cell types,
#             "multi"  finitely many cell types per pixel, e.g. 3 or 4

#RCTD <- run.RCTD(RCTD, doublet_mode = "full")
#doublet.test <- run.RCTD(doublet.test, doublet_mode = "doublet")


results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
barcodes <- colnames(RCTD@spatialRNA@counts)


resultsdir <- './RCTD_Plots' ## you may change this to a more accessible directory on your computer.
try(dir.create(resultsdir), silent = T)
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)

#plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
#                     results$results_df)

plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
#plot_all_cell_types(results_df = results$weights, spatialRNA@coords, cell_type_names, resultsdir) 

plot_puck_continuous(spatialRNA, barcodes, norm_weights[,'Denate'], ylimit = c(0,0.5), 
                     title ='plot of Dentate weights')



### get cell types

ct_exp2 <- vector("list", length(RCTD.reps2@RCTD.reps))
names(ct_exp2) <- paste0("LBD_",c(1,3:9,11,12))
celltypes <- RCTD.reps2@RCTD.reps[[1]]@cell_type_info$info[[2]]


### generate cell type expression
for(i in 1:length(ct_exp2)){
  print(paste0("sample ", i))
  ct_renorm <- RCTD.reps2@RCTD.reps[[i]]@cell_type_info$renorm[[1]]
  ct_prop <- RCTD.reps2@RCTD.reps[[i]]@results$weights
  ct_mat <- lapply(1:length(celltypes), function(k){
    mat1 <- matrix(ct_renorm[,colnames(ct_renorm) == celltypes[k]], nrow = nrow(ct_renorm), ncol = 1)
    mat2 <- matrix(ct_prop[,colnames(ct_prop) == celltypes[k]], nrow = 1, ncol = nrow(ct_prop))
    mat1 %*% mat2 
  })
  ct_sum <- Reduce("+", ct_mat)
  ct_exp2[[i]] <- vector("list", length(celltypes))
  names(ct_exp2[[i]]) <- celltypes
  for(j in 1:length(celltypes)){
    print(paste0("celltype ", celltypes[j]))
    ct_exp2[[i]][[j]] <- RCTD.reps2@RCTD.reps[[i]]@spatialRNA@counts * (ct_mat[[j]] / ct_sum)
  }
}
saveRDS(ct_exp2, "ct_exp2_SCT_spacexr.rds")
ct_exp2 <- readRDS("spacexr/ct_exp2_SCT_spacexr.rds")
#library(SpatialExperiment)
#spatial.obj <- AddSp

### get normalized weights (total ~=1) added to spatial data
#weightmat <- do.call(rbind, lapply(1:length(RCTD.reps@RCTD.reps), function(i) RCTD.reps@RCTD.reps[[i]]@results$weights))
weightmat2 <- do.call(rbind, lapply(1:length(RCTD.reps2@RCTD.reps), function(i) normalize_weights(RCTD.reps2@RCTD.reps[[i]]@results$weights)))
saveRDS(object = weightmat2, "weightmat2_SCT_normalized.rds")


spatial.obj <- AddMetaData(spatial.obj, metadata = weightmat2)

#### merge sparce matrix and plot
for (celltype in celltypes) {
  dfName <- paste0("adjExp_",celltype)
  tmp <- RowMergeSparseMatrices(ct_exp2$LBD_1[[celltype]], list(ct_exp2$LBD_3[[celltype]],
                                                               ct_exp2$LBD_4[[celltype]],
                                                               ct_exp2$LBD_5[[celltype]],
                                                               ct_exp2$LBD_6[[celltype]],
                                                               ct_exp2$LBD_7[[celltype]],
                                                               ct_exp2$LBD_8[[celltype]],
                                                               ct_exp2$LBD_9[[celltype]],
                                                               ct_exp2$LBD_11[[celltype]],
                                                               ct_exp2$LBD_12[[celltype]]))
  tmpAssay <- CreateAssayObject(data = tmp)
  spatial.obj[[dfName]] <- tmpAssay
}

saveRDS(spatial.obj, file = "../bothCohorts_spatial_withAdjExpAssays_083023.rds")


imageOrder <- c("LBD_1","LBD_5", "LBD_2", "LBD_6", "LBD_3", "LBD_4")
imageNo3 <- c("LBD_1","LBD_5", "LBD_2", "LBD_6", "LBD_4")
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

lapply(celltypes, function(cell) {
  limits <- c(0,max(dataMat_spatial[[cell]]))
  p1<- SpatialFeaturePlot(dataMat_spatial, features = cell, 
                          pt.size.factor = 1.6, 
                          #alpha = c(0.1,0.9),
                          images = imageNo3, 
                          stroke = 0.01,
                          #min.cutoff = "q1", 
                          #max.cutoff = "q99", 
                          combine=F
  )
  names(p1) <- imageNo3
  p2 <- SpatialFeaturePlot(dataMat_spatial, features = cell, 
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
  p1 <- lapply(1:length(p1),function(i) p1[[i]] + labs(title = imageOrder[[i]]) + theme(plot.title = element_text(hjust = 0.5))) 
  #p1 <- unlist(p1,recursive = F)
  m <- grid.arrange(grobs=p1,ncol=length(dataMat_spatial@images),nrow=length(p1)/length(dataMat_spatial@images))
  ggsave(paste0("./figures/spacexrDecon_",cell,"_spatialDistribution_041023.png"), m, limitsize = F, dpi=300, width = 12, height = 8)
})




gps <- VlnPlot(dataMat_spatial, features = colnames(weightmat),group.by = "seurat_clusters",split.by = "disease_state",combine = F)

plotCorrelationMatrix(as.matrix(weightmat))
plotInteractions(as.matrix(weightmat), which = "heatmap",metric = "jaccard",min_prop = 0.02)


### add radius (radial approach to GEX analysis of spots)
coord <- lapply(names(dataMat_spatial@images), function(i) {
  tmp <- GetTissueCoordinates(dataMat_spatial,image = i)
  tmp$slideID <- i
  tmp
})
coord <- do.call(rbind, coord)

### I wrote a function plotInteractionsNeighbours, extend the function across several neighour spots
### neighbours =5, means 5 most adjacent spots
plotInteractionsNeighbours(as.matrix(weightmat), which = "heatmap",metric = "jaccard",min_prop = 0.02,coordmat = coord,neighbours = 5)


