### visualization scripts

library(Seurat)
library(openxlsx)
library(gridExtra)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(VennDiagram)
library(SPOTlight, lib.loc = "./Rlib4.2.2/")
#library(SpatialExperiment)
library(ggplot2)
library(ggridges)
library(cowplot)
library(aplot)
library(ggpubr)
library(patchwork, lib.loc = "./Rlib4.2.2/")
source("./scripts/helper.r")
source("./tools/helper.r")
options(warn = 1, Seurat.object.assay.version = "v5")
library(future)
library(doFuture)
library(grid)
options(future.globals.maxSize=2097152000,future.seed=TRUE)
future::plan("multisession", workers = 10)
plan("sequential")
library(CellChat, lib.loc = "./Rlib4.2.2/")

timeStamp <- format(Sys.time(), "%m%d%y")

setwd("./")
### no Unknown
data2 <- readRDS("bothCohorts_spatial_withAdjExpAssays_dropUnknowns_correctedWeight_fixedImageRotation_051024.rds")
library(imager)
#library(RStoolbox)
library(raster)
imager::grayscale()
### ONLY RUN ONCE
### fix all rotations permenantly
# flip_angle %in% c(180, "R90", "L90", "Hf", "Vf")
data2=rotateSeuratImage(data2,rotation = "L90", slide = "LBD_9")
data2=rotateSeuratImage(data2,rotation = "L90", slide = "LBD_6")
data2=rotateSeuratImage(data2,rotation = "L90", slide = "LBD_8")
data2=rotateSeuratImage(data2,rotation = "L90", slide = "LBD_11")
data2=rotateSeuratImage(data2,rotation = "R90", slide = "LBD_12")
data2=rotateSeuratImage(data2,rotation = "Hf", slide = "LBD_3")


saveRDS(data2, "bothCohorts_spatial_withAdjExpAssays_dropUnknowns_correctedWeight_fixedImageRotation_051024.rds")
source("./scripts/seurat_grayscale_dimplots.r")
gssdp <- SpatialDimPlot.grayscale_imagescales(data2, group.by="CorrectedLayers",cols = layerCols2, imagescales = imagescales)
gs <- wrap_plots(gssdp[imageOrder], ncol = 5, guides = "collect")
saveImages("./figures/051024_Grayscale_spatialDimPlot_CorrectedLayers.png", gs, width = 10, height = 6, dpi = 300, units = "in")

imagescales <- getImagePointSizes(data2)
Idents(data2) <- "LBannotation"
lbSpots <-  WhichCells(data2, idents = "LB(+)")
lbSpotsOnly <- subset(data2, idents = "LB(+)")
lbSurroundOnly <- subset(data2, idents= "LB Surround")
lbNeg <- subset(data2, idents="LB(-)")
lbConNeg <- subset(data2, idents="CLB(-)")

layerCols2 <- c(#"Unknown" = "#33A02C",
                "Layer1" ="#8D405C", 
                "Layer23"= "#E7BDE1",
                "Layer4"= "#CF8CA4",
                "Layer5"= "#9F6E80",
                "Layer6"= "#CDADB9",
                "WM"= "#67A9D8")
layers <- names(layerCols2)
layerColsAPOE <- c(layerCols2,layerCols2)
names(layerColsAPOE) <- c(paste0(names(layerCols2),"_","E3/3"),
                       paste0(names(layerCols2),"_","E3/4"))
layerCols3 <- c(layerCols2,layerCols2,layerCols2,layerCols2,layerCols2,layerCols2)
names(layerCols3) <- c(names(layerCols2),
                       paste0(names(layerCols2),"_","Ctrl"),
                       paste0(names(layerCols2),"_","LBD"),
                       paste0(names(layerCols2),"_","Tri"),
                       paste0(names(layerCols2),"_","E3/3"),
                       paste0(names(layerCols2),"_","E3/4"))

layerCols4 <- c(layerCols2,layerCols2)
names(layerCols4) <- c(paste0(names(layerCols2),".","CON"),
                       paste0(names(layerCols2),".","LBD"))


apoeColors <- c("E3"="#2367AC", 
                "E4"="#B21F2C")

lbAnnoCols <- c("LB(+)" ="green", 
                "LB Surround"= "purple",
                "LB(-)"= "grey",
                "CLB(-)"= "grey40")
lbLayerAnnoCols <- rep(c("LB(+)" ="green", 
                     "LB Surround"= "purple",
                     "LB(-)"= "grey",
                     "CLB(-)"= "grey40"),6)
names(lbLayerAnnoCols) <- unlist(lapply(layers, function(layer) paste0(layer,"_", names(lbAnnoCols))))

data2$layers.anno <- factor(data2$layers.anno, levels = unlist(lapply(layers, function(layer) paste0(layer,"_", names(lbAnnoCols)))))

layerColsAnno <- c(layerCols2,layerCols2,layerCols2,layerCols2)
names(layerColsAnno) <- c(paste0(names(layerCols2),"_",names(lbAnnoCols)[1]),
                          paste0(names(layerCols2),"_",names(lbAnnoCols)[2]),
                          paste0(names(layerCols2),"_",names(lbAnnoCols)[3]),
                          paste0(names(layerCols2),"_",names(lbAnnoCols)[4]))
imageOrder <- c("LBD_1","LBD_7","LBD_6","LBD_8","LBD_3",
                "LBD_5","LBD_9","LBD_11","LBD_12","LBD_4")
lbd_images <- c("LBD_3","LBD_4", "LBD_6","LBD_8","LBD_11", "LBD_12")
sdp <- SpatialDimPlot_imagescales(lbSpotsOnly, imagescales = imagescales, group.by = "CorrectedLayers", imageAlpha = 0.5)
sdp <- wrap_plots(sdp[lbd_images], ncol = 3, guides = "collect")
saveImages("./distance/figures/LBannotation_spatialDimPlot_CorrectedLayers.png", sdp, width = 10, height = 6, dpi = 300, units = "in")

p1 <- SpatialDimPlot_imagescales(seurat.obj = data2, imagescales = imagescales,
                                 group.by="LBannotation",
                                 cols = lbAnnoCols, stroke=0
)
p1 <- p1[imageOrder]
p1$LBD_9 <- p1$LBD_9 + coord_flip() + scale_y_reverse()
p1$LBD_6 <- p1$LBD_6 + coord_flip() + scale_y_reverse()
p1$LBD_8 <- p1$LBD_8 + coord_flip() + scale_y_reverse()
p1$LBD_11 <- p1$LBD_11 + coord_flip() + scale_y_reverse()
p1$LBD_12 <- p1$LBD_12 + coord_flip() + scale_x_reverse()
p1$LBD_3 <- p1$LBD_3 + scale_x_reverse()
p <- wrap_plots(p1, guides = 'collect', nrow = 2) 
saveImages("./figures/051024_LBannotation_spatialDimPlot.png", p, width = 16, height = 8, dpi = 300, units = "in")

Idents(data2) <- "CorrectedLayers"
noWM <- subset(data2, idents=names(layerCols2)[-6])
p1 <- SpatialDimPlot_imagescales(seurat.obj = noWM, imagescales = imagescales,
                                 group.by="LBannotation",
                                 cols = lbAnnoCols, stroke=0
)



p1 <- p1[imageOrder]

p <- wrap_plots(p1, guides = 'collect', nrow = 2) 
saveImages("./figures/051024_LBannotation_noWM_spatialDimPlot.png", p, width = 16, height = 8, dpi = 300, units = "in")

#noWM <- subset(data2, subset=CorrectedLayers != "WM")
p1 <- SpatialDimPlot_imagescales(seurat.obj = lbSpotsOnly, imagescales = imagescales[lbd_images],
                                 group.by="CorrectedLayers",
                                 cols = lbAnnoCols, stroke=0
)
p1 <- p1[imageOrder]
p <- wrap_plots(p1, guides = 'collect', nrow = 2) 
saveImages("./figures/051024_LBspotsOnly_LayerAnnotation_spatialDimPlot.png", p, width = 16, height = 8, dpi = 300, units = "in")


### layer figures
Idents(data2) <- "CorrectedLayers"
cells <- c()
for (name in layers) {
  Idents(data2) <- "CorrectedLayers"
  cells[[name]] <- WhichCells(data2, idents = name)
}
for (i in names(cells)) {
  layer <- i
  cell <- cells[[layer]]
  layerList <- list(tmp=cell)
  names(layerList) <- i
  #colList <- c("1"="green","Unselected"="grey")
  #names(colList) <- c(names(out3)[[i]], "Unselected")
  p1 <- SpatialDimPlot_imagescales(seurat.obj = data2, imagescales = imagescales,
                                   cells.highlight = layerList,
                                   cols.highlight = layerCols2,
  )
  #png(paste0("./figures/",layer, "_byClusterDEGs.png"), width = 12, height = 6, res = 600, units = "in")
  #print(wrap_plots(p1, guides = 'collect', ncol = 2))
  #dev.off()
  p <- wrap_plots(p1, guides = 'collect', ncol = 2) 
  saveImages(paste0("./figures/CorrectedLayers_",layer, ".png"), p, width = 12, height = 5*length(p)/2, dpi = 300, units = "in")
}

######### cell composition
options(warn = 1)
Idents(data2) <- "LBannotation"
DefaultAssay(data2) <- "SCT"
lbSpots <-  WhichCells(data2, idents = "LB(+)")
lbSpotsOnly <- subset(data2, idents = "LB(+)")
layerCols2 <- c("Layer1" ="#8D405C", 
                "Layer23"= "#E7BDE1",
                "Layer4"= "#CF8CA4",
                "Layer5"= "#9F6E80",
                "Layer6"= "#CDADB9",
                "WM"= "#67A9D8")
### barplots/pie for cell composition
library(ggpubr)
library(patchwork)
celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
layers <- names(layerCols2)
weightmat <- lbSpotsOnly@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
write.xlsx(cellPercent, "./outputData/lbSpotsOnly_cellCompositionPerSpot_CorrectedLayerAndSample.xlsx")
#expects a matrix
bp <- c()
for (layerName in layers){
  layerMat <- weightmat[weightmat$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])*100
  sdev <- sapply(layerMat[celltypes], sd)*100
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp,aes(x =celltypes, y=means))+ 
    geom_bar(stat="identity", show.legend = F, color=layerCols2[[layerName]], fill=layerCols2[[layerName]]) + 
    geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    ggtitle(paste0(layerName,"_allSamples")) + ylim(-5,80) +
    ylab("") + xlab("") + theme(axis.text = element_text(size = 20))
  outName <- paste0("./distance/figures/lbSpotOnly_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}
bp<- c()
weightmat <- lbSpotsOnly@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
max(cellPercent[celltypes])
min(cellPercent[celltypes])
for (layerName in layers[-6]){
  layerMat <- cellPercent[cellPercent$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])
  sdev <- sapply(layerMat[celltypes], sd)
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp,aes(x =celltypes, y=means))+ 
    geom_bar(stat="identity", show.legend = F, color=layerCols2[[layerName]], fill=layerCols2[[layerName]]) + 
    geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    #ggtitle(paste0(layerName,"_allSamples")) + 
    ylim(0,100) +
    ylab(layerName) + xlab("") + theme(axis.text = element_text(size = 0))
  #outName <- paste0("./distance/figures/lbSpotOnly_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  #saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}
p <- wrap_plots(bp, ncol = 1) + theme(axis.text.x = element_text(size = 16))
p
saveImages("figures/lbSpotsOnly_barPlot_percentageCelltypesPerLayer.pdf",p, width = 8, height = 6)



#### lb surround
### barplots/pie for cell composition

library(ggpubr)
library(patchwork)
celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
layers <- names(layerCols2)
lsSurr <- subset(data2, idents = "LB Surround")
weightmat <- lsSurr@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
write.xlsx(cellPercent, "./outputData/lbSurroundOnly_cellCompositionPerSpot_layerAndSample.xlsx")
#expects a matrix
bp <- c()
for (layerName in layers){
  layerMat <- weightmat[weightmat$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])*100
  sdev <- sapply(layerMat[celltypes], sd)*100
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp,aes(x =celltypes, y=means))+ 
    geom_bar(stat="identity", show.legend = F, color=layerCols2[[layerName]], fill=layerCols2[[layerName]]) + 
    geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    ggtitle(paste0(layerName,"_allSamples")) + ylim(-5,80) +
    ylab("") + xlab("") + theme(axis.text = element_text(size = 20))
  outName <- paste0("./distance/figures/lbSurroundOnly_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}

#### all spots in layers
### barplots/pie for cell composition
library(ggpubr)
library(patchwork)
celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
layers <- names(layerCols2)
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
write.xlsx(cellPercent, "./outputData/allSpots_cellCompositionPerSpot_CorrectedLayerAndSample.xlsx")
#expects a matrix
bp <- c()
for (layerName in layers){
  layerMat <- weightmat[weightmat$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])*100
  sdev <- sapply(layerMat[celltypes], sd)*100
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp,aes(x =celltypes, y=means))+ 
    geom_bar(stat="identity", show.legend = F, color=layerCols2[[layerName]], fill=layerCols2[[layerName]]) + 
    geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    ggtitle(paste0(layerName,"_allSamples")) + ylim(-5,80) +
    ylab("") + xlab("") + theme(axis.text = element_text(size = 20))
  outName <- paste0("./distance/figures/allSpots_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}

###LB annotations by layer
tmp <- as.data.frame(table(data2$CorrectedLayers, data2$LBannotation, data2$sampleid))
tmp <- split(tmp, ~Var3)
tmp$all <- as.data.frame(table(data2$CorrectedLayers, data2$LBannotation))

layerOrder <- levels(data2$CorrectedLayers)
for (i in names(tmp)){
  tmp2 <- tmp[[i]]
  total <- sum(tmp2$Freq)
  tmp2$percTotalSpots <- (tmp2$Freq/total)*100
  for (x in unique(tmp2$Var1)){
    idx <- rownames(tmp2[tmp2$Var1 == x,])
    tmp2[idx,"percLayerSpots"] <- (tmp2[idx,"Freq"]/sum(tmp2[idx,"Freq"]))*100
  }
  tmp[[i]] <- tmp2
  tmp2 <- tmp2[tmp2$Var2 == "LB(+)",]
  p1 <- ggplot(data = tmp2,aes(x =Var1, y=Freq, fill=Var1))+ 
    geom_bar(stat="identity", show.legend = F) +
    scale_fill_manual(values=layerCols2) +
    ggtitle(paste0("LB(+) by layer- ",i)) +
    ylab("Counts") + xlab("") + theme(axis.text = element_text(size = 14),
                                      axis.title.y = element_text(size = 20), 
                                      title = element_text(size = 20)) + rotate_x_text()
  p2 <- ggplot(tmp2, aes(x=Var1, y=percTotalSpots, fill=Var1))+ 
    geom_bar(stat="identity", show.legend = F) +
    scale_fill_manual(values=layerCols2) +
    #ggtitle(paste0("LB(+) by layer- ",i)) +
    ylab("% total") + xlab("") + theme(axis.text = element_text(size = 14), 
                                       axis.title.y = element_text(size = 20), 
                                       title = element_text(size = 20)) + rotate_x_text()
  p3 <- ggplot(tmp2, aes(x=Var1, y=percLayerSpots, fill=Var1))+ 
    geom_bar(stat="identity", show.legend = F) +
    scale_fill_manual(values=layerCols2) +
    #ggtitle(paste0("LB(+) by layer- ",i)) +
    ylab("% layer") + xlab("") + theme(axis.text = element_text(size = 14), 
                                       axis.title.y = element_text(size = 20), 
                                       title = element_text(size = 20)) + rotate_x_text()
  p <- wrap_plots(list(p1,p2,p3), guides='collect', ncol=3)
  saveImages(paste0("./distance/figures/LBannotation_CorrectedLayers_",i,".pdf"), p, width = 20, height = 5, dpi = 300, units = "in")
}


p <- SpatialDimPlot_imagescales(lbSpotsOnly, group.by="CorrectedLayers", imagescales = imagescales,
                                cols = layerCols2)

p <- wrap_plots(p[lbd_images], ncol = 3, guides = "collect")

saveImages("./distance/figures/LBspots_byCorrectedLayer.png", p)



### lb vs control layer comparison
data2$layerDisease <- factor(paste0(data2$CorrectedLayers,"_",data2$disease_state), 
                             levels=c(paste0(c("Layer1","Layer23","Layer4","Layer5","Layer6","WM"),"_","Ctrl"),
                                      paste0(c("Layer1","Layer23","Layer4","Layer5","Layer6","WM"),"_","LBD"),
                                      paste0(c("Layer1","Layer23","Layer4","Layer5","Layer6","WM"),"_","Tri")))


## layer split by disease state
#Idents(data2) <- "SCT_LayerRes_0.6"
options(warn = 1)
Idents(data2) <- "CorrectedLayers"
noWM <- subset(data2, idents=levels(data2$CorrectedLayers)[-6])
#data2 <- MetaFeature(data2, features = gliosis, meta.name = "Gliosis_metaFeature")
noWM <- AddModuleScore(noWM, features = list(gliosis), name = "GliosisMS")

DefaultAssay(data2) <- "SCT"
sfp <- SpatialFeaturePlot_imagescales(noWM, imagescales = imagescales,
                                      features = "Gliosis_metaFeature")
sfp <- wrap_plots(sfp, ncol = 5, guides = 'collect')
saveImages("./figures/Spatial_GliosisMetaFeature.pdf", sfp)
SpatialFeaturePlot(data2, images = "LBD_3", features = "GliosisMS1", keep.scale = "all")
sfp <- SpatialFeaturePlot_imagescales(noWM, imagescales = imagescales,
                                      features = "GliosisMS1", keep.scale = "all")
sfp <- wrap_plots(sfp, ncol = 5, guides = 'collect')
saveImages(fileName = "./figures/Spatial_GliosisModuleScore_noWM.pdf", plotToSave = sfp)

data2$layersAPOE <- paste0(data2$CorrectedLayers, "_", data2$apoe_geno)

### Violin plots
library(ggpubr)
# + stat_compare_means(label = "p.format", paired = TRUE)?
test_sign <- rstatix::get_comparisons(data2@meta.data, variable = CorrectedLayers)
y_max <- 9

goi_vp <- VlnPlot(data2, 
                  features = "SNCA", 
                  group.by = "CorrectedLayers",
                  #split.by = "apoe_geno",
                  cols = layerCols2, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0, 
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = F, step.increase = .05,p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_SNCA_spotLevel_layers_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)

y_max <- 9.5
goi_vp <- VlnPlot(data2, 
                  features = "APOE", 
                  group.by = "CorrectedLayers",
                  #split.by = "apoe_geno",
                  cols = layerCols2, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_APOE_spotLevel_layers_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)

test_sign <- rstatix::get_comparisons(data2@meta.data, variable = layersAPOE, )
test_sign <- test_sign[c(1,22,39,52,61,66)]
y_max <- 7

goi_vp <- VlnPlot(data2, 
                  features = "SNCA", 
                  group.by = "layersAPOE",
                  #split.by = "apoe_geno",
                  cols = layerCols3, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05, vjust = .8, p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_SNCA_spotLevel_layers_APOEsplit_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)

y_max <- 7
goi_vp <- VlnPlot(data2, 
                  features = "APOE", 
                  group.by = "layersAPOE",
                  #split.by = "apoe_geno",
                  cols = layerCols3, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,vjust = .8, p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_APOE_spotLevel_layers_APOEsplit_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)

test_sign <- rstatix::get_comparisons(data2@meta.data, variable = CorrectedLayers)
goi_vp <- VlnPlot(data2, 
                  features = goi, 
                  group.by = "CorrectedLayers",
                  #split.by = "apoe_geno",
                  cols = layerCols2, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = F)
outName <- "./figures/ViolinPlot_SNCA_APOE_spotLevel_layers.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 9, width = 12)

data2$disease_state <- factor(data2$disease_state)
data2$layers.disease <- factor(paste0(data2$CorrectedLayers, "_", data2$disease_state))

######################## 
### gene expression plots
#### 
#lb annotation split violins with pvals
# lb+
data2$layersAPOE <- factor(data2$layersAPOE, levels = c(paste0(layers,"_E3/3"), paste0(layers,"_E3/4")))
for (dfName in c("lbSpotsOnly", "lbSurroundOnly", "lbNeg", "lbConNeg")){
  df <- get(dfName)
  df$layersAPOE <- droplevels(df$layersAPOE)
  pres <- levels(df$layersAPOE)
  test_sign <- c()
  for (layer in layers){
    if(length(grep(layer, pres)) == 2){
    test_sign[[layer]] <- pres[grep(layer, pres)] 
    }
    }
  y_max <- 7
  goi_vp <- VlnPlot(df, 
                    features = "SNCA", 
                    group.by = "layersAPOE",
                    #split.by = "apoe_geno",
                    cols = layerColsAPOE, 
                    log = F, 
                    same.y.lims = F,
                    assay = "Spatial", 
                    pt.size = 0,
                    y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = F, step.increase = .05, vjust = .6, p.adjust.methods="BH") +
    ggtitle(paste0("SNCA_", dfName))
  outName <- paste0("./figures/",dfName,"_ViolinPlot_SNCA_spotLevel_layers_APOEsplit_pvalues_032624.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)
  
  y_max <- 7
  goi_vp <- VlnPlot(df, 
                    features = "APOE", 
                    group.by = "layersAPOE",
                    #split.by = "apoe_geno",
                    cols = layerColsAPOE, 
                    log = F, 
                    same.y.lims = F,
                    assay = "Spatial", 
                    pt.size = 0,
                    y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = F, step.increase = .06, vjust = .6, p.adjust.methods="BH") +
    ggtitle(paste0("APOE_", dfName))
  outName <- paste0("./figures/",dfName,"_ViolinPlot_APOE_spotLevel_layers_APOEsplit_pvalues_032624.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)
}

droplevels(lbSpotsOnly@meta.data$layersAPOE)
test_sign <- rstatix::get_comparisons(lbSpotsOnly@meta.data, variable = layersAPOE, )
test_sign <- test_sign[c(1,22,39,52,61,66)]
y_max <- 7

goi_vp <- VlnPlot(data2, 
                  features = "SNCA", 
                  group.by = "layersAPOE",
                  #split.by = "apoe_geno",
                  cols = layerCols3, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05, vjust = .8, p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_SNCA_spotLevel_layers_APOEsplit_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)

y_max <- 7
goi_vp <- VlnPlot(data2, 
                  features = "APOE", 
                  group.by = "layersAPOE",
                  #split.by = "apoe_geno",
                  cols = layerCols3, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,vjust = .8, p.adjust.methods="BH")
outName <- "./figures/ViolinPlot_APOE_spotLevel_layers_APOEsplit_pvalues_122123.pdf"
saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)


### Dotplots
dp <- DotPlot(data2, features = rev(unname(unlist(celltypeGenes))), cluster.idents = F) + coord_flip()
saveImages(outName,goi_rp, dpi = 300, units = "in", height = 6, width = 8)

##################
dp <- dp + rotate_x_text(45)
dp$data$Group <- unlist(lapply(dp$data$features.plot, function(x) names(celltypeGenes)[grep(x, celltypeGenes)]))
labels <- ggplot(dp$data, 
                 aes(x = 1, y = features.plot, fill = Group)) + 
  geom_tile() + 
  scale_fill_brewer(palette = 'Set1') + 
  theme_nothing() +
  ylim2(dp)

legend <- plot_grid(get_legend(labels + theme(legend.position = "bottom")))

final <- dp + labels + guide_area() +
  legend + plot_spacer() + plot_spacer() +
  plot_layout(ncol = 3, widths = c(4, -0.1, 1.8), heights = c(4, .5), guides = 'collect')

saveImages("./figures/Dotplot_cellMarkers_spatial.png", final, width = 8, height=7, units='in')


###############################
### split by sample
splitObj <- SplitObject(data2, split.by = "sampleid")

#spatial feature plot LB+ spots by WM
#SpatialDimPlot(data2, image="LBD_1", cells.highlight)


### compare E3/4 vs E3/3 LB+ spot
table(lbSpotsOnly$apoe_geno)

### cell to cell communication between LB+ vs LB surruound
#devtools::install_github("jinworks/CellChat", lib = "./Rlib4.2.2/")
library(CellChat, lib.loc = "./Rlib4.2.2/")
library(patchwork)
Idents(data2) <- "LBannotation"
lbpos_lbsurr <- subset(data2, idents = c("LB(+)","LB Surround"))
color.use <- scPalette(nlevels(lbpos_lbsurr)); names(color.use) <- levels(lbpos_lbsurr)
Seurat::SpatialDimPlot(lbpos_lbsurr, label = T, label.size = 3, cols = color.use, images = "LBD_8")
#get data
data.input = Seurat::GetAssayData(lbpos_lbsurr, layer = "data", assay = "SCT")
#get cell labels
meta = data.frame(labels = Idents(lbpos_lbsurr), row.names = names(Idents(lbpos_lbsurr)))
unique(meta$labels)
# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = Seurat::GetTissueCoordinates(lbpos_lbsurr, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
jsons <- list.files(c("./10xVisiumCytassist_LBDbrain/mrnaseq/230125-NZ/",
                      "/research/bsi/archive/PI/Zhao_Na_m134064/secondary/s309550.10xVisiumCytassist_LBD_2ndcohort/mrnaseq/230726-NZ/"),
                    pattern = 'scalefactors_json.json', recursive = T, full.names = T)
jsons <- jsons[c(-2,-7)]
scale.factorList <- c()
for (json in jsons){
  sample <- gsub("-","_",strsplit(json,'/')[[1]][12])
  scale.factors = jsonlite::fromJSON(txt = json)
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )
  scale.factorList[[sample]] <- scale.factors
}
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
cellchatList <- c()
for (sample in lbd_images) {
  message(sample)
  barcode <- barcodeKey[barcodeKey$id == gsub("LBD_","",x = sample),"Barcode"]
  idx <- grep(paste0("_",barcode,"$"),colnames(data.input))
  data.input.tmp <- data.input[,idx,drop=F]
  #get cell labels
  meta.tmp <- meta[idx,,drop=F]
  spatial.locs.tmp = Seurat::GetTissueCoordinates(lbpos_lbsurr, image=sample, 
                                              scale = NULL, cols = c("imagerow", "imagecol")) 
  cellchat <- createCellChat(object = data.input.tmp, meta = meta.tmp, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs.tmp, 
                             scale.factors = scale.factorList[[sample]])
  print(cellchat)
  cellchat <- updateCellChat(cellchat)
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse
  CellChatDB.use <- CellChatDB
  # use a subset of CellChatDB for cell-cell communication analysis
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
  
  cellchat@DB <- CellChatDB.use
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- updateCellChat(cellchat)
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- updateCellChat(cellchat)
  future::plan("multisession", workers = 10) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- updateCellChat(cellchat)
  if ( dim(cellchat@LR$LRsig)[1] != 0 ) {
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", #"triMean"
                                  trim = 0.1, 
                                  distance.use = TRUE, 
                                  interaction.length = 250, 
                                  scale.distance = 0.01
    )
    
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    cellchat <- computeCommunProbPathway(cellchat)
    
    cellchat <- aggregateNet(cellchat)
    
  }
  
  cellchatList[[sample]] <- cellchat
}
#cellchatList <- readRDS("cellchatList.rds")
saveRDS(cellchatList, "cellchatList.rds")

# Define the cell labels to lift up
#group.new = levels(cellchatList$LBD_3@idents)
#cellchatList <- lapply(cellchatList[-1], function(cellchat){
#  cellchat <- liftCellChat(cellchat, group.new)
#})
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...

data2$layers.anno <- factor(paste0(data2$CorrectedLayers, "_",data2$LBannotation), levels = c(
  paste0(names(layerCols2), "_LB(+)"),
  paste0(names(layerCols2), "_LB Surround"),
  paste0(names(layerCols2), "_LB(-)"),
  paste0(names(layerCols2), "_CLB(-)")))
#saveRDS(data2, "bothCohorts_spatial_withAdjExpAssays_dropUnknowns_110923.rds")

splitObj <- SplitObject(data2, split.by = "sampleid")
namesList <- names(splitObj)

###set up data.input
data.input <- lapply(names(splitObj), function(x) {
  df <- splitObj[[x]]
  df <- Seurat::GetAssayData(df, layer = "data", assay = "SCT")
  df
  })
names(data.input) <- namesList

### get genes common to all 
genes.common <- Reduce(intersect, lapply(data.input, rownames))

### set correct colnames
data.input <- lapply(1:length(data.input), function(x) {
  df <- data.input[[x]]
  tmp <- names(data.input)[[x]]
  colnames(df) <- paste0(tmp,"_", colnames(df))
  df
  })
names(data.input) <- namesList

### bind all together
data.input.red <- lapply(data.input, function(x) x[genes.common,])

data.input.red <- cbind(data.input.red$LBD_1,
                        data.input.red$LBD_3,
                        data.input.red$LBD_4,
                        data.input.red$LBD_5,
                        data.input.red$LBD_6,
                        data.input.red$LBD_7,
                        data.input.red$LBD_8,
                        data.input.red$LBD_9,
                        data.input.red$LBD_11,
                        data.input.red$LBD_12)

# define the meta data
# a column named `slices` should be provided for spatial transcriptomics analysis, which is useful for analyzing cell-cell communication by aggregating multiple slices/replicates. Of note, for comparison analysis across different conditions, users still need to create a CellChat object seperately for each condition.  
meta <- lapply(1:length(splitObj), function(x) {
  df <- splitObj[[x]]
  Idents(df) <- "layers.anno"
  tmp <- names(splitObj)[[x]]
  df<- data.frame(labels = Idents(df), slices = names(splitObj)[[x]])
  rownames(df) <- paste0(tmp,"_", rownames(df))
  df
})

names(meta) <- namesList
meta <- rbind(meta$LBD_1,
              meta$LBD_3,
              meta$LBD_4,
              meta$LBD_5,
              meta$LBD_6,
              meta$LBD_7,
              meta$LBD_8,
              meta$LBD_9,
              meta$LBD_11,
              meta$LBD_12)

#rownames(meta) <- colnames(data.input.red)
#meta$labels <- factor(meta$labels, levels = levels(Idents(seu1)))
meta$slices <- factor(meta$slices, levels = levels(data2$sampleid))
unique(meta$labels) # check the cell labels
unique(meta$slices)

### get scalefactors and whatnot
jsons <- list.files(c("./10xVisiumCytassist_LBDbrain/mrnaseq/230125-NZ/",
                      "/research/bsi/archive/PI/Zhao_Na_m134064/secondary/s309550.10xVisiumCytassist_LBD_2ndcohort/mrnaseq/230726-NZ/"),
                    pattern = 'scalefactors_json.json', recursive = T, full.names = T)
jsons <- jsons[c(-2,-7)]
scale.factorList <- c()
for (json in jsons){
  sample <- gsub("-","_",strsplit(json,'/')[[1]][12])
  scale.factors = jsonlite::fromJSON(txt = json)
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )
  scale.factorList[[sample]] <- scale.factors
}
conversion.factors <- lapply(1:length(scale.factorList), function(x){
  df <- scale.factorList[[x]]
  cf <- df$spot.diameter/df$spot
  data.frame(ratio = cf, tol = df$spot.diameter/2)
}) 
names(conversion.factors) <- namesList
scale.factors <- rbind(conversion.factors$LBD_1,
                            conversion.factors$LBD_3,
                            conversion.factors$LBD_4,
                            conversion.factors$LBD_5,
                            conversion.factors$LBD_6,
                            conversion.factors$LBD_7,
                            conversion.factors$LBD_8,
                            conversion.factors$LBD_9,
                            conversion.factors$LBD_11,
                            conversion.factors$LBD_12)
rownames(scale.factors) <- namesList

### get spatial locs
spatial.locs = lapply(1:length(splitObj), function(x) {
  df <- splitObj[[x]]
  tmp <- names(splitObj)[[x]]
  df <- Seurat::GetTissueCoordinates(df, scale = NULL, image = tmp,
                              cols = c("imagerow", "imagecol")
                               )
  rownames(df) <- paste0(tmp,"_", rownames(df))
  df
  })

names(spatial.locs) <- namesList

spatial.locs <- rbind(spatial.locs$LBD_1,
                      spatial.locs$LBD_3,
                      spatial.locs$LBD_4,
                      spatial.locs$LBD_5,
                      spatial.locs$LBD_6,
                      spatial.locs$LBD_7,
                      spatial.locs$LBD_8,
                      spatial.locs$LBD_9,
                      spatial.locs$LBD_11,
                      spatial.locs$LBD_12)

apoe33 <- c("LBD_1", "LBD_3", "LBD_6", "LBD_7", "LBD_8")
apoe34 <- c("LBD_4", "LBD_5", "LBD_9", "LBD_11", "LBD_12")
idx <- c()
for (i in apoe33){
  i <- paste0("^",i,"_")
  tmp <- grep(i,rownames(spatial.locs))
  idx <- c(idx,tmp)
}
spatial_a33 <- spatial.locs[idx,]
meta_a33 <- meta[idx,]
meta_a33$labels = droplevels(meta_a33$labels, exclude = setdiff(levels(meta_a33$labels),unique(meta_a33$labels)))
meta_a33$slices = droplevels(meta_a33$slices, exclude = setdiff(levels(meta_a33$slices),unique(meta_a33$slices)))
scale_a33 <- scale.factors[apoe33,]
data_a33 <- data.input.red[,idx]


idx2 <- c()
for (i in apoe34){
  i <- paste0("^",i,"_")
  tmp <- grep(i,rownames(spatial.locs))
  idx2 <- c(idx2,tmp)
}
spatial_a34 <- spatial.locs[idx2,]
meta_a34 <- meta[idx2,]
meta_a34$labels = droplevels(meta_a34$labels, exclude = setdiff(levels(meta_a34$labels),unique(meta_a34$labels)))
meta_a34$slices = droplevels(meta_a34$slices, exclude = setdiff(levels(meta_a34$slices),unique(meta_a34$slices)))
scale_a34 <- scale.factors[apoe34,]
data_a34 <- data.input.red[,idx2]

cellchat <- createCellChat(object = data.input.red, meta = meta, group.by = "labels", 
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

cellchat33 <- createCellChat(object = data_a33, meta = meta_a33, group.by = "labels",
                             datatype = "spatial", coordinates = spatial_a33, scale.factors = scale_a33)
cellchat34 <- createCellChat(object = data_a34, meta = meta_a34, group.by = "labels",
                             datatype = "spatial", coordinates = spatial_a34, scale.factors = scale_a34)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat33@DB <- CellChatDB
cellchat34@DB <- CellChatDB

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat33 <- subsetData(cellchat33) # This step is necessary even if using the whole database
cellchat34 <- subsetData(cellchat34) # This step is necessary even if using the whole database


options(future.globals.maxSize=2097152000)
future::plan("multisession", workers = 10)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
saveRDS(cellchat,"cellchat_all.layerAnno.rds")
### apoe split
cellchat33 <- identifyOverExpressedGenes(cellchat33)
cellchat33 <- identifyOverExpressedInteractions(cellchat33)
saveRDS(cellchat33,"cellchat_apoe33.layerAnno.rds")
cellchat34 <- identifyOverExpressedGenes(cellchat34)
cellchat34 <- identifyOverExpressedInteractions(cellchat34)
saveRDS(cellchat34,"cellchat_apoe34.layerAnno.rds")

#cellchat <- readRDS("cellchat_all.layerAnno.rds")
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.knn = TRUE, contact.knn.k = 6, )
saveRDS(cellchat,"cellchat_all_probCalc.layerAnno.rds")

### apoe split
cellchat33 <- readRDS("cellchat_apoe33.layerAnno.rds")
cellchat33 <- computeCommunProb(cellchat33, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.knn = TRUE, contact.knn.k = 6, )
saveRDS(cellchat33,"cellchat_apoe33_probCalc.layerAnno.rds")

cellchat34 <- readRDS("cellchat_apoe34.layerAnno.rds")
cellchat34 <- computeCommunProb(cellchat34, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.knn = TRUE, contact.knn.k = 6, )
saveRDS(cellchat34,"cellchat_apoe34_probCalc.layerAnno.rds")

cellchat33 <- readRDS("cellchat_apoe33_probCalc.layerAnno.rds")
cellchat34 <- readRDS("cellchat_apoe34_probCalc.layerAnno.rds")

#cellchat <- readRDS("cellchat_all_probCalc.layerAnno.rds")
### filter communication
minCells <- reshape2::melt(table(data2$layers.anno))
minCell <- min(minCells[minCells$value > 0, "value"]) -1
#minCell=4
cellchat <- filterCommunication(cellchat, min.cells = minCell)
cellchat33 <- filterCommunication(cellchat33, min.cells = 5)
cellchat34 <- filterCommunication(cellchat34, min.cells = 5)

### get pathway level data
cellchat <- computeCommunProbPathway(cellchat)
cellchat33 <- computeCommunProbPathway(cellchat33)
cellchat34 <- computeCommunProbPathway(cellchat34)


c33 <- unique(cellchat33@LR$LRsig$pathway_name)
c34 <- unique(cellchat34@LR$LRsig$pathway_name)

write.csv(cellchat33@LR$LRsig[cellchat33@LR$LRsig$pathway_name %in% cellchat33@netP$pathways,], "apoe33_significantPathways_layer.LBanno_010924.csv", quote = F)
write.csv(cellchat34@LR$LRsig[cellchat34@LR$LRsig$pathway_name %in% cellchat34@netP$pathways,], "apoe34_significantPathways_layer.LBanno_010924.csv", quote = F)


### aggregate it
cellchat <- aggregateNet(cellchat)

cellchat33 <- aggregateNet(cellchat33)
cellchat34 <- aggregateNet(cellchat34)
# Compute the network centrality scores
cellchat33 <- netAnalysis_computeCentrality(cellchat33, slot.name = "netP")# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat34 <- netAnalysis_computeCentrality(cellchat34, slot.name = "netP")

saveRDS(cellchat33, "cellchat_apoe33_centrality.layerAnno.rds")
saveRDS(cellchat34, "cellchat_apoe34_centrality.layerAnno.rds")

cellchat33 <- readRDS("cellchat_apoe33_centrality.layerAnno.rds")
cellchat34 <- readRDS("cellchat_apoe34_centrality.layerAnno.rds")


df.net33 <- subsetCommunication(cellchat33)
df.net34 <- subsetCommunication(cellchat34)

write.csv(df.net33, "apoe33_significantLigandReceptor_layer.LBanno_010924.csv", quote = F)
write.csv(df.net34, "apoe34_significantLigandReceptor_layer.LBanno_010924.csv", quote = F)


df.netP33 <- subsetCommunication(cellchat33, slot.name = "netP")
df.netP34 <- subsetCommunication(cellchat34, slot.name = "netP")

write.csv(df.netP33, "apoe33_significantPathways_layer.LBanno_010924.csv", quote = F)
write.csv(df.netP34, "apoe34_significantPathways_layer.LBanno_010924.csv", quote = F)


pdf("pathwayBubbleplot_APOE33.pdf", height = 60, width = 14)
netVisual_bubble(cellchat33, thresh = 0.05)
dev.off()
pdf("pathwayBubbleplot_APOE34.pdf", height = 60, width = 14)
netVisual_bubble(cellchat34, thresh = 0.05)
dev.off()


cellchat33.merge <- updateCellChat(cellchat33)
cellchat34.merge <- updateCellChat(cellchat34)

### filtering Layer1 lb /lb surrournd since not present in both
cellchat33.merge <- filterCommunication(cellchat33.merge, min.cells = 20)
cellchat34.merge <- filterCommunication(cellchat34.merge, min.cells = 20)
# Define the cell labels to lift up
group.new = levels(cellchat33.merge@idents)
cellchat34.merge <- liftCellChat(cellchat34.merge, group.new)

cellchat33.merge <- netAnalysis_computeCentrality(cellchat33.merge, slot.name = "netP")
cellchat34.merge <- netAnalysis_computeCentrality(cellchat34.merge, slot.name = "netP")

object.list <- list(APOE33 = cellchat33.merge, APOE34 = cellchat34.merge)
cellchatAPOE2 <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE, merge.data = T)
cellchatAPOE

### computed centrality again...
saveRDS(cellchatAPOE,"cellchatAPOE.merged.041024.rds")

gg1 <- compareInteractions(cellchatAPOE, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchatAPOE, show.legend = F, group = c(1,2), measure = "weight")
pdf("InteractionComparison.pdf")
gg1 + gg2
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
#groupsToAnalyze <- levels(cellchat33.merge@meta$labels)[levels(cellchat33.merge@meta$labels) %in% levels(cellchat34.merge@meta$labels)]
netVisual_diffInteraction(cellchatAPOE, 
                          weight.scale = T, 
                          comparison = c(1,2))

netVisual_diffInteraction(cellchatAPOE, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchatAPOE)
gg2 <- netVisual_heatmap(cellchatAPOE, measure = "weight")
plotList <- gg1 + gg2 
pdf("Differential_signaling_heatmap.pdf",height = 10, width = 12)
draw(plotList, column_title="APOE3/4(red) relative to APOE3/3(blue)")
dev.off()

### control for max weights
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

### show differential interactions between any two cell types
group.cellType <- c(levels(cellchatAPOE@meta$labels))
group.cellType <- factor(group.cellType, levels(cellchatAPOE@meta$labels))
object.list2 <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchatLayers <- mergeCellChat(object.list2, add.names = names(object.list2))
weight.max <- getMaxWeight(object.list2, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_circle(object.list2[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list2)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatLayers, weight.scale = T, measure = "count.merged", label.edge = T, ) 
netVisual_diffInteraction(cellchatLayers, weight.scale = T, measure = "weight.merged", label.edge = T,)
gg1 <- recordPlot()
pdf("Differential_interactions_circleplot_strength.pdf", width = 12 , height = 8)
gg1
title("APOE3/4 (red), APOE3/3 (blue)")
dev.off()

num.link <- sapply(object.list2, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list2)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list2[[i]], title = names(object.list2)[i], weight.MinMax = weight.MinMax)
}
saveImages("Differential_sources_targets_strength.pdf", patchwork::wrap_plots(plots = gg))

### identify signalling changes between conditions
gg <- c()
for (i in levels(droplevels(cellchatAPOE@idents$joint))) {
  print(i)
  try(gg[[i]] <- netAnalysis_signalingChanges_scatter(cellchatAPOE, idents.use = i, xlabel = "", ylabel = "")+ ggtitle(i)) 
}
gg1 <- netAnalysis_signalingChanges_scatter(cellchatAPOE, idents.use = "Layer5_LB(+)")
gg2 <- netAnalysis_signalingChanges_scatter(cellchatAPOE, idents.use = "Layer6_LB(+)")

ylabel <- "Differential incoming interaction strength"
xlabel <- "Differential outgoing interaction strength"
library(ggpubr)
p <- ggplot() + labs(x = xlabel, y = ylabel)
x_axis <- cowplot::get_plot_component(p, "xlab-b") 
y_axis <- cowplot::get_plot_component(p, "ylab-l")
design = "
BA
#C
"
for (i in c("_LB(-)","_LB(+)","_LB Surround","_CLB(-)")){
  ggs <- wrap_plots(gg[grep(i, names(gg), fixed = T)])
  plot <- list(
    ggs,
    y_axis,# B
    x_axis # C
  ) %>%
    wrap_plots() + 
    plot_layout(heights = c(50, 1), widths = c(1, 50), design = design, guides = 'collect')
  out <- gsub(pattern = " ", replacement = "", i)
  saveImages(paste0("Differential_scatter",out,".pdf"), plot, width = 12, height = 8, units = "in", dpi = 300)
  }
patchwork::wrap_plots(plots = gg[grep("_LB(-)", names(gg), fixed = T)], ncol = 3, guides = 'collect') 
patchwork::wrap_plots(plots = gg[grep("_LB(+)", names(gg), fixed = T)], ncol = 2, guides = 'collect') +xlab(xlabel) + ylab(ylabel)
patchwork::wrap_plots(plots = gg[grep("_LB Surround", names(gg), fixed = T)], ncol = 3, guides = 'collect')
patchwork::wrap_plots(plots = gg[grep("_CLB(-)", names(gg), fixed = T)], ncol = 3, guides = 'collect')

### diff interactions between specific "celltypes"
cellchatAPOE <- readRDS("cellchatAPOE.merged.rds")
### identify signaling groups based on functional similarity
cellchatAPOE <- computeNetSimilarityPairwise(cellchatAPOE, type = "functional")
library(umap)
cellchatAPOE <- netEmbedding(cellchatAPOE, type = "functional", umap.method = "uwot") ### usually used umap-learn but can't find it
cellchatAPOE <- netClustering(cellchatAPOE, type = "functional")
pdf("UMAP_SignalingSimilarity_APOEsplit.pdf", height = 10, width = 12)
netVisual_embeddingPairwise(cellchatAPOE, type = "functional", label.size = 3.5)
dev.off()

### show networks with larger or smaller differences between pathways found in both
rankSimilarity(cellchatAPOE, type = "functional")

gg1 <- rankNet(cellchatAPOE, mode = "comparison", stacked = T, do.stat = TRUE, return.data = T, thresh = 0.05)
gg2 <- rankNet(cellchatAPOE, mode = "comparison", stacked = F, do.stat = TRUE)
gg <- gg1 + gg2
saveImages(fileName <- "rankSimilarity_APOE33vsAPOE34.pdf", gg, height = 12, width = 12)
### side by side comparison between datasets
library(ComplexHeatmap)
i = 1
### combining all the identified signaling pathways from different datasets 
#outgoing
pathway.union <- union(object.list$APOE33@netP$pathways, object.list$APOE34@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20)
pdf("SignalingComparison_outgoing_APOEsplit.pdf", height = 10, width = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "GnBu")
pdf("SignalingComparison_incoming_APOEsplit.pdf", height = 10, width = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "OrRd")
pdf("SignalingComparison_all_APOEsplit.pdf", height = 10, width = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf("test.pdf", height = 8, width = 8)
netVisual_embedding(cellchatAPOE, type = "functional", slot.name = "netP")
dev.off()

gg1 <- netVisual_bubble(cellchatAPOE, sources.use = "Layer5_LB(+)", targets.use = "Layer5_LB Surround", comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in APOE34", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchatAPOE, sources.use = "Layer5_LB(+)",targets.use = "Layer5_LB Surround", comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in APOE33", angle.x = 45, remove.isolate = T)
gg1 + gg2

#sources.use = , targets.use = c(5:11),
netVisual_bubble(cellchatAPOE, sources.use = idents.use, targets.use = idents.use,  comparison = c(1, 2), angle.x = 45, 
                 sort.by.source.priority = T, remove.isolate = F, max.dataset = 2, signaling = "ApoE")

pathways.show <- "ApoE"
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
i <-1
plan("sequential")
idents.use <- grep("LB(+)", levels(cellchatAPOE@idents$joint), fixed = T)
levels(cellchatAPOE@idents$joint)[7]
gg1 <- netAnalysis_signalingChanges_scatter(cellchatAPOE, idents.use = "Layer5_LB(+)", signaling.exclude = "Glutamate")

cellchatAPOE <- netAnalysis_computeCentrality(cellchatAPOE, slot.name = "netP")
netAnalysis_signalingRole_network(cellchatAPOE, signaling = "Glutamate", )


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

cellchatAPOE <- computeNetSimilarityPairwise(cellchatAPOE, type = "structural")
cellchatAPOE <- netEmbedding(cellchatAPOE, type = "structural", umap.method = "uwot")
cellchatAPOE <- netClustering(cellchatAPOE, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchatAPOE, type = "structural")
netVisual_embeddingPairwiseZoomIn(cellchatAPOE, type = "structural", nCol = 2)


#netVisual_bubble(cellchatAPOE, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#cellchatAPOE@meta$datasets = factor(cellchat@meta$datasets, levels = c("APOE33", "APOE34")) # set factor level

### modify cellchatAPOE to match expected
for (name in c("coordinates", "scale.factors")){
  cellchatAPOE2@images[[name]] <- rbind(cellchatAPOE2@images$APOE33[[name]],cellchatAPOE2@images$APOE34[[name]])
}

imageOrder <- c("LBD_1","LBD_7","LBD_6","LBD_8","LBD_3",
                "LBD_5","LBD_9","LBD_11","LBD_12","LBD_4")
greenPathways <- c( "APP", "RELN", "SPP1", "MIF", "ApoE", "CSF", "CXCL", "CX3C", "Cholesterol") #"Glutamate",


for (pathways.show in greenPathways){
  #netAnalysis_contribution(cellchatAPOE, signaling = pathways.show)
  #vp <- plotGeneExpression(cellchatAPOE, enriched.only = T, signaling = pathways.show,
   #                split.by = "datasets", colors.ggplot = T, type = "violin", color.use = layerColsAnno)
  interactionName <- extractEnrichedLR(cellchatAPOE, signaling = pathways.show, geneLR.return = T)
  gg <- c()
  for (slice in imageOrder){
    for (pair in interactionName$pairLR[[1]]) {
      gg <- spatialFeaturePlot_groups(cellchatAPOE2, 
                                     #signaling = pathways.show, 
                                     ncol = 1,
                                     #n.colors = 25,
                                     pairLR.use = pair,
                                     slice.use = slice, point.size = 1.3, 
                                     do.binary = TRUE, cutoff = 0.05, show.legend.combined = F,
                                     enriched.only = F, color.heatmap = "Reds", direction = 1)
      command <- paste0("gg + ", adjustment[[slice]])
      gg <- eval(parse(text = command))
      fileName <- paste0("./figuresBySample/BinarySpatialFeaturePlot_cellChat_LR_receptorComplex_",
                         pathways.show,"_",slice,"_",pair,"_041024.pdf")
      try(saveImages(fileName, gg, width = 8, height = 6))
    }
  }
}
  
  
  #gg$LBD_12 <- wrap_plots(gg$LBD_12, guides = 'collect') & coord_flip() & scale_x_reverse() & scale_y_reverse()
  for (slice in names(gg)){
    fileName <- paste0("./figuresBySample/BinarySpatialFeaturePlot_cellChat_LR_receptorComplex_",
                     pathways.show,"_",slice,"_041024.pdf")
    saveImages(fileName, gg[[slice]][[1]], width = 8, height = 8)
  }
  #bfp <- wrap_plots(gg, guides = "collect", ncol = 5)
  #fileName <- paste0("./figures/BinarySpatialFeaturePlot_cellChat_LR_receptorComplex_",
  #                   pathways.show,"_allSamples_041024.pdf")
  #saveImages(fileName, gg[[slice]], width = 4*le, height = 8)
  
  genes <- unique(unlist(lapply(c("APOE33","APOE34"), function(x){
    c(unique(unlist(strsplit(cellchatAPOE@LR[[x]]$LRsig[cellchatAPOE@LR[[x]]$LRsig$pathway_name == pathways.show,"receptor"], 
           split = "_"))),
           unique(unlist(strsplit(cellchatAPOE@LR[[x]]$LRsig[cellchatAPOE@LR[[x]]$LRsig$pathway_name == pathways.show,"ligand"], 
                                  split = "_"))))
  })))
  genes <- genes[genes %in% rownames(data2@assays$SCT$data)]
  for (gene in genes) {
    sfps <- SpatialFeaturePlot_imagescales(data2, imagescales = imagescales, features = gene, 
                                          imageAlpha = 0, alpha=c(.1,1),
                                          stroke=0)
    sfps <- sfps[imageOrder]
    ### correct image rotation
    sfps$LBD_9 <- sfps$LBD_9 + coord_flip() + scale_y_reverse()
    sfps$LBD_6 <- sfps$LBD_6 + coord_flip() + scale_y_reverse()
    sfps$LBD_8 <- sfps$LBD_8 + coord_flip() + scale_y_reverse()
    sfps$LBD_11 <- sfps$LBD_11 + coord_flip() + scale_y_reverse()
    sfps$LBD_12 <- sfps$LBD_12 + coord_flip() + scale_x_reverse()
    sfps$LBD_3 <- sfps$LBD_3 + scale_x_reverse()
    sfp <- wrap_plots(sfps, ncol = 5, guides = "collect")
    fileName <- paste0("./figures/SpatialFeaturePlot_cellChat_LR_",pathways.show,"_",gene,"_041024.pdf")
    saveImages(fileName, sfp, width = 16, height = 8)
    rm(sfp)
    rm(sfps)
  }
  rm(genes)
}

bigSpots <- lapply(imagescales, function(x) { x*1.5 })
gene <- "SNCA"
sfps <- SpatialFeaturePlot_imagescales(data2, imagescales = bigSpots, features = gene, 
                                       imageAlpha = 0, alpha=c(.1,1),
                                       stroke=0)
sfps <- sfps[imageOrder]
### correct image rotation
sfps$LBD_9 <- sfps$LBD_9 + coord_flip() + scale_y_reverse()
sfps$LBD_6 <- sfps$LBD_6 + coord_flip() + scale_y_reverse()
sfps$LBD_8 <- sfps$LBD_8 + coord_flip() + scale_y_reverse()
sfps$LBD_11 <- sfps$LBD_11 + coord_flip() + scale_y_reverse()
sfps$LBD_12 <- sfps$LBD_12 + coord_flip() + scale_x_reverse()
sfps$LBD_3 <- sfps$LBD_3 + scale_x_reverse()
sfp <- wrap_plots(sfps, ncol = 5, guides = "collect")
fileName <- paste0("./figures/SpatialFeaturePlot_cellChat_SNCA_bigSpot_041024.pdf")
saveImages(fileName, sfp, width = 16, height = 8)


################################################
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- "ApoE"
par(mfrow=c(1,1))

netAnalysis_signalingRole_network(cellchat33, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

groupSize <- as.numeric(table(cellchat@idents))
pdf("circleplots_layer.anno.pdf", width = 15, height = 8)
par(mfrow = c(1,2), xpd=TRUE)
nc1 <- netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
                        weight.scale = T, label.edge= F, title.name = "Number of interactions", 
                        arrow.size = .5, vertex.label.cex = 0.7)
nc2 <- netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                        weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", 
                        arrow.width = .5, vertex.label.cex = 0.7)
nc1 + nc2
dev.off()

#saveImages("circleplots_layer.anno.pdf", plots, width = 8, height = 4)
pdf("cellchat_heatmap_count.layer.anno.pdf", width = 8, height = 6)
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
dev.off()
pdf("cellchat_heatmap_weight.layer.anno.pdf", width = 8, height = 6)
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
dev.off()


cellchat@netP$pathways
### write pathways to file
write.csv(cellchat@LR$LRsig[cellchat@LR$LRsig$pathway_name %in% cellchat@netP$pathways,], "significantPathways_layer.LBanno_120723.csv", quote = F)


for (pathways.show in cellchat@netP$pathways) {
  fileName <- paste0(pathways.show, "_circle.pdf")
  par(mfrow=c(1,1))
  pdf(fileName, width = 4, height = 4)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  for (slice in levels(cellchat@meta$slices)) {
    fileName <- paste0(slice,"_",pathways.show, "_spatial.pdf")
    par(mfrow=c(1,1))
    pdf(fileName, width = 12, height = 8)
    nc<- netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = slice, layout = "spatial", 
                        edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0, 
                        #sources.use = "Layer5_LB Surround", targets.use = "Layer5_LB(+)", 
                        #idents.use = c("Layer5_LB(+)","Layer5_LB Surround"), 
                        thresh = 0.001, remove.isolate = F
    )
    print(nc)
    dev.off()
  }
  dev.off()
  
  fileName <- paste0(pathways.show, "_signalingNetwork.pdf")
  pdf(fileName, width = 10, height = 4)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 18, height = 2.5, font.size = 10)
  dev.off()
  
  tmp <- netAnalysis_contribution(cellchat, signaling = pathways.show, return.data = T)
  top.contrib <- rownames(tmp$LR.contribution)[1]
  geneLength <- length(strsplit(top.contrib, "_")[[1]])
  contribLength <- dim(tmp$LR.contribution)[1]
  
  fileName <- paste0(pathways.show, "_LR_contributions.pdf")
  pdf(fileName, width = 6, height = contribLength/2)
  na <- netAnalysis_contribution(cellchat, signaling = pathways.show)
  print(na)
  dev.off()
}

tmp <- netAnalysis_contribution(cellchat, signaling = pathways.show, return.data = T)
top.contrib <- rownames(tmp$LR.contribution)[1]
geneLength <- length(strsplit(top.contrib, "_")[[1]])
par(mfrow=c(1,1))
for (slice in rownames(cellchat@images$scale.factors)){
  fileName <- paste0(slice,"_",pathways.show, "_topContributors_spatialFeaturePlot_LR.pdf")
  pdf(fileName, width = geneLength*4, height = 6)
  spatialFeaturePlot(cellchat, pairLR.use =top.contrib, slice.use = slice, 
                   point.size = 0.3, color.heatmap = "Reds", direction = 1, 
                   do.group = T, enriched.only = T)
  dev.off()
}
  
# Spatial plot
par(mfrow=c(1,1))
pathways.show <-  "TGFb"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = "LBD_4", layout = "spatial", 
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0, 
                    #sources.use = "LB Surround", targets.use = "LB(+)"
                    )

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, slice.use = "LBD_12", layout = "spatial", 
                    edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", 
                    vertex.size.max = 6, vertex.label.cex = 0)
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("TGFB2","TGFBR2","ACVR1B"), slice.use = "LBD_3", 
                   point.size = 0.8, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("TGFB2","TGFBR2","ACVR1B"), slice.use = "LBD_4",
                   point.size = 0.8, color.heatmap = "Reds", direction = 1)
# Take an input of a ligand-receptor pair
spatialFeaturePlot(cellchat, pairLR.use = paste("TGFB2","ACVR1B","TGFBR2",sep="_"), slice.use = "LBD_3", point.size = 0.5,
                   do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = paste("TGFB2","ACVR1B","TGFBR2",sep="_"), slice.use = "LBD_8", point.size = 1.5, 
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

saveRDS(cellchat, file = "cellchat_all_complete.rds")
cellchat <- readRDS("cellchat_all_probCalc.layerAnno.rds")

### pathway vis 
orangePathways <- c("Cholesterol", "MIF", "CSF", "SPP1", "ANGPTL", "APP", "ApoE")
yellowPathways <- c("2-AG", "Glutamate", "TFGb", "NRG", "GRN", "PSAP", "CADM", "CNTN", "L1CAM", "NCAM", "NECTIN", "NOTCH",
                    "SEMA6", "UNC5", "Netrin")
greenPathways <- c("Glutamate", "APP", "RELN", "SPP1", "MIF", "ApoE", "CSF", "CXCL", "CX3C", "Cholesterol")
pathways.show <- "ApoE"
noclb <- levels(data2$layers.anno)[grep("CLB", levels(data2$layers.anno),invert=T)]
for ( pathways.show in greenPathways ){
  ### possible layouts "hierarchy" "chord" "circle" 
  fileName <- paste0("./figures/",pathways.show, "_circlePlot_individual.pdf")
  pdf(fileName, width = 4, height = 6)
  nv <- netVisual_individual(cellchatAPOE, signaling = pathways.show, layout = "hierarchy", top = 1, 
                       sources.use = noclb, targets.use = noclb, remove.isolate = F)
  print(nv)
  dev.off()
  #netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show,  idents.use = "Inflam. DC")
  #netVisual_individual(cellchat, 
                   #sources.use = 4, targets.use = c(5:11),  #comparison = c(1, 2), 
                   #angle.x = 45)
  #pairLR.use = paste("APOE","TREM2",sep="_")
  for (slice in lbd_images) {
    fileName <- paste0("./figures/",slice,"_",pathways.show, "_binary_spatialFeaturePlot_LR.pdf")
    pdf(fileName, width = 4, height = 6)
    sfp <- spatialFeaturePlot(cellchatAPOE, signaling = pathways.show, slice.use = slice, point.size = 1.5, 
                     do.binary = TRUE, cutoff = 0.05, enriched.only = T, color.heatmap = "Reds", direction = 1, )
    print(sfp)
    dev.off()
  }
  
}

#################################################
### cell type genes dotplots
Idents(data2) <- "CorrectedLayers"
for (ass in names(data2@assays)) {
  dp <- DotPlot(data2, features = celltypeGenes, group.by = "CorrectedLayers", 
                assay = "SCT", cluster.idents = F, scale = F, 
                #split.by = "CorrectedLayers" ,
  ) 
  dp <- dp + theme(axis.text.x = element_text(angle = 90)) + ggtitle(ass) + coord_flip() + RestoreLegend()
  dp$data$id <- gsub(pattern = "^.*_", replacement = "", dp$data$id)
  dp$data$colors <- sapply(dp$data$id, function(x) layerCols2[[x]] )   
  print(dp)
  outName <- paste0("./figures/CorrectedLayers_Dotplot_",ass,"_celltypeGenes.pdf")
  saveImages(outName,dp, width = 15, height = 7, dpi=300)
}


###### layer markers
layerMarkers <- c("RELN","SPARC", 
                  "LAMP5", "MGP", 
                  "NEFM", "CARTPT",
                  "CLSTN2","STMN2",
                  "NPY", "SEMA3E", "DIRAS2",
                  "MBP", "PLP1")

common_correctedLayer_markers <- readRDS("DEGs/common_correctedLayer_markers.rds")
commonLayerMarkers <- lapply(common_correctedLayer_markers, function(x) {
  head(x[order(x$combinedRank),"gene.bulk"])
  })


DotPlot(data2, layerMarkers, group.by = "CorrectedLayers") + rotate_x_text()
hm <-DoHeatmap(data2, features = layerMarkers, group.by = "CorrectedLayers", assay = "SCT")
hm <-DoHeatmap(data2, features = unlist(commonLayerMarkers, use.names = F), group.by = "CorrectedLayers", assay = "SCT")

#hm <- DoHeatmap(tmp, features = c(broadClass, broadSecondaryMarker), group.by = "CorrectedLayers")
saveImages("./figures/Heatmap_ourLayerMakers_CorrectedLayers_common_combinedRank.png", hm, width = 10, height = 8, dpi = 300, units = "in")

commonLayerMarkers <- lapply(common_correctedLayer_markers, function(x) {
  head(x[order(x$spot.LFC.order),"gene.bulk"])
})

hm <-DoHeatmap(data2, features = unlist(commonLayerMarkers, use.names = F), group.by = "CorrectedLayers", assay = "SCT")
saveImages("./figures/Heatmap_ourLayerMakers_CorrectedLayers_common_spotLFCrank.png", hm, width = 10, height = 8, dpi = 300, units = "in")

commonLayerMarkers <- lapply(common_correctedLayer_markers, function(x) {
  head(x[order(x$bulk.LFC.order),"gene.bulk"])
})

hm <-DoHeatmap(data2, features = unlist(commonLayerMarkers, use.names = F), group.by = "CorrectedLayers", assay = "SCT")
saveImages("./figures/Heatmap_ourLayerMakers_CorrectedLayers_common_bulkLFCrank.png", hm, width = 10, height = 8, dpi = 300, units = "in")

combined_pub_ours <- c(unlist(commonLayerMarkers, use.names = F),unname(broadClass), broadSecondaryMarker) 
hm <- DoHeatmap(data2, features = combined_pub_ours, group.by = "CorrectedLayers", assay = "SCT")
saveImages("./figures/Heatmap_ourLayerMakers_CorrectedLayers_common_spotLFCrank_publishedMarkers.png", hm, width = 10, height = 8, dpi = 300, units = "in")
rm(hm)

######## 
topLBMarkers <- list("Layer 2/3"= c("PLOD3","G3BP2"),
                  "Layer 4" = c("SLAIN2","MPC2"),
                  "Layer 5" = c("BHMT","SST"),
                  "Layer 6" = c("TRIM59","RPH3A"))



###selected layer markers YJ 032024
layerMarkers <- c("RELN","FOS", "ID3", "AQP4", "GFAP","VIM",
                  "LAMP5", "PCDH8", "HPCAL1", "HOPX", "KCNIP2", "PTK2B",
                  "NEFM", "CARTPT", "SCN1B", "MET", "NEFL", "PVALB",
                  "CLSTN2","SYT1", "SMYD2", "HS3ST2", "CPNE4", "SNCA",
                  "SEMA3E", "HS3ST4", "SST", "FILIP1", "SLC8A1", "OPALIN",
                  "MBP", "PLP1", "SPP1", "MOBP", "CLDN11", "NKX6-2")

hm <- DoHeatmap(data2, features = layerMarkers, group.by = "CorrectedLayers", assay = "SCT")
saveImages("./figures/Heatmap_ourLayerMakers_CorrectedLayers_common_selectedLayerMarkers.png", 
           hm, width = 10, height = 8, dpi = 300, units = "in")

#######################################

### spatial dim plot 
ctl33 <- c("LBD_1", "LBD_7")
ctl34 <- c("LBD_5", "LBD_9")
lbd33 <- c("LBD_6", "LBD_8")
lbd34 <- c("LBD_11", "LBD_12")
tri34 <- c("LBD_3", "LBD_4")
sampleOrder <- c(ctl33, lbd33, tri34[[1]], ctl34, lbd34, tri34[[2]])


Idents(data2) <- "CorrectedLayers"
plots <- SpatialDimPlot_imagescales(data2, imagescales = imagescales, cols = layerCols2)
plots <- plots[sampleOrder]

p <- wrap_plots(plots, ncol = 5, guides = 'collect')
saveImages("figures/LayerSpatialPlot_CorrectedLayers_dropUnknowns_121923.png", p, width = 12, height = 6, dpi = 600 )

### stacked barplot of deconvoluted cell types by layers (also check celltrek script)
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100

toPlotCellPerc <- reshape2::melt(weightmat)
toPlotCellPerc <- aggregate(toPlotCellPerc$value, list(toPlotCellPerc$CorrectedLayers,toPlotCellPerc$variable), FUN=mean)
names(toPlotCellPerc) <- c("Layers","Celltypes","Proportion")
p <- ggplot(toPlotCellPerc, aes(x=Celltypes, y=Proportion, fill=Layers)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=layerCols2) + ggtitle("Mean (Cell type per spot)")
p
saveImages("figures/CellProportion_mean_CorrectedLayers_dropUnknowns_122023.png", p, width = 6, height = 4, dpi = 600 )

toPlotCellPerc <- reshape2::melt(weightmat)
toPlotCellPerc <- aggregate(toPlotCellPerc$value, list(toPlotCellPerc$CorrectedLayers,toPlotCellPerc$variable), FUN=median)
names(toPlotCellPerc) <- c("Layers","Celltypes","Proportion")
p <- ggplot(toPlotCellPerc, aes(x=Celltypes, y=Proportion, fill=Layers)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=layerCols2) + ggtitle("Median (Celltype per spot)")
p
saveImages("figures/CellProportion_median_CorrectedLayers_dropUnknowns_122023.png", p, width = 6, height = 4, dpi = 600 )


### get greyscale
## set up scanpy images
### barplots/pie for cell composition
library(ggpubr)
library(patchwork)
celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
layers <- levels(data2$CorrectedLayers)
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
write.xlsx(cellPercent, "./outputData/allspots_cellCompositionPerSpot_CorrectedLayerAndSample.xlsx")
#expects a matrix
bp <- c()
for (layerName in layers){
  layerMat <- weightmat[weightmat$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])*100
  sdev <- sapply(layerMat[celltypes], sd)*100
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp,aes(x =celltypes, y=means))+ 
    geom_bar(stat="identity", show.legend = F, color=layerCols2[[layerName]], fill=layerCols2[[layerName]]) + 
    geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    ggtitle(paste0(layerName,"_allSamples")) + ylim(-5,80) +
    ylab("") + xlab("") + theme(axis.text = element_text(size = 20))
  outName <- paste0("./distance/figures/lbSpotOnly_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  ggsave(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}



vlList <- list()
yMax <- max(FetchData(data2, vars = "SNCA", slot = 'counts')) +30
lbCol <- c("LB+"="red","LB_surround"="blue","NoLB"="green","NegControl"="yellow")


for (celltype in celltypes){
  vlList[[celltype]] <- VlnPlot(data2, assay = paste0("adjExp_",celltype), 
                                features = "SNCA", split.by = "apoe_geno",
                                group.by = "LBannotation",  
                                layer = "data", cols = lbAnnoCols, 
                                pt.size = 0, split.plot = T
  )
}
vlList <- lapply(celltypes, function(i){
  vlList[[i]] + FontSize(x.title = 0, main = 0, x.text = 0, y.title = 8, y.text = 4) + ylab(i)
})
p <- wrap_plots(vlList, ncol=1, guides = 'collect', widths = 5, heights = 1) + FontSize(x.text = 8)
pdf("./figures/ViolinPlot_layers_byCelltype_compressed.pdf",width = 8, height = 10)
p
dev.off()

library(ggpubr)
vlList <- list()
layers <- names(layerCols2)
Idents(data2) <- "CorrectedLayers"
for ( celltype in celltypes ){
  #tmp <- subset(data2, idents = layer)
  vlList[[celltype]] <- VlnPlot(data2, 
                             features = celltype,
                             cols = layerCols2, 
                             pt.size = 0,
                             y.max = 1.3
                             ) + NoLegend()
}
vlList2 <- lapply(celltypes, function(i){
  vlList[[i]] + FontSize(x.title = 0, main = 0, 
                         y.title = 20, y.text = 4) + ylab(i)  + theme(axis.title.y = element_text(angle = 0))
})
p <- wrap_plots(vlList2, ncol=1, widths = 5, heights = 2) & FontSize(x.text = 0)
p
vp <- VlnPlot(data2, split.by = "CorrectedLayers", 
        features = celltypes, 
        stack = T, 
        cols = layerCols2, 
        same.y.lims = T, 
        flip = T) + stat_summary(fun = mean, geom='point', shape=95, size=10) +
  NoLegend() + 
  xlab(NULL) + ylab("Percentage per spot")

saveImages("figures/ViolinPlot_percentageCelltypesPerLayer.pdf",vp, width = 8, height = 6) 

library(reshape2)
celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
layers <- names(layerCols2)
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
write.xlsx(cellPercent, "./outputData/allSpots_cellCompositionPerSpot_CorrectedLayerAndSample.xlsx")
#expects a matrix

bp <- c()
for (layerName in layers){
  layerMat <- cellPercent[cellPercent$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  tmp <- melt(layerMat)
  colnames(tmp) <- c("Celltypes", "Percentage")
  tmp$Percentage <- tmp$Percentage/max(tmp$Percentage)
  #tmp <- t(layerMat)
  #means <- colMeans(layerMat[celltypes])*100
  #sdev <- sapply(layerMat[celltypes], sd)*100
  #tmp <- as.data.frame(list(means=means,sd=sdev))
  #tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp, aes(x = Celltypes, y = Percentage))+ 
    geom_violin(#color=layerCols2[[layerName]], 
                fill=layerCols2[[layerName]], draw_quantiles = 0.5, trim = T, scale = "width") + 
    #geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    #ggtitle(paste0(layerName,"_allSamples")) + 
    ylim(0,1) +
    ylab(layerName) + xlab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
  #utName <- paste0("./distance/figures/allSpots_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  #saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}

p <- wrap_plots(bp, ncol = 1) + theme(axis.text.x = element_text(size = 16))
saveImages("figures/ViolinPlot_percentageCelltypesPerLayer_version2.pdf",p, width = 8, height = 6)


celltypes <- c("Ast","Endo","Ex","In","Mic","OPC","Olig","Peri")
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "LBannotation", "apoe_geno")]
cellPercent <- weightmat[,celltypes]*100

cellPercent <- cbind(cellPercent, weightmat[,c("CorrectedLayers", "LBannotation", "apoe_geno")])
head(colMaxs(as.matrix(cellPercent[,celltypes]), useNames = F, rows = T))
bp <- c()
annos <- levels(data2$LBannotation)
for (lbStatus in annos){
  layerMat <- cellPercent[cellPercent$CorrectedLayers == lbStatus,]
  layerMat <- layerMat[-c(9,10)]
  tmp <- melt(layerMat)
  colnames(tmp) <- c("Celltypes", "Percentage")
  tmp$Percentage <- tmp$Percentage/max(tmp$Percentage)
  #tmp <- t(layerMat)
  #means <- colMeans(layerMat[celltypes])*100
  #sdev <- sapply(layerMat[celltypes], sd)*100
  #tmp <- as.data.frame(list(means=means,sd=sdev))
  #tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp, aes(x = Celltypes, y = Percentage))+ 
    geom_violin(#color=layerCols2[[layerName]], 
      fill=layerCols2[[layerName]], draw_quantiles = 0.5, trim = T, scale = "width") + 
    #geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    #ggtitle(paste0(layerName,"_allSamples")) + 
    ylim(0,1) +
    ylab(layerName) + xlab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
  #utName <- paste0("./distance/figures/allSpots_BarPlot_CellComposition_spotLevel_",layerName,"_CorrectedLayers.pdf")
  #saveImages(outName, bp[[layerName]], units = "in", height = 4, width = 6)
}


weightmat <- lbSpotsOnly@meta.data[c(celltypes,"CorrectedLayers", "sampleid", "apoe_geno", "LBannotation")]
cellPercent <- weightmat
cellPercent[celltypes] <- cellPercent[celltypes]*100
max(cellPercent[celltypes])
min(cellPercent[celltypes])
for (layerName in layers[-6]){
  layerMat <- cellPercent[cellPercent$CorrectedLayers == layerName,]
  layerMat <- layerMat[-c(9,10)]
  #tmp <- t(layerMat)
  means <- colMeans(layerMat[celltypes])
  sdev <- sapply(layerMat[celltypes], sd)
  tmp <- as.data.frame(list(means=means,sd=sdev))
  tmp$celltypes <- rownames(tmp)
  bp[[layerName]] <- ggplot(data = tmp, aes(x = Celltypes, y = Percentage))+ 
    geom_violin(#color=layerCols2[[layerName]], 
      fill=layerCols2[[layerName]], draw_quantiles = 0.5, trim = T, scale = "width") + 
    #geom_errorbar( aes(x=celltypes, ymin=means-sd, ymax=means+sd), width=0.4, color="black", alpha=0.9, linewidth=1.3) +
    #ggtitle(paste0(layerName,"_allSamples")) + 
    ylim(0,1) +
    ylab(layerName) + xlab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
}
p <- wrap_plots(bp, ncol = 1) + theme(axis.text.x = element_text(size = 16))
p
saveImages("figures/lbSpotsOnly_barPlot_percentageCelltypesPerLayer.png",p, width = 8, height = 6)


weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid", "apoe_geno", "LBannotation", "disease_state")]
#tmp1 <- weightmat[,celltypes]
#tmp1[tmp1 < 0.05] <- 0
#weightmat[,celltypes] <- tmp1

cellPercent <- weightmat
cellPercent[,celltypes] <- cellPercent[,celltypes]*100
#tmp1 <- cellPercent[,celltypes]
#tmp1[tmp1 > 100] <- 100
#cellPercent[,celltypes] <- tmp1

AnnoCols <- c("LB(+)" ="green", 
                "LB Surround"= "purple",
                "LB(-)"= "grey",
                "CLB(-)"= "grey30")
apoeCols <- c("E3/3" ="blue", 
                "E3/4"= "orange")
diseaseCols <- c("Ctrl" ="red", 
              "Tri"= "blue",
              "LBD" = "purple")

bp <- c()
bp_anno <- c()
celltypes <- c("Ex","In", "Olig","OPC","Ast","Mic","Endo","Peri")
for (celltype in celltypes){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation")]
  #cellMat <- cellMat[!cellMat$LBannotation == "CLB(-)",]
  #cellMat <- cellMat[!cellMat$CorrectedLayers == "WM",]
  #cellMat <- cellMat[cellMat[[celltype]] > 1, ]
  bp <- ggplot(cellMat, aes(x=CorrectedLayers, y=cellMat[,celltype])) + 
    geom_boxplot(aes(fill=LBannotation), alpha=1, outlier.shape = NA)  + facet_wrap(facets = ~apoe_geno)  +
    ylim(0,100) +
    scale_fill_manual(values = lbAnnoCols) + 
    ylab(celltype) + xlab("") #+ theme(axis.text = element_blank(), axis.ticks = element_blank())
  filename <- paste0("./figures/Boxplot_deconvolutedSpots_", celltype,"_allAnnotations_splitAPOE_byLayers_032624.pdf")
  saveImages(fileName = filename, bp, width = 8, height = 4)
  
  bp_anno[[celltype]] <- lapply(names(lbAnnoCols), function(anno){
    annoMat <- cellMat[cellMat$LBannotation == anno,]
    bp_anno[[celltype]][[anno]] <- ggplot(annoMat, aes(x=CorrectedLayers, y=annoMat[,1])) + 
      geom_boxplot(aes(fill=apoe_geno), alpha=1, outlier.shape = NA)  + 
      #ylim(0,100) +
      scale_fill_manual(values = apoeCols) + 
      ylab(paste0(celltype, "_", anno)) + xlab("") #+ theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_blank())
  })
  names(bp_anno[[celltype]]) <- names(lbAnnoCols)
}

lbAnnoCols <- c("LB(+)" ="green", 
                "LB Surround"= "purple",
                "LB(-)"= "grey",
                "CLB(-)"= "grey40")
layers <- c("L1","L23","L4","L5","L6","WM")
### change level order to change plot ordering
cellPercent$LBannotation <- factor(cellPercent$LBannotation, levels = c("LB(+)", "LB Surround", "LB(-)","CLB(-)"))
bp <- c()
cellPercent$CorrectedLayers <- factor(gsub("ayer","",cellPercent$CorrectedLayers), levels = layers)
cellPercent$layer.apoe <- paste0(cellPercent$CorrectedLayers, ".", cellPercent$apoe_geno)
cellPercent$layer.apoe <- factor(cellPercent$layer.apoe, levels = paste0(levels(cellPercent$CorrectedLayers), ".",c(rep("E3/3",6),rep("E3/4",6))))
cellPercent$layer.apoe.anno <- paste0(cellPercent$layer.apoe, ".", cellPercent$LBannotation)
#cellPercent$layer.apoe.anno <- factor(cellPercent$layer.apoe.anno, levels = c(paste0(levels(cellPercent$layer.apoe),".",("LB(+)", "LB Surround", "LB(-)","CLB(-)"))

### plot too many versions of deconvoluted cell type figures
for (celltype in celltypes){
  print(celltype)
  name <- paste0("bp_", celltype)
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "layer.apoe", "layer.apoe.anno")]
  #limitList <- quantile(cellMat[,celltype], c(0, 0.99))
  #test_sign <- rstatix::get_comparisons(cellMat,layer.apoe.anno)
  #test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = ".")[[1]][1] == strsplit(x[2], split = ".")[[1]][1]))]
  #command <- eval(parse(text= paste0(celltype, " ~ layer.apoe.anno")))
  #stats <- compare_means(command ,data = cellMat, p.adjust.method = "BH")
  #stats <- stats[gsub("\\..*", "", stats$group1, perl = T) == gsub("\\..*", "", stats$group2, perl = T) ,]
  #cellMat <- cellMat[!cellMat$LBannotation == "CLB(-)",]
  #cellMat <- cellMat[!cellMat$CorrectedLayers == "WM",]
  #cellMat <- cellMat[cellMat[[celltype]] > 1, ]
  #command <- eval(parse(text= paste0("scale_y_continuous(limits =c(0,",unname(limitList)[2],"))")))
  bp <- ggplot(cellMat, aes(x=CorrectedLayers, y=eval(parse(text=celltype)))) + 
    geom_boxplot(aes(fill=LBannotation), alpha=1, 
                 outlier.shape = NA)  + 
    #command + 
    scale_y_continuous(limits = c(0,unname(quantile(cellMat[,celltype],.99)))) +
    facet_wrap(facets = ~apoe_geno)  +
    #ylim(0,100) +
    scale_fill_manual(values = lbAnnoCols) + rotate_x_text()+
    ylab(celltype) + xlab("") + 
    theme(panel.grid = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line())
    #theme(axis.text.x = element_blank(), axis.ticks.x =element_blank(), strip.text.x = element_blank())
  filename <- paste0("./figures/",timeStamp,"_Boxplot_deconvolutedSpots_", celltype,"_allAnnotations_statsAPOE_byLayers.pdf")
  saveImages(fileName = filename, bp, width = 8, height = 4)
  assign(name, bp)
}


bps <- lapply(celltypes, function(celltype){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "layer.apoe", "layer.apoe.anno")]
  ggplot(cellMat, aes(x=CorrectedLayers, y=eval(parse(text=celltype)))) + 
    geom_boxplot(aes(fill=LBannotation), alpha=1, 
                 outlier.shape = NA)  + 
    facet_wrap(facets = ~apoe_geno)  +
    scale_y_continuous(limits = c(0,unname(quantile(cellMat[[celltype]],.99)))) +
    scale_fill_manual(values = lbAnnoCols) + rotate_x_text()+
    ylab(celltype) + xlab("") + theme(panel.grid = element_blank(),
                                      panel.background = element_blank(), 
                                      axis.line = element_line(),
                                      axis.text.x = element_blank(), 
                                      axis.ticks.x =element_blank(), 
                                      strip.text.x = element_blank())
  
})
#bp <- lapply(bp, function(x) x +theme(axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text = element_blank()))

p <- wrap_plots(bps, ncol = 1, guides = "collect") + theme(axis.text.x = element_text(size = 16))
p[[1]] <- p[[1]] + theme(strip.text.x = element_text(size=16))

filename <- paste0("./figures/",timeStamp,"_Boxplot_deconvolutedSpots_allAnnotations_splitAPOE_byLayers.pdf")
saveImages(fileName = filename, p, width = 8, height = 10)

statsList <- c()
for (celltype in celltypes){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "layer.apoe", "layer.apoe.anno","sampleid")]
  cellMat$apoe.anno <- paste0(cellMat$apoe_geno, ".", cellMat$LBannotation)
  cellMat$layer.anno <- paste0(cellMat$CorrectedLayers, ".", cellMat$LBannotation)
  tmp <- reshape2::melt(cellMat)
  colnames(tmp)[grep("variable", colnames(tmp))] <- "Celltypes"
  tmp2 <- tmp %>% group_by(sampleid, CorrectedLayers, Celltypes, LBannotation) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  pvals <- tmp2 %>%
    group_by(CorrectedLayers) %>%
    wilcox_test(Mean_Frequency ~ LBannotation) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  #names(yposition) <- celltypes
  #pvals <- pvals %>% add_xy_position(fun = "max", x = "Celltypes", group = "apoe_geno", dodge = 0.8, step.increase = 0.05)
  #pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position)
  statsList[[celltype]][["layer"]] <- pvals
  print(paste0(celltype," Layer: ",length(grep("ns", pvals$p.signif, invert = T))))
  ### within apoe annotation comparison
  comparisons <- lapply(split(tmp, tmp$layer.apoe), function(x) if(any(length(table(x$LBannotation)) < 2) || any(table(x$LBannotation) < 2)) 
    NULL else table(x$apoe_geno)
  )
  comparisons <- names(rlist::list.clean(comparisons))
  tmp2 <- tmp %>% group_by(sampleid, layer.apoe, Celltypes, LBannotation) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  pvals <- tmp2[tmp2$layer.apoe %in% comparisons,] %>%
    group_by(layer.apoe) %>%
    wilcox_test(Mean_Frequency ~ LBannotation) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  statsList[[celltype]][["layer.apoe"]] <- pvals
  print(paste0(celltype," LayerApoe: ",length(grep("ns", pvals$p.signif, invert = T))))
  
  ### across apoe within anno
  ### drop missing comparison 
  library(dplyr)
  comparisons <- lapply(split(tmp, tmp$layer.anno), function(x) if(any(length(table(x$apoe_geno)) < 2) || any(table(x$apoe_geno) < 2)) 
    NULL else table(x$apoe_geno)
    )
  comparisons <- names(rlist::list.clean(comparisons))
  tmp2 <- tmp %>% group_by(sampleid, layer.anno, Celltypes, apoe_geno) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  pvals <- tmp2[tmp2$layer.anno %in% comparisons,] %>%
    group_by(layer.anno) %>%
    wilcox_test(Mean_Frequency ~ apoe_geno) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p") 
  statsList[[celltype]][["layer.anno"]] <- pvals
  print(paste0(celltype," LayerAnno: ",length(grep("ns", pvals$p.signif, invert = T))))
  
}

wb <- createWorkbook()
for (celltype in names(statsList)) {
  for (comparison in names(statsList[[celltype]])){
    name <- paste0(celltype, ".", comparison)
    addWorksheet(wb, name)
    writeData(wb, sheet = name,statsList[[celltype]][[comparison]],rowNames=TRUE,colNames=TRUE)
  }
}

saveWorkbook(wb, "./figures/092624_pseudobulkStats_cellPercentageComparisons_wilcoxtest_layer.apoe.anno.xlsx", overwrite = T)

condCols <- c("CON" ="orange", 
               "LBD"= "purple")
data2$condition <- ifelse(data2$disease_state == "Ctrl", "CON", "LBD")
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid", "apoe_geno", "LBannotation", "condition")]

cellPercent <- weightmat
cellPercent[,celltypes] <- cellPercent[,celltypes]*100
cellPercent$CorrectedLayers <- gsub("ayer","",cellPercent$CorrectedLayers)
bps <- lapply(celltypes, function(celltype){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "condition")]
  ggplot(cellMat, aes(x=CorrectedLayers, y=eval(parse(text=celltype)))) + 
    geom_boxplot(aes(fill=condition), alpha=1, 
                 outlier.shape = NA)  + 
    facet_wrap(facets = ~apoe_geno)  +
    scale_y_continuous(limits = c(0,unname(quantile(cellMat[[celltype]],.99)))) +
    scale_fill_manual(values = condCols) + rotate_x_text()+
    ylab(celltype) + xlab("") + theme(panel.grid = element_blank(),
                                      panel.background = element_blank(), 
                                      axis.line = element_line(),
                                      axis.text.x = element_blank(), 
                                      axis.ticks.x =element_blank(), 
                                      strip.text.x = element_blank())
})

#bp <- lapply(bp, function(x) x +theme(axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text = element_blank()))

p <- wrap_plots(bps, ncol = 1, guides = "collect") + theme(axis.text.x = element_text(size = 16))
p[[1]] <- p[[1]] + theme(strip.text.x = element_text(size=16))

filename <- paste0("./figures/",timeStamp,"_Boxplot_deconvolutedSpots_condition_splitAPOE.pdf")
saveImages(fileName = filename, p, width = 8, height = 10)

condCols <- c("CON" ="orange", 
              "LBD"= "purple")
data2$condition <- ifelse(data2$disease_state == "Ctrl", "CON", "LBD")
weightmat <- data2@meta.data[c(celltypes,"CorrectedLayers", "sampleid", "apoe_geno", "LBannotation", "condition")]

cellPercent <- weightmat
cellPercent[,celltypes] <- cellPercent[,celltypes]*100
cellPercent$CorrectedLayers <- gsub("ayer","",cellPercent$CorrectedLayers)
library(ggpubr)
bps <- lapply(celltypes, function(celltype){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "condition")]
  ggplot(cellMat, aes(x=CorrectedLayers, y=eval(parse(text=celltype)))) + 
    geom_boxplot(aes(fill=condition), alpha=1, 
                 outlier.shape = NA)  + 
    #facet_wrap(facets = ~apoe_geno)  +
    scale_y_continuous(limits = c(0,unname(quantile(cellMat[[celltype]],.99)))) +
    scale_fill_manual(values = condCols) + rotate_x_text()+
    ylab(celltype) + xlab("") + theme(panel.grid = element_blank(),
                                      panel.background = element_blank(), 
                                      axis.line = element_line(),
                                      axis.text.x = element_blank(), 
                                      axis.ticks.x =element_blank(), 
                                      strip.text.x = element_blank())
})

#bp <- lapply(bp, function(x) x +theme(axis.text.x = element_blank(), axis.ticks = element_blank(), strip.text = element_blank()))

p <- wrap_plots(bps, ncol = 1, guides = "collect") + theme(axis.text.x = element_text(size = 16))
#p[[1]] <- p[[1]] + theme(strip.text.x = element_text(size=16))

filename <- paste0("./figures/",timeStamp,"_Boxplot_deconvolutedSpots_condition.pdf")
saveImages(fileName = filename, p, width = 8, height = 10)


bps <- lapply(celltypes, function(celltype){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "condition")]
  bp <- ggplot(cellMat, aes(x=CorrectedLayers, y=eval(parse(text=celltype)))) + 
    geom_boxplot(aes(fill=condition), alpha=1, 
                 outlier.shape = NA)  + 
    facet_wrap(facets = ~apoe_geno)  +
    scale_y_continuous(limits = c(0,unname(quantile(cellMat[[celltype]],.99)))) +
    scale_fill_manual(values = condCols) + rotate_x_text()+
    ylab(celltype) + xlab("") + theme(panel.grid = element_blank(),
                                      panel.background = element_blank(), 
                                      axis.line = element_line()) + rotate_x_text()
  filename <- paste0("./figures/",timeStamp,"_Boxplot_deconvolutedSpots_", celltype,"_condition_splitAPOE_byLayers.pdf")
  saveImages(fileName = filename, bp, width = 8, height = 4)
})


ggplot(cellMat, aes(x=layer.con, y=eval(parse(text=celltype)))) + 
  geom_boxplot(aes(fill=apoe_geno), alpha=1, #scale = "width", na.rm = T, draw_quantiles = c(.5)
               #outlier.shape = 21
              )  + 
  #facet_wrap(facets = ~apoe_geno)  +
  scale_y_continuous(limits = c(0,unname(quantile(cellMat[[celltype]],.99)))) +
  scale_fill_manual(values = condCols) + rotate_x_text()+
  ylab(celltype) + xlab("") + theme(panel.grid = element_blank(),
                                    panel.background = element_blank(), 
                                    axis.line = element_line()) + rotate_x_text()

library(tidyr)
library(tidyverse)
library(rstatix)
statsList <- c()
for (celltype in celltypes){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "condition","sampleid")]
  cellMat$apoe.con <- paste0(cellMat$apoe_geno, ".", cellMat$condition)
  cellMat$layer.con <- paste0(cellMat$CorrectedLayers, ".", cellMat$condition)
  cellMat$layer.apoe <- paste0(cellMat$CorrectedLayers, ".", cellMat$apoe_geno)
  tmp <- reshape2::melt(cellMat)
  colnames(tmp)[grep("variable", colnames(tmp))] <- "Celltypes"
  ### get means from each sample for stats
  #tmp3 <- tmp %>% group_by(sampleid, disease_state, CorrectedLayers, Celltypes) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  tmp2 <- tmp %>% group_by(sampleid, CorrectedLayers, Celltypes, condition) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  pvals <- tmp2 %>%
    group_by(CorrectedLayers) %>%
    wilcox_test(Mean_Frequency ~ condition) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  #names(yposition) <- celltypes
  #pvals <- pvals %>% add_xy_position(fun = "max", x = "Celltypes", group = "apoe_geno", dodge = 0.8, step.increase = 0.05)
  #pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position)
  statsList[[celltype]][["layer"]] <- pvals
  
  ### within apoe annotation comparison
  tmp2 <- tmp %>% group_by(sampleid, layer.apoe, Celltypes, condition) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  pvals <- tmp2 %>%
    group_by(layer.apoe) %>%
    wilcox_test(Mean_Frequency ~ condition) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  statsList[[celltype]][["layer.apoe"]] <- pvals
  
  ### across apoe within anno
  tmp2 <- tmp %>% group_by(sampleid, layer.con, apoe_geno, Celltypes) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  pvals <- tmp2 %>%
    group_by(layer.con) %>%
    wilcox_test(Mean_Frequency ~ apoe_geno) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p") 
  statsList[[celltype]][["layer.con"]] <- pvals
}

wb <- createWorkbook()
for (celltype in names(statsList)) {
  for (comparison in names(statsList[[celltype]])){
    name <- paste0(celltype, ".", comparison)
    addWorksheet(wb, name)
    writeData(wb, sheet = name,statsList[[celltype]][[comparison]],rowNames=TRUE,colNames=TRUE)
  }
}

saveWorkbook(wb, "./figures/092624_pseudobulkStats_cellPercentageComparisons_wilcoxTest_layer.apoe.condition.xlsx", overwrite = T)

statsList <- c()
for (celltype in celltypes){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation", "condition","sampleid")]
  cellMat$apoe.con <- paste0(cellMat$apoe_geno, ".", cellMat$condition)
  cellMat$layer.con <- paste0(cellMat$CorrectedLayers, ".", cellMat$condition)
  cellMat$layer.apoe <- paste0(cellMat$CorrectedLayers, ".", cellMat$apoe_geno)
  tmp <- reshape2::melt(cellMat)
  colnames(tmp)[grep("variable", colnames(tmp))] <- "Celltypes"
  ### get means from each sample for stats
  #tmp3 <- tmp %>% group_by(sampleid, disease_state, CorrectedLayers, Celltypes) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  tmp2 <- tmp %>% group_by(sampleid, CorrectedLayers, Celltypes, condition, layer.con) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  
  pvals <- tmp %>%
    #group_by(CorrectedLayers) %>%
    wilcox_test(value ~ layer.con) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  statsList[[celltype]][["layer.con.spot.W"]] <- pvals
  
  pvals <- tmp2 %>% ungroup() %>%
    #group_by(CorrectedLayers) %>%
    wilcox_test(Mean_Frequency ~ layer.con) %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  #names(yposition) <- celltypes
  #pvals <- pvals %>% add_xy_position(fun = "max", x = "Celltypes", group = "apoe_geno", dodge = 0.8, step.increase = 0.05)
  #pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position)
  statsList[[celltype]][["layer.con.bulk.W"]] <- pvals
  
  pvals <- tmp2 %>% ungroup() %>%
    #group_by(CorrectedLayers) %>%
    t_test(Mean_Frequency ~ layer.con, p.adjust.method = "none") %>%
    #adjust_pvalue(method = "BH") %>%
    add_significance("p")
  statsList[[celltype]][["layer.con.bulk.T"]] <- pvals
  
}

wb <- createWorkbook()
for (celltype in names(statsList)) {
  for (comparison in names(statsList[[celltype]])){
    name <- paste0(celltype, ".", comparison)
    addWorksheet(wb, name)
    writeData(wb, sheet = name,statsList[[celltype]][[comparison]],rowNames=TRUE,colNames=TRUE)
  }
}

saveWorkbook(wb, "./figures/092624_pseudobulkStats_cellPercentageComparisons_wilcox_TTest_condition.xlsx", overwrite = T)

#### 

dfe33 <- subset(data2, apoe_geno == "E3/3")
dfe34 <- subset(data2, apoe_geno == "E3/4")

test_sign <- rstatix::get_comparisons(dfe33@meta.data, layers.anno)
test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
tmp <- FetchData(dfe33, c("SNCA","layers.anno"))
stats <- compare_means(SNCA ~ layers.anno ,data = tmp, p.adjust.method = "BH")
stats <- stats[stats$p.adj < 0.05,]
stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
test_sign <- lapply(test_sign, function(x){
  if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
    x
  }
})
test_sign <- rlist::list.clean(test_sign)

for (celltype in celltypes){
  goi_vp <- c()
  assay <- paste0("adjExp_", celltype)
  DefaultAssay(dfe33) <- assay
  test_sign <- rstatix::get_comparisons(dfe33@meta.data, layers.anno)
  test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
  tmp <- FetchData(dfe33, c("SNCA","layers.anno"))
  stats <- compare_means(SNCA ~ layers.anno ,data = tmp, p.adjust.method = "BH")
  stats <- stats[stats$p.adj < 0.05,]
  stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
  test_sign <- lapply(test_sign, function(x){
    if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
      x
    }
  })
  test_sign <- rlist::list.clean(test_sign)
  
  #y_max <- 9
  goi_vp[["e33"]] <- VlnPlot(dfe33, 
                    features = "SNCA", 
                    group.by = "layers.anno",
                    #split.by = "LBannotation",
                    cols = lbLayerAnnoCols, 
                    log = F, 
                    same.y.lims = F,
                    assay = assay, 
                    pt.size = 0, 
                    y.max = 2, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/3")+CenterTitle() + ylab("SNCA Expression Level")  + 
    stat_pvalue_manual(stats, y.position = 1)
    stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp[["e34"]] <- VlnPlot(dfe34, 
                             features = "SNCA", 
                             group.by = "CorrectedLayers",
                             split.by = "LBannotation",
                             cols = lbAnnoCols, 
                             log = F, 
                             same.y.lims = F,
                             assay = assay, 
                             pt.size = 0, 
                             #y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/4")+CenterTitle()  + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
                                             axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # + stat_compare_means(label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp <- wrap_plots(goi_vp, ncol = 2, guides = "collect")
  outName <- paste0("./figures/ViolinPlot_SNCA_deconvoluted_",celltype,"_layers_032724.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 6, width = 10)
  
  goi_vp <- c()
 # y_max <- 9.5
  goi_vp[["e33"]] <- VlnPlot(dfe33, 
                    features = "APOE", 
                    group.by = "CorrectedLayers",
                    split.by = "LBannotation",
                    cols = lbAnnoCols, 
                    log = F, 
                    same.y.lims = F,
                    assay = assay, 
                    pt.size = 0,
                   # y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/3")+CenterTitle() + ylab("APOE Expression Level")
  # + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp[["e34"]] <- VlnPlot(dfe34, 
                    features = "APOE", 
                    group.by = "CorrectedLayers",
                    split.by = "LBannotation",
                    cols = lbAnnoCols, 
                    log = F, 
                    same.y.lims = F,
                    assay = assay, 
                    pt.size = 0,
                    # y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/4")+CenterTitle() + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
                                              axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  # + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp <- wrap_plots(goi_vp, ncol = 2, guides = "collect")
  outName <- paste0("./figures/ViolinPlot_APOE_deconvoluted_",celltype,"_layers_032724.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 6, width = 10)
}

for (celltype in celltypes){
  for (anno in names(lbAnnoCols)){
    filename <- paste0("./figures/Boxplot_deconvolutedSpots_", celltype,"_",anno,"_splitAPOE_byLayers.pdf")
    filename <- gsub(" ","",filename)
    saveImages(fileName = filename, bp_anno[[celltype]][[anno]], width = 8, height = 4)
    
  }
}

bp <- c()
for ( i in names(lbAnnoCols)){
  bp[[i]] <- lapply(1:length(celltypes), function(idx) {
    celltype <- celltypes[idx]
    cellMat <- cellPercent[,c(celltype,"CorrectedLayers","apoe_geno", "LBannotation")]
    cellMat <- cellMat[cellMat$LBannotation == i,]
    ggplot(cellMat, aes(x=CorrectedLayers, y=cellMat[,celltype])) + 
      geom_boxplot(aes(fill=CorrectedLayers), alpha=1)  +
      ylim(0,100) +
      scale_fill_manual(values = layerCols2) + 
      ylab(celltype) + xlab("") + theme(axis.text = element_blank(), 
                                        axis.ticks = element_blank(),
                                        legend.position = "none")
  })
  names(bp[[i]]) <- celltypes
}

p <- wrap_plots(bp[[i]], ncol = 1) + theme(axis.text.x = element_text(size = 16))
p 
saveImages("figures/lbSpotsOnly_BoxPlot_percentageCelltypesPerLayer_032624.pdf",p, width = 8, height = 6)


comparisonAPOE <- lapply(levels(data2$CorrectedLayers), function(x) paste0(x,c("_E3/3","_E3/4")))
comparisonDisease <- lapply(levels(data2$CorrectedLayers), function(x) paste0(x,c("_Ctrl","_LBD","Tri")))

library(rstatix)
library(ggpubr)
library(tidyr)
library(tidyverse)
bp <- c()
bp2 <- c()
for (layer in levels(data2$CorrectedLayers)){
  cellMat <- cellPercent[cellPercent$CorrectedLayers == layer,]
  tmp <- reshape2::melt(cellMat)
  colnames(tmp)[grep("variable", colnames(tmp))] <- "Celltypes"
   pvals <- tmp %>%
     group_by(Celltypes) %>%
     wilcox_test(value ~ disease_state) %>%
     adjust_pvalue(method = "BH") %>%
     add_significance("p.adj")
  #names(yposition) <- celltypes
  #pvals <- compare_means(data = tmp, formula = value ~ disease_state, group.by = "Celltypes", p.adjust.method = "BH")
  pvals <- pvals %>% add_xy_position(fun = "max", x = "Celltypes", group = "disease_state", dodge = 0.8, step.increase = 0.05)
  pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position) 
  bp[[layer]] <- ggplot(tmp, aes(x=Celltypes, y=value)) + 
    geom_boxplot(aes(fill=disease_state), outlier.shape = NA)  +
    scale_fill_manual(values = diseaseCols) +
    ylim(0,100) +
    ylab(layer) + 
    xlab("") + 
    theme(#axis.text = element_blank(), 
                                      axis.ticks = element_blank(),
                                      legend.position = "right")  +
  stat_pvalue_manual(data = pvals, label = "p.adj.signif", hide.ns = F, p.adjust.methods="BH", tip.length = 0, bracket.nudge.y = -0.6)
    #stat_compare_means(aes(group=group), label = "p.signif", step.increase = 0.5)
    #stat_pvalue_manual(pvals, label = "p.signif", hide.ns = F, y.position = yposition, step.group.by = "variable")
  pvals <- tmp %>%
    group_by(Celltypes) %>%
    wilcox_test(value ~ apoe_geno) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  #names(yposition) <- celltypes
  pvals <- pvals %>% add_xy_position(fun = "max", x = "Celltypes", group = "apoe_geno", dodge = 0.8, step.increase = 0.05)
  pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position) 
  bp2[[layer]] <- ggplot(tmp, aes(x=Celltypes, y=value)) + 
    geom_boxplot(aes(fill=apoe_geno), outlier.shape = NA)  +
    scale_fill_manual(values = apoeCols) + 
    ylim(0,100) +
    ylab(layer) + 
    xlab("") + 
    theme(#axis.text = element_blank(), 
      axis.ticks = element_blank(),
      legend.position = "right") + 
    stat_pvalue_manual(data = pvals, label = "p.adj.signif", hide.ns = F, tip.length = 0, bracket.nudge.y = -0.6)
  rm(tmp)
  rm(cellMat)
  saveImages(paste0("figures/cellPercentage_BoxPlot_CelltypesPerLayer_diseaseState_",layer,"_042524.png"),bp[[layer]], width = 8, height = 6)
  saveImages(paste0("figures/cellPercentage_BoxPlot_CelltypesPerLayer_apoeGeno_",layer,"_042524.png"),bp2[[layer]], width = 8, height = 6)
}


for (i in levels(data2$CorrectedLayers)[-6]){
  bp[[i]] <- bp[[i]] +  theme(axis.text.x = element_blank())
  bp2[[i]] <- bp2[[i]] +  theme(axis.text.x = element_blank())
}
bp <- wrap_plots(bp, ncol = 1, guides = "collect")
bp2 <- wrap_plots(bp2, ncol = 1, guides = "collect")
saveImages("figures/cellPercentage_BoxPlot_CelltypesPerLayer_diseaseState_032624.pdf",bp, width = 8, height = 6)
saveImages("figures/cellPercentage_BoxPlot_CelltypesPerLayer_apoeGeno_032624.pdf",bp2, width = 8, height = 6)

bp <- c()
bp2 <- c()
pval <- c()
for (celltype in celltypes){
  cellMat <- cellPercent[,c(celltype,"CorrectedLayers", "sampleid", "apoe_geno", "LBannotation", "disease_state")]
  tmp <- reshape2::melt(cellMat)
  colnames(tmp)[grep("variable", colnames(tmp))] <- "Celltypes"
  ### get means from each sample for stats
  tmp3 <- tmp %>% group_by(sampleid, disease_state, CorrectedLayers, Celltypes) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  pvals <- tmp %>%
    group_by(CorrectedLayers) %>%
    wilcox_test(value ~ disease_state) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p")
  #names(yposition) <- celltypes
  #pvals <- compare_means(data = tmp, formula = value ~ disease_state, group.by = "Celltypes", p.adjust.method = "BH")
  pvals <- pvals %>% add_xy_position(fun = "max", x = "CorrectedLayers", group = "disease_state", dodge = 0.8, step.increase = 0.05)
  pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position) 
  pval[[celltype]][["disease_state"]] <- pvals
  bp[[celltype]] <- ggplot(tmp, aes(x=CorrectedLayers, y=value)) + 
    geom_boxplot(aes(fill=disease_state), outlier.shape = NA)  +
    scale_fill_manual(values = diseaseCols) + 
    ylim(0,100) +
    ylab(celltype) + 
    xlab("") + 
    theme(#axis.text = element_blank(), 
      axis.ticks = element_blank(),
      legend.position = "right")+
    stat_pvalue_manual(data = pvals, label = "p.adj.signif", hide.ns = F,
                       tip.length = 0, 
                       #bracket.nudge.y = -0.6, 
                       label.size = 2)
  
  ### get means from each sample for stats
  tmp2 <- tmp %>% group_by(sampleid, apoe_geno, CorrectedLayers, Celltypes) %>% summarise_at(vars(value), list(Mean_Frequency = mean))
  pvals <- tmp2 %>%
    group_by(CorrectedLayers) %>%
    wilcox_test(Mean_Frequency ~ apoe_geno) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  #names(yposition) <- celltypes
  pvals <- pvals %>% add_xy_position(fun = "max", x = "CorrectedLayers", group = "apoe_geno", dodge = 0.8, step.increase = 0.05)
  pvals$y.position <- ifelse(pvals$y.position > 100, 100 - (pvals$y.position -100 ),pvals$y.position) 
  pval[[celltype]][["apoe_geno"]] <- pvals
  bp2[[celltype]] <- ggplot(tmp, aes(x=CorrectedLayers, y=value)) + 
    geom_boxplot(aes(fill=apoe_geno), outlier.shape = NA)  +
    scale_fill_manual(values = apoeCols) + 
    ylim(0,100) +
    ylab(celltype) + 
    xlab("") + 
    theme(#axis.text = element_blank(), 
      axis.ticks = element_blank(),
      legend.position = "right") + 
    stat_pvalue_manual(data = pvals, label = "p.adj.signif", hide.ns = F,
                       tip.length = 0, 
                       #bracket.nudge.y = -0.6, 
                       label.size = 2)
  rm(tmp)
  rm(cellMat)
  #saveImages(paste0("figures/cellPercentage_BoxPlot_CelltypesPerLayer_diseaseState_",celltype,"_042524_layersOnX.png"),bp[[celltype]], width = 8, height = 6)
  #saveImages(paste0("figures/cellPercentage_BoxPlot_CelltypesPerLayer_apoeGeno_",celltype,"_042524_layersOnX.png"),bp2[[celltype]], width = 8, height = 6)
}


for (i in celltypes[-8]){
  bp[[i]] <- bp[[i]] +  theme(axis.text.x = element_blank())
  bp2[[i]] <- bp2[[i]] +  theme(axis.text.x = element_blank())
}
bp_out <- wrap_plots(bp, ncol = 1, guides = "collect")
bp2_out <- wrap_plots(bp2, ncol = 1, guides = "collect")
saveImages("figures/cellPercentage_BoxPlot_CelltypesPerLayer_diseaseState_042524_layersOnX.pdf",bp_out, width = 8, height = 8)
saveImages("figures/cellPercentage_BoxPlot_CelltypesPerLayer_apoeGeno_042524_layersOnX.pdf",bp2_out, width = 8, height = 8)


for (celltype in celltypes){
  for (anno in names(lbAnnoCols)){
    filename <- paste0("./figures/Boxplot_deconvolutedSpots_", celltype,"_",anno,"_splitAPOE_byLayers.pdf")
    filename <- gsub(" ","",filename)
    saveImages(fileName = filename, bp_anno[[celltype]][[anno]], width = 8, height = 4)
  }
}

## plot MIF signals in spatial data
mif <- c( "MIF", "CD74", "CD44", "CXCR4")
sfp <- SpatialFeaturePlot_imagescales(data2, imagescales = imagescales, features = mif)
names(sfp) <- names(imagescales)
saveImages("./figures/Spatial_GliosisMetaFeature.pdf", sfp)
imageAlpha <- 0.3
for (image in 1:length(names(imagescales)) ){
  imageName <- names(imagescales)[image]
  sfp <- SpatialFeaturePlot(data2, 
                     features = mif, 
                     images = imageName, 
                     stroke = 0,
                     pt.size.factor = imagescales[[image]]*2, 
                     alpha = c(0.1,1), 
                     image.alpha = imageAlpha, combine = T)
  saveImages(paste0("./figures/Spatial_MIF_",imageName,".pdf"), sfp, width = 4, height = 5)
}
for (image in 1:length(names(imagescales)) ){
  imageName <- names(imagescales)[[image]]
  sfp <- p[[image]]
  saveImages(paste0("./figures/Spatial_MIF_",imageName,".pdf"), sfp)
}

imageName <- "LBD_5"
imageName <- "LBD_3"

FeaturePlot(lbSpotsOnly, 
            features = c("APOE", "TREM2"), 
            split.by = "LBannotation", 
            pt.size = 1, 
            order = T, 
            blend.threshold = 0.1,
            blend=T, cols = c("grey","red","blue"))

saveImages(paste0("./figures/Spatial_apoe_",imageName,".pdf"), p)


SpatialFeaturePlot(data2, 
                   features = mif, 
                   images = imageName, 
                   pt.size.factor = imagescales[[image]], 
                   alpha = c(0.1,1), 
                   image.alpha = imageAlpha)
SpatialPlot(data2, images = "LBD_3", facet.highlight = T, 
            cells.highlight = CellsByIdentities(data2, idents = c("LB(+)","LB Surround")),
            #features = "APOE", 
            #cells.highlight = c("LB(+)", "LB Surround")
            ) + ggplot2::theme_minimal()


### spatial plots of genes
goi <- c("BIN1", "TMEM175", "GBA","PRKN")
imageAlpha <- 0.3
for (image in 1:length(names(imagescales)) ){
  imageName <- names(imagescales)[image]
  sfp <- lapply(goi, function(gene) { SpatialFeaturePlot(data2, 
                            features = gene, 
                            images = imageName, 
                            stroke = 0,
                            pt.size.factor = imagescales[[image]]*2, 
                            alpha = c(0.1,1), 
                            image.alpha = imageAlpha, combine = T)
  #saveImages(paste0("./figures/Spatial_",gene,"_",imageName,".pdf"), sfp, width = 4, height = 5)
  })
}
for (image in 1:length(names(imagescales)) ){
  imageName <- names(imagescales)[[image]]
  sfp <- p[[image]]
  #saveImages(paste0("./figures/Spatial_",gene,"_",imageName,".pdf"), sfp)
}

for (gene in goi){
  sfp <- SpatialFeaturePlot_imagescales(data2, features = gene, imagescales = imagescales, 
                                        stroke = 0, imageAlpha = 0.1)
  p <- wrap_plots(sfp, ncol = 5, guides = 'collect')
  filename <- paste0("SpatialPlot_",gene,".pdf")
  saveImages(filename, p, width = 10, height = 6)
}


### violin plots
goi_vp <- VlnPlot(data2, features = goi, group.by = "layersAPOE", cols = layerColsAPOE, pt.size = 0, )
deconAss <- names(data2@assays)[-c(1,2)]
for (gene in goi) {
  yMax[[gene]] <- max(unlist(lapply(deconAss, function(assay) FetchData(data2, vars = gene, layer = 'data', assay=assay))))
}

for (gene in goi) {
  vlList <- list()
  command <- paste0("subset(data2, ", gene,"  > 0)")
  df <- eval(parse(text = command))
  for (celltype in celltypes){
    try(vlList[[celltype]] <- VlnPlot(df, assay = paste0("adjExp_",celltype), 
                                  features = gene, 
                                  group.by = "CorrectedLayers", 
                                  cols = layerCols2, 
                                  pt.size = 0, y.max = yMax[[gene]]
    ))
  }
  try(vlList <- lapply(celltypes, function(i){
    vlList[[i]] + FontSize(x.title = 0, main = 0, x.text = 0, y.title = 8, y.text = 4) + ylab(i)
  }))
  try(p <- wrap_plots(vlList, ncol=1, guides = 'collect', widths = 5, heights = 1) + FontSize(x.text = 8))
  pdf(paste0("./figures/ViolinPlot_layers_byCelltype_compressed_",gene,".pdf"),width = 8, height = 10)
  p
  dev.off()
}

for (gene in goi) {
  command <- paste0("subset(data2, ", gene,"  > 0)")
  df <- eval(parse(text = command))
  vp <- VlnPlot(df, assay = "Spatial", 
          features = gene, 
          group.by = "CorrectedLayers", 
          cols = layerCols2, 
          pt.size = 0.1
  )
  filename <- paste0("./figures/VlnPlot_",gene,"_spotLevel_CorrectedLayers_dropZeros.pdf")
  saveImages(filename,vp, dpi = 300, units = "in", height = 6, width = 8)
}
for (gene in goi) {
  command <- paste0("subset(data2, ", gene,"  > 0)")
  df <- eval(parse(text = command))
  vp <- VlnPlot(df, assay = "Spatial", 
                features = gene, 
                group.by = "layersAPOE", 
                cols = layerColsAPOE, 
                pt.size = 0.1
  )
  filename <- paste0("./figures/VlnPlot_",gene,"_spotLevel_layersAPOE_dropZeros.pdf")
  saveImages(filename,vp, dpi = 300, units = "in", height = 6, width = 8)
}


### cellchat binary plots if possible
#SELENOP-LRP8, APOE-LRP8 and SELENOP-APOE

genePairs <- c("SELENOP_LRP8", "APOE_LRP8", "SELENOP_APOE")

for (genePair in genePairs) {
  genePair <- strsplit(genePair, split = "_")[[1]]
  
  sfp <- FeaturePlot(data2, 
              features = genePair, 
              split.by = "LBannotation", 
              pt.size = 3, 
              #blend.threshold = 0.1,
              blend=T, 
              cols = c("grey","red","blue"), 
              order = T, 
              coord.fixed = T)

  #sfp <- sfp[imageOrder]
  #sfp <- wrap_plots(sfp, ncol = 5, guides = 'collect')
  filename <- paste0("./figures/061424_BinaryFeaturePlot_",paste0(genePair, collapse = "_"),".pdf")
  saveImages(fileName = filename, plotToSave = sfp)  

}

### modify cellchatAPOE to match expected
for (name in c("coordinates", "scale.factors")){
  cellchatAPOE2@images[[name]] <- rbind(cellchatAPOE2@images$APOE33[[name]],cellchatAPOE2@images$APOE34[[name]])
}

imageOrder <- c("LBD_1","LBD_7","LBD_6","LBD_8","LBD_3",
                "LBD_5","LBD_9","LBD_11","LBD_12","LBD_4")

for (slice in imageOrder){
  for (genePair in genePairs) {
    genePair <- strsplit(genePair, split = "-")[[1]]
    gg <- spatialFeaturePlot_groups(cellchatAPOE, 
                                    #signaling = pathways.show, 
                                    ncol = 1,
                                    #n.colors = 25,
                                    pairLR.use = genePair,
                                    slice.use = slice, point.size = 1.3, 
                                    do.binary = TRUE, cutoff = 0.05, show.legend.combined = F,
                                    enriched.only = F, color.heatmap = "Reds", direction = 1)
    command <- paste0("gg + ", adjustment[[slice]])
    gg <- eval(parse(text = command))
    fileName <- paste0("./figuresBySample/BinarySpatialFeaturePlot_cellChat_LR_receptorComplex_",
                       pathways.show,"_",slice,"_",pair,"_041024.pdf")
    try(saveImages(fileName, gg, width = 8, height = 6))
  }
}

for (genePair in genePairs){
  genes <- strsplit(genePair, split = "_")[[1]]
  allSpots <- colnames(data2@assays$SCT@data)
  
  spotsL <- data2@assays$SCT@data[genes[1],]
  cutoff <- quantile(spotsL[spotsL>0], 0.25)
  spotsL <- names(spotsL)[spotsL > cutoff]
  
  spotsR <- data2@assays$SCT@data[genes[2],]
  cutoff <- quantile(spotsR[spotsR>0], 0.25)
  spotsR <- names(spotsR)[spotsR > cutoff]
  
  spotsBoth <- spotsL[spotsL %in% spotsR] 
  
  idlist <- allSpots
  allSpots <- ifelse(allSpots %in% spotsBoth, "Both",
                     ifelse(allSpots %in% spotsL, genes[1], 
                     ifelse(allSpots %in% spotsR, genes[2], "Neither")))
  tmp <- as.data.frame(x = allSpots, idlist)
  data2 <- AddMetaData(data2, metadata = tmp, col.name = "Gene Pair")
  colorList <- list("firebrick2",
                    "forestgreen",
                    "dodgerblue2",
                    "grey")
  names(colorList) <- c(genes[[1]], genes[[2]], "Both", "Neither")
  ggs <- SpatialDimPlot_imagescales(data2, imagescales = imagescales, 
                             stroke = 0, 
                             group.by="Gene Pair", 
                             cols=colorList)
  ggs <- ggs[imageOrder]
  gg <- wrap_plots(ggs, ncol=5, guides = 'collect')
  rm(ggs)
  fileName <- paste0("./figures/BinarySpatialFeaturePlot_genePairs_",
                     genePair,"_061724.pdf")
  saveImages(fileName, gg, width = 15, height = 8)
  rm(gg, tmp, idlist, allSpots, spotsL, spotsR, genes)
}


genes <- c("GPX4", "SELENOP", "SELENOS", "SELENOK", "LRP8")
for ( gene in genes ){
    test_sign <- rstatix::get_comparisons(data2@meta.data, variable = CorrectedLayers, )
  #test_sign <- test_sign[c(1,22,39,52,61,66)]
  y_max <- 8
  
  goi_vp <- VlnPlot(data2, 
                    features = gene, 
                    group.by = "CorrectedLayers",
                    #split.by = "apoe_geno",
                    cols = layerCols2, 
                    log = F, 
                    same.y.lims = F,
                    assay = "Spatial", 
                    pt.size = 0,
                    y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + stat_compare_means(comparisons = test_sign, 
                         label = "p.signif", hide.ns = T, 
                         step.increase = .05, vjust = .8, 
                         p.adjust.methods="BH")
  outName <- paste0("./figures/061724_ViolinPlot_",gene,"_spotLevel_layers_pvalues.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)
  
  goi_vp <- VlnPlot(subset(data2, cells= spots), 
                    features = gene, 
                    group.by = "CorrectedLayers",
                    #split.by = "apoe_geno",
                    cols = layerCols2, 
                    log = F, 
                    same.y.lims = F,
                    assay = "Spatial", 
                    pt.size = 0,
                    y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
  ) + stat_compare_means(comparisons = test_sign, 
                         label = "p.signif", hide.ns = T, 
                         step.increase = .05, vjust = .8, 
                         p.adjust.methods="BH")
  outName <- paste0("./figures/061724_ViolinPlot_",gene,"_spotLevel_layers_pvalues_dropZeros.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 4, width = 6)
}


dp <- DimPlot(data2, group.by = "CorrectedLayers", cols = layerCols2, shuffle = T, raster = F) + 
  xlab("UMAP_1") + ylab("UMAP_2")
saveImages("./figures/082224_UMAP_CorrectedLayers.png", dp, width = 12, height = 10, dpi = 300, units = "in")


### violinplots snca and apoe
library(ggpubr)
gene <- "SNCA"
test_sign <- rstatix::get_comparisons(data2@meta.data, variable = layers.anno)
test_sign <- lapply(test_sign, function(pair){
  p1 <- gsub("_.*","",pair[1])
  p2 <- gsub("_.*","",pair[2])
  if( p1 == p2 ){
    pair
  } else { NULL }
})
test_sign <- rlist::list.clean(test_sign)
test_sign <- test_sign[-c(1,2,4,7,8,10,13,19,20,25,26)]
data2$LBannotation <- factor(data2$LBannotation, levels = c("LB(+)",
                                                            "LB Surround",
                                                            "LB(-)",
                                                            "CLB(-)"))
data2$layers.anno <- factor(data2$layers.anno, levels=names(lbLayerAnnoCols))

goi_vp <- VlnPlot(data2,  
                  features = gene, 
                  group.by = "layers.anno",
                  #split.by = "apoe_geno",
                  cols = lbLayerAnnoCols, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = 14,# add the y-axis maximum value - otherwise p-value hidden
) + NoLegend() + stat_compare_means(comparisons = test_sign, 
                       label = "p.signif", hide.ns = T, 
                       step.increase = .08, vjust = .8, 
                       p.adjust.methods="BH"); goi_vp
outName <- paste0("./figures/090524_ViolinPlot_",gene,"_spotLevel_layers_anno_pvalues.pdf")
saveImages(outName, goi_vp, dpi = 300, height = 6, width = 8)

gene <- "APOE"
test_sign <- rstatix::get_comparisons(data2@meta.data, variable = layers.anno)
test_sign <- lapply(test_sign, function(pair){
  p1 <- gsub("_.*","",pair[1])
  p2 <- gsub("_.*","",pair[2])
  if( p1 == p2 ){
    pair
  } else { NULL }
})
test_sign <- rlist::list.clean(test_sign)
test_sign <- test_sign[-c(1:5,7,8,10,13,14,16,19,25)]
goi_vp <- VlnPlot(data2,  
                  features = gene, 
                  group.by = "layers.anno",
                  #split.by = "apoe_geno",
                  cols = lbLayerAnnoCols, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = 14,# add the y-axis maximum value - otherwise p-value hidden
) + NoLegend() + stat_compare_means(comparisons = test_sign, 
                                    label = "p.signif", hide.ns = T, 
                                    step.increase = .08, vjust = .8, 
                                    p.adjust.methods="BH"); goi_vp
outName <- paste0("./figures/090524_ViolinPlot_",gene,"_spotLevel_layers_anno_pvalues.pdf")
saveImages(outName, goi_vp, dpi = 300, height = 6, width = 8)

##combined e3 and e4 both by layers and all together
gene <- "SNCA"
test_sign <- rstatix::get_comparisons(data2@meta.data, variable = layers.anno)
test_sign <- lapply(test_sign, function(pair){
  p1 <- gsub("_.*","",pair[1])
  p2 <- gsub("_.*","",pair[2])
  if( p1 == p2 ){
    pair
  } else { NULL }
})
test_sign <- rlist::list.clean(test_sign)
test_sign <- test_sign[-c(1,2,4,7,8,10,13,19,20,25,26)]
data2$LBannotation <- factor(data2$LBannotation, levels = c("LB(+)",
                                                            "LB Surround",
                                                            "LB(-)",
                                                            "CLB(-)"))
data2$layers.anno <- factor(data2$layers.anno, levels=names(lbLayerAnnoCols))

goi_vp <- VlnPlot(data2,  
                  features = gene, 
                  group.by = "layers.anno",
                  #split.by = "apoe_geno",
                  cols = lbLayerAnnoCols, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = 14,# add the y-axis maximum value - otherwise p-value hidden
) + NoLegend() + stat_compare_means(comparisons = test_sign, 
                                    label = "p.signif", hide.ns = T, 
                                    step.increase = .08, vjust = .8, 
                                    p.adjust.methods="BH"); goi_vp
outName <- paste0("./figures/090524_ViolinPlot_",gene,"_spotLevel_layers_anno_pvalues.pdf")
saveImages(outName, goi_vp, dpi = 300, height = 6, width = 8)

gene <- "APOE"
test_sign <- rstatix::get_comparisons(data2@meta.data, variable = layers.anno)
test_sign <- lapply(test_sign, function(pair){
  p1 <- gsub("_.*","",pair[1])
  p2 <- gsub("_.*","",pair[2])
  if( p1 == p2 ){
    pair
  } else { NULL }
})
test_sign <- rlist::list.clean(test_sign)
test_sign <- test_sign[-c(1:5,7,8,10,13,14,16,19,25)]
goi_vp <- VlnPlot(data2,  
                  features = gene, 
                  group.by = "layers.anno",
                  #split.by = "apoe_geno",
                  cols = lbLayerAnnoCols, 
                  log = F, 
                  same.y.lims = F,
                  assay = "Spatial", 
                  pt.size = 0,
                  y.max = 14,# add the y-axis maximum value - otherwise p-value hidden
) + NoLegend() + stat_compare_means(comparisons = test_sign, 
                                    label = "p.signif", hide.ns = T, 
                                    step.increase = .08, vjust = .8, 
                                    p.adjust.methods="BH"); goi_vp
outName <- paste0("./figures/090524_ViolinPlot_",gene,"_spotLevel_layers_anno_pvalues.pdf")
saveImages(outName, goi_vp, dpi = 300, height = 6, width = 8)


### boxplot of LB counts per layer/ apoe
tmp <- reshape2::melt(table(data2$layers.anno.APOE, data2$sampleid))
genoSamp <- reshape2::melt(table(data2$apoe_geno, data2$sampleid))
genoSamp <- genoSamp[genoSamp$value != 0,]
genoSamp$sample.geno <- paste0(genoSamp$Var2, ".", genoSamp$Var1)
tmp <- tmp[grep("LB(+)",tmp$Var1, fixed = T),]
tmp$geno <- str_split_i(tmp$Var1, pattern = "_", i = 3)
tmp$sample.geno <- paste0(tmp$Var2, ".", tmp$geno)
tmp <- tmp[tmp$sample.geno %in% genoSamp$sample.geno,]
tmp <- tmp[tmp$Var2 %in% lbd_images,]
tmp$Layer <- str_split_i(tmp$Var1, "_", i=1)
tmp$Layer.Geno <- paste0(tmp$Layer, ".", tmp$geno)

test_sign <- rstatix::get_comparisons(tmp, Layer.Geno)
test_sign <- lapply(test_sign, function(pair){
  p1 <- gsub("\\..*","",pair[1], perl = T)
  p2 <- gsub("\\..*","",pair[2], perl=T)
  if( p1 == p2 ){
    pair
  } else { NULL }
})
test_sign <- rlist::list.clean(test_sign)
stats <- compare_means(value ~ Layer.Geno ,data = tmp, method = "t.test")
stats <- stats[gsub("\\..*","",stats$group1, perl=T) == gsub("\\..*","",stats$group2, perl=T),]


bp <- tmp %>% ggplot(aes(x=Layer.Geno, y=value, fill=geno)) + 
  geom_boxplot() +
  scale_fill_manual(values = unname(apoeColors)) +
  geom_jitter()+ rotate_x_text() +
  stat_compare_means(comparisons =test_sign, 
                     #method = "t.test",
                     label = "p.signif", hide.ns = F, 
                     step.increase = .08, vjust = .8, 
                     p.adjust.methods="BH"); bp

outName <- paste0("./figures/090524_Boxplot_LB+_spotLevel_layers_anno_splitAPOE_pvalues.pdf")
saveImages(outName, bp, dpi = 300, height = 6, width = 8)

tmp <- reshape2::melt(table(data2$layers.anno, data2$sampleid))
tmp <- tmp[grep("LB(+)",tmp$Var1, fixed = T),]
tmp$Var1 <- factor(gsub("_LB(+)","",tmp$Var1, fixed = T), levels = names(layerCols2))

tmp <- tmp[tmp$Var2 %in% lbd_images,]
test_sign <- rstatix::get_comparisons(tmp, Var1)
test_sign <- test_sign[-c(4,5,6,8,11)]
bp <- tmp %>% ggplot(aes(x=Var1, y=value, fill=Var1)) + 
  geom_boxplot(outliers = F) + 
  scale_fill_manual(values = layerCols2) +
  geom_jitter()+ rotate_x_text() +
  stat_compare_means(comparisons= test_sign,
                     label = "p.signif", 
                     #method = "t.test",
                     hide.ns = T, 
                     step.increase = .08, vjust = .8, 
                     p.adjust.methods="BH"); bp
outName <- paste0("./figures/090524_Boxplot_LB+_spotLevel_layers_anno_pvalues.pdf")
saveImages(outName, bp, dpi = 300, height = 6, width = 8)

### celltype snca/apoe
dfe33 <- subset(data2, apoe_geno == "E3/3")
dfe34 <- subset(data2, apoe_geno == "E3/4")

for (celltype in celltypes){
  goi_vp <- c()
  assay <- paste0("adjExp_", celltype)
  DefaultAssay(dfe33) <- assay
  DefaultAssay(data2) <- assay
  
  test_sign <- rstatix::get_comparisons(dfe33@meta.data, layers.anno)
  test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
  tmp <- FetchData(dfe33, c("SNCA","layers.anno"))
  stats <- compare_means(SNCA ~ layers.anno ,data = tmp, p.adjust.method = "BH")
  stats <- stats[stats$p.adj < 0.05,]
  stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
  test_sign1 <- lapply(test_sign, function(x){
    if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
      x
    }
  })
  test_sign1 <- rlist::list.clean(test_sign1)
  #ymax <- max(tmp$SNCA)*(nrow(stats)/4.5)
  
  DefaultAssay(dfe34) <- assay
  test_sign <- rstatix::get_comparisons(dfe34@meta.data, layers.anno)
  test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
  tmp <- FetchData(dfe34, c("SNCA","layers.anno"))
  stats <- compare_means(SNCA ~ layers.anno ,data = tmp, p.adjust.method = "BH")
  stats <- stats[stats$p.adj < 0.05,]
  stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
  test_sign2 <- lapply(test_sign, function(x){
    if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
      x
    }
  })
  test_sign2 <- rlist::list.clean(test_sign2)
  
  ymax <- max(FetchData(data2, "SNCA"))*(max(length(test_sign1), length(test_sign2))/6.5)
  
  #y_max <- 9
  goi_vp[["e33"]] <- VlnPlot(dfe33, 
                             features = "SNCA", 
                             group.by = "layers.anno",
                             #split.by = "LBannotation",
                             cols = lbLayerAnnoCols, 
                             log = F, 
                             same.y.lims = F,
                             assay = assay, add.noise = F,
                             pt.size = 0, 
                             y.max = ymax, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/3")+CenterTitle() + ylab("SNCA Expression Level")  + NoLegend() +rotate_x_text() + 
    stat_compare_means(comparisons = test_sign1, 
                         label = "p.signif", hide.ns = T, 
                         step.increase = .05, vjust = .8, 
                         p.adjust.methods="BH"); goi_vp[["e33"]]


  #ymax <- max(tmp$SNCA)*(nrow(stats)/4.5)
  goi_vp[["e34"]] <- VlnPlot(dfe34, 
                             features = "SNCA", 
                             group.by = "layers.anno",
                             #split.by = "LBannotation",
                             cols = lbLayerAnnoCols, 
                             log = F,
                             same.y.lims = F,
                             assay = assay, add.noise = F,
                             pt.size = 0, 
                             y.max = ymax, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/4")+ CenterTitle()  +NoLegend() + rotate_x_text() + 
    stat_compare_means(comparisons = test_sign2, 
                       label = "p.signif", hide.ns = T, 
                       step.increase = .05, vjust = .8, 
                       p.adjust.methods="BH") +
  theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
                                            axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # + stat_compare_means(label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp <- wrap_plots(goi_vp, ncol = 2, guides = "collect")
  outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA_deconvoluted_",celltype,"_layers.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 6, width = 10)
  
  goi_vp <- c()

  test_sign <- rstatix::get_comparisons(dfe33@meta.data, layers.anno)
  test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
  tmp <- FetchData(dfe33, c("APOE","layers.anno"))
  stats <- compare_means(APOE ~ layers.anno ,data = tmp, p.adjust.method = "BH")
  stats <- stats[stats$p.adj < 0.05,]
  stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
  test_sign1 <- lapply(test_sign, function(x){
    if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
      x
    }
  })
  test_sign1 <- rlist::list.clean(test_sign1)
  
  #ymax <- max(tmp$APOE)*(nrow(stats)/4.5)
  
  test_sign <- rstatix::get_comparisons(dfe34@meta.data, layers.anno)
  test_sign <- test_sign[unlist(lapply(test_sign, function(x) strsplit(x[1], split = "_")[[1]][1] == strsplit(x[2], split = "_")[[1]][1]))]
  tmp <- FetchData(dfe34, c("APOE","layers.anno"))
  stats <- compare_means(APOE ~ layers.anno ,data = tmp, p.adjust.method = "BH")
  stats <- stats[stats$p.adj < 0.05,]
  stats <- stats[gsub("_.*", "", stats$group1, perl = T) == gsub("_.*", "", stats$group2, perl = T) ,]
  test_sign2 <- lapply(test_sign, function(x){
    if(nrow(stats[stats$group1 == x[[1]] & stats$group2 == x[[2]],]) != 0){
      x
    }
  })
  test_sign2 <- rlist::list.clean(test_sign2)
  #ymax <- max(tmp$APOE)*(nrow(stats)/4.5)
  ymax <- max(FetchData(data2, "APOE"))*(max(length(test_sign1), length(test_sign2))/6.5)
  
  goi_vp[["e33"]] <- VlnPlot(dfe33, 
                             features = "APOE", 
                             group.by = "layers.anno",
                             #split.by = "LBannotation",
                             cols = lbLayerAnnoCols, 
                             log = F, 
                             same.y.lims = F,
                             assay = assay, add.noise = F,
                             pt.size = 0,
                             y.max = ymax, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/3")+CenterTitle() + ylab("APOE Expression Level") + NoLegend() + rotate_x_text() + 
    stat_compare_means(comparisons = test_sign1, 
                       label = "p.signif", hide.ns = T, 
                       step.increase = .05, vjust = .8, 
                       p.adjust.methods="BH")
  
  
  goi_vp[["e34"]] <- VlnPlot(dfe34, 
                             features = "APOE", 
                             group.by = "layers.anno",
                             #split.by = "LBannotation",
                             cols = lbLayerAnnoCols, 
                             log = F, 
                             same.y.lims = F,
                             assay = assay, add.noise = F,
                             pt.size = 0,
                             y.max = ymax, # add the y-axis maximum value - otherwise p-value hidden
  ) + ggtitle("E3/4")+CenterTitle() + NoLegend() + rotate_x_text() + 
    stat_compare_means(comparisons = test_sign2, 
                       label = "p.signif", hide.ns = T, 
                       step.increase = .05, vjust = .8, 
                       p.adjust.methods="BH")+
    theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # + stat_compare_means(comparisons = test_sign, label = "p.signif", hide.ns = T, step.increase = .05,p.adjust.methods="BH")
  goi_vp <- wrap_plots(goi_vp, ncol = 2, guides = "collect"); goi_vp
  outName <- paste0("./figures/",timeStamp,"_ViolinPlot_APOE_deconvoluted_",celltype,"_layers.pdf")
  saveImages(outName,goi_vp, dpi = 300, units = "in",  height = 6, width = 10)
}

### violin plots for SNCA split
data2$layer.snca <- paste0(data2$CorrectedLayers, ".", data2$SNCAexpr)
include <- grep("NA", unique(data2$layer.snca), invert = T, value = T)
Idents(data2) <- "layer.snca"
vp <- VlnPlot(data2, features = "SNCA",
        group.by = "layer.snca", 
        pt.size = 0, 
        idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_SNCA_AllSpots_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(data2, features = "APOE",
              group.by = "layer.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_APOE_AllSpots_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

lbpos <- subset(data2, subset= LBannotation == "LB(+)")
include <- grep("NA", unique(lbpos$layer.snca), invert = T, value = T)
vp <- VlnPlot(lbpos, features = "SNCA",
              group.by = "layer.snca", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_SNCA_LBpos_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbpos, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_APOE_LBpos_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

lbsur <- subset(data2, subset= LBannotation == "LB Surround")
include <- grep("NA", unique(lbsur$layer.snca), invert = T, value = T)
vp <- VlnPlot(lbsur, features = "SNCA",
              group.by = "layer.snca", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_SNCA_LBsur_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbsur, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_APOE_LBsur_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

lbneg <- subset(data2, subset= LBannotation == "LB(-)")
include <- grep("NA", unique(lbneg$layer.snca), invert = T, value = T)
vp <- VlnPlot(lbneg, features = "SNCA",
              group.by = "layer.snca", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_SNCA_LBneg_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbneg, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA50_APOE_LBneg_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)


### snca 25
data2$layer.snca25 <- paste0(data2$CorrectedLayers, ".", data2$SNCAexpr25)
include <- grep("NA", unique(data2$layer.snca25), invert = T, value = T)
Idents(data2) <- "layer.snca25"
vp <- VlnPlot(data2, features = "SNCA", 
        group.by = "layer.snca25", 
        pt.size = 0, 
        idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_SNCA_AllSpots_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(data2, features = "APOE",
              group.by = "layer.snca25", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_APOE_AllSpots_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

#lbpos <- subset(data2, subset= LBannotation == "LB(+)")
Idents(lbpos) <- "layer.snca25"
include <- grep("NA", unique(lbpos$layer.snca25), invert = T, value = T)
vp <- VlnPlot(lbpos, features = "SNCA",
              group.by = "layer.snca25", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_SNCA_LBpos_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbpos, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_APOE_LBpos_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

#lbsur <- subset(data2, subset= LBannotation == "LB Surround")
Idents(lbsur) <- "layer.snca25"
include <- grep("NA", unique(lbsur$layer.snca25), invert = T, value = T)
vp <- VlnPlot(lbsur, features = "SNCA",
              group.by = "layer.snca25", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_SNCA_LBsur_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbsur, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_APOE_LBsur_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

#lbneg <- subset(data2, subset= LBannotation == "LB(-)")
Idents(lbneg) <- "layer.snca25"
include <- grep("NA", unique(lbneg$layer.snca25), invert = T, value = T)
vp <- VlnPlot(lbneg, features = "SNCA",
              group.by = "layer.snca25", 
              #split.by = "LBannotation", 
              #split.plot = F,
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_SNCA_LBneg_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

vp <- VlnPlot(lbneg, features = "APOE",
              group.by = "layer.anno.snca", 
              pt.size = 0, 
              idents = include) + NoLegend()

outName <- paste0("./figures/",timeStamp,"_ViolinPlot_SNCA25_APOE_LBneg_byLayers.pdf")
saveImages(outName,vp, dpi = 300, units = "in",  height = 8, width = 8)

### spatial plots for marker genes:
markerGenes <- c("RELN","SPARC","CXCL14","FOS","ID3","AEBP1", 
                 "LAMP5", "MGP", "PCDH8", "HPCAL1", "HOPX", "TESPA1",
                 "NEFM", "CARTPT", "COL5A2", "HAPLN4", "SCN1B", "SNCG",
                 "CLSTN2", "STMN2", "SYT1", "PCP4", "SMYD2", "EFHD2", 
                 "NPY", "SEMA3E", "DIRAS2", "CPLX3", "HS3ST4", "B3GALT2",
                 "MBP", "PLP1", "SPP1", "MOBP", "TF","CLDN11")
imagescalesBig <- lapply(imagescales, function(x) x*10000)
adjustedImagescales <- imagescalesBig
adjustedImagescales$LBD_3 <- (adjustedImagescales$LBD_3)/5.6
adjustedImagescales$LBD_1 <- (adjustedImagescales$LBD_1)/3.5
adjustedImagescales$LBD_6 <- (adjustedImagescales$LBD_6)/3.5
adjustedImagescales$LBD_5 <- (adjustedImagescales$LBD_5)/4
adjustedImagescales$LBD_4 <- (adjustedImagescales$LBD_4)/3.9
adjustedImagescales$LBD_7 <- (adjustedImagescales$LBD_7)/2.5
adjustedImagescales$LBD_12 <- (adjustedImagescales$LBD_12)/1.05
adjustedImagescales$LBD_11 <- (adjustedImagescales$LBD_11)/1.05

DefaultAssay(data2) <- "SCT"

for (gene in markerGenes){
  sfp <- SpatialFeaturePlot_imagescales(data2, 
                                        features = gene, 
                                        imagescales = adjustedImagescales, 
                                        imageAlpha = 0,
                                        alpha = 1)
  sfp <- wrap_plots(sfp[imageOrder], ncol = 5, guides = "collect"); sfp
  ggsave(paste0("./figures/",timeStamp,"_MarkerGenes_spatialDimPlot_",gene,".pdf"), sfp, width = 12, height = 6, dpi = 300, units = "in")
}

