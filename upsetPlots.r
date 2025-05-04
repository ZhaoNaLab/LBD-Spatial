### Upset plots 
### E3 vs E4 across each layers/LB+, LBsur, LB- spots. 

library(UpSetR)
library(colorRamps)
library(tidyverse)
library(venn)
library(openxlsx)
setwd("./")
timeStamp <- format(Sys.time(), "%m%d%y")

outdir <- "./figures/"
  
layerCols <- c( 
  "L1" ="#8D405C",
  "L23"= "#E7BDE1",
  "L4"= "#CF8CA4",
  "L5"= "#9F6E80",
  "L6"= "#CDADB9",
  "WM"= "#67A9D8")

apoeColors <- c("E3"="#2367AC", 
                "E4"="#B21F2C")

layerAPOEcols <- rep(apoeColors,6)
names(layerAPOEcols) <- paste0(sort(rep(names(layerCols), 2)),"_",names(layerAPOEcols))

### read in common markers
lbpos <- readRDS("lbSpotsOnly_commonmarkers.rds")
lbsur <- readRDS("lbSurroundOnly_commonmarkers.rds")
lbneg <- readRDS("lbNeg_commonmarkers.rds")
lbcon <- readRDS("lbConNeg_commonmarkers.rds")

### create named lists and adjust for E3 and E4 labels
for (name in c("lbpos", "lbsur","lbneg","lbcon")){
  input <- get(name)
  outName <- paste0(name,"_split")
  split <- c()
  for (layer in names(input)) {
    i <- gsub("ayer","", layer)
    df <- input[[layer]]
    E3 <- df[df$avg_log2FC.spot > 0.25,]
    E4 <- df[df$avg_log2FC.spot < -0.25,]
    split[[paste0(i,"_E3")]] <- E3
    split[[paste0(i,"_E4")]] <- E4
  }
  assign(outName, split)
  rm(split, input, outName, layer, E3, E4)
}

for (name in  ls(pattern = "split$")) {
  df <- get(name)
  genesOnly <- lapply(df, function(x) x$genes)
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  p <- upset(fromList(genesOnly), 
        nintersects = NA, mainbar.y.label = paste0(gsub("_split","",name), " Intersection Size"),
        sets.bar.color = layerAPOEcols[names(genesOnly)], 
        sets = rev(names(genesOnly)), 
        keep.order = T,
        order.by = "freq", 
        point.size = 3, 
        main.bar.color = "#0c6b58", 
        text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
        #boxplot.summary = c("L1_E3","L2_E3")
        )
  numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
  numLabels <- length(p$labels)
  fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_",name,".pdf")
  pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
  print(p)
  dev.off()
  rm(p, fileName, genesOnly, df)
}

### get the venn gene matrix:
wb <- createWorkbook()
for (name in  ls(pattern = "split$")) {
  df <- get(name)
  genesOnly <- lapply(df, function(x) x$genes)
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  elements <- unique(unlist(genesOnly))
  data <- unlist(lapply(genesOnly, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(genesOnly)
  rownames(data) <- elements
  addWorksheet(wb, name)
  writeData(wb = wb, sheet = name, x = data, , rowNames = T)
}
fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_byAnnotation_genes.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

mergedGenes <- c()
for (name in  ls(pattern = "split$")) {
  df <- get(name)
  names(df) <- paste0(gsub("_split","",name),"_",names(df))
  genesOnly <- lapply(df, function(x) x$genes)
  mergedGenes <- c(mergedGenes, genesOnly)
}

### for each layer
for (layer in names(layerCols)) {
  idx <- grep(layer, names(mergedGenes))
  genesOnly <- mergedGenes[idx]
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  idx2 <- paste0('^.*_', layer)
  p <- upset(fromList(genesOnly), 
             nintersects = NA, 
             mainbar.y.label = paste0(layer, " Intersection Size"),
             sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
             sets = rev(names(genesOnly)), 
             keep.order = T,
             order.by = "freq", 
             point.size = 3, 
             main.bar.color = "#0c6b58", 
             text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
             #boxplot.summary = c("L1_E3","L2_E3")
  )
  numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
  numLabels <- length(p$labels)
  fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_",layer,".pdf")
  pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
  print(p)
  dev.off()
}

### get the venn gene matrix:
wb <- createWorkbook()
for (layer in names(layerCols)) {
  idx <- grep(layer, names(mergedGenes))
  genesOnly <- mergedGenes[idx]
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  elements <- unique(unlist(genesOnly))
  data <- unlist(lapply(genesOnly, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(genesOnly)
  rownames(data) <- elements
  addWorksheet(wb, layer)
  writeData(wb = wb, sheet = layer, x = data, rowNames = T)
}
fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_byLayer_genes.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)



### deconvoluted ex 
file <- "DEGs/DESeq2_common_deconvolutedCelltypes_byLBAnnotation_APOE_Layers_marker_genes_LFC-consistent.xlsx"
sheets <- getSheetNames("DEGs/DESeq2_common_deconvolutedCelltypes_byLBAnnotation_APOE_Layers_marker_genes_LFC-consistent.xlsx")
selectedSheets <- grep("_Ex_",sheets, value = T)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(selectedSheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- selectedSheets
for (name in names(df)){
  if(is.null(df[[name]])){
    df[[name]] <- NULL
  }
}
lbpos <- df[grep("LB(+)", names(df), fixed = T, value=T)]
names(lbpos) <- gsub("adjExp_Ex_|_LB\\(\\+\\)", "", names(lbpos))
lbsur <- df[grep("LB Surround", names(df), fixed = T, value=T)]
names(lbsur) <- gsub("adjExp_Ex_|_LB\ Surround", "", names(lbsur))
lbneg <- df[grep("_LB(-)", names(df), fixed = T, value=T)]
names(lbneg) <- gsub("adjExp_Ex_|_LB\\(-\\)", "", names(lbneg))
lbcon <- df[grep("CLB(-)", names(df), fixed = T, value=T)]
names(lbcon) <- gsub("adjExp_Ex_|_CLB\\(-\\)", "", names(lbcon))



for (name in c("lbpos", "lbsur","lbneg","lbcon")){
  input <- get(name)
  outName <- paste0(name,"_split")
  split <- c()
  for (layer in names(input)) {
    i <- gsub("ayer","", layer)
    df <- input[[layer]]
    E3 <- df[df$avg_log2FC.spatial > 0.25,]
    E4 <- df[df$avg_log2FC.spatial < -0.25,]
    split[[paste0(i,"_E3")]] <- E3
    split[[paste0(i,"_E4")]] <- E4
  }
  assign(outName, split)
  rm(split, input, outName, layer, E3, E4)
}

for (name in  ls(pattern = "split$")) {
  df <- get(name)
  genesOnly <- lapply(df, function(x) x$gene)
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  p <- upset(fromList(genesOnly), 
             nintersects = NA, mainbar.y.label = paste0(gsub("_split","",name), " Intersection Size - Ex deconvoluted"),
             sets.bar.color = layerAPOEcols[names(genesOnly)], 
             sets = rev(names(genesOnly)), 
             keep.order = T,
             order.by = "freq", 
             point.size = 3, 
             main.bar.color = "#0c6b58", 
             text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
             #boxplot.summary = c("L1_E3","L2_E3")
  )
  numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
  numLabels <- length(p$labels)
  fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_Ex_deconvoluted_",name,".pdf")
  pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
  print(p)
  dev.off()
  rm(p, fileName, genesOnly, df)
}

### get the venn gene matrix:
wb <- createWorkbook()
for (name in  ls(pattern = "split$")) {
  df <- get(name)
  genesOnly <- lapply(df, function(x) x$gene)
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  elements <- unique(unlist(genesOnly))
  data <- unlist(lapply(genesOnly, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(genesOnly)
  rownames(data) <- elements
  addWorksheet(wb, name)
  writeData(wb = wb, sheet = name, x = data, , rowNames = T)
}
fileName <- paste0(outdir,timeStamp,"_Ex_Deconvoluted_DESeq2_common_UpsetPlot_byAnnotation_genes.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

mergedGenes <- c()
for (name in  ls(pattern = "split$")) {
  df <- get(name)
  names(df) <- paste0(gsub("_split","",name),"_",names(df))
  genesOnly <- lapply(df, function(x) x$gene)
  mergedGenes <- c(mergedGenes, genesOnly)
}

### for each layer
for (layer in names(layerCols)) {
  idx <- grep(layer, names(mergedGenes))
  genesOnly <- mergedGenes[idx]
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  idx2 <- paste0('^.*_', layer)
  p <- upset(fromList(genesOnly), 
             nintersects = NA, 
             mainbar.y.label = paste0(layer, " Intersection Size"),
             sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
             sets = rev(names(genesOnly)), 
             keep.order = T,
             order.by = "freq", 
             point.size = 3, 
             main.bar.color = "#0c6b58", 
             text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
             #boxplot.summary = c("L1_E3","L2_E3")
  )
  numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
  numLabels <- length(p$labels)
  fileName <- paste0(outdir,timeStamp,"_DESeq2_common_UpsetPlot_Ex_deconvoluted_",layer,".pdf")
  pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
  print(p)
  dev.off()
}

### get the venn gene matrix:
wb <- createWorkbook()
for (layer in names(layerCols)) {
  idx <- grep(layer, names(mergedGenes))
  genesOnly <- mergedGenes[idx]
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  elements <- unique(unlist(genesOnly))
  data <- unlist(lapply(genesOnly, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(genesOnly)
  rownames(data) <- elements
  addWorksheet(wb, layer)
  writeData(wb = wb, sheet = layer, x = data, rowNames = T)
}
fileName <- paste0(outdir,timeStamp,"_Ex_Deconvoluted_DESeq2_common_UpsetPlot_byLayer_genes.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)


### pathway upsets
file <- "2405 DEseq LB+ spots L2-6 E3 E4-cano combined.xlsx"
sheets <- getSheetNames("2405 DEseq LB+ spots L2-6 E3 E4-cano combined.xlsx")
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df$Summary <- NULL
df$Backward <- NULL
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
)
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_CanoCombined_UpsetPlot_LB+_L2-6_E3_E4.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)
fileName <- paste0(outdir,timeStamp,"_CanoCombined_UpsetPlot_LB+_L2-6_E3_E4.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

file <- "240820 DEseq LBsur spots L2-6 E3 E4-cano sum.xlsx"
sheets <- getSheetNames(file)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
)
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBsur_L2-6_E3_E4.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBsur_L2-6_E3_E4.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

file <- "240823 DEseq LB- spots L2-6 E3 E4-cano sum.xlsx"
sheets <- getSheetNames(file)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
)
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBneg_L2-6_E3_E4.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBneg_L2-6_E3_E4.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

file <- "240824 DEseq CLB- spots L2-6 E3 E4-cano sum.xlsx"
sheets <- getSheetNames(file)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
)
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_CLBneg_L2-6_E3_E4.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_CLBneg_L2-6_E3_E4.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)

### pathway upsets
file <- "pathways/E3 vs E4 combined Ctrl and LBD L1-6, WM-E3, E4-YJ.xlsx"
sheets <- getSheetNames(file)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df$Summary <- NULL
df$Backward <- NULL
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
); p
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_E3vsE4_combinedCtrlandLBD.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)

patterns <- unname(unlist(unique(unite(p$New_data, sep = ".", col = "pattern"))))
tmp <- unite(data, sep = ".", col = "pattern")
data$intersection <- tmp$pattern

outpaths <- lapply(patterns, function(pattern) {
  pathways <- data.frame(rownames(data[data$intersection == pattern,]))
  colnames(pathways) <- paste(colnames(data)[c(
    as.logical(
      as.numeric(
        unlist(
          strsplit(pattern,".", fixed = T)
        )
      )
    ),
    FALSE)], 
    collapse = ".")
  
  pathways
  addWorksheet(wb, pattern)
  writeData(wb, sheet = pattern ,pathways,rowNames=F,colNames=T)
})
names(outpaths) <- patterns


fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_E3vsE4_combinedCtrlandLBD.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)


### lbd vs control 
file <- "pathways/240906 LBD vs Ctrl combined genotype canonical pathway.xlsx"
sheets <- getSheetNames(file)
outname <- basename(gsub(".xlsx","",file))
df <- lapply(sheets, function(x) {
  read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
})
names(df) <- sheets
df$Summary <- NULL
df$Backward <- NULL
df <- lapply(df, function(x) x$Ingenuity.Canonical.Pathways)

p <- upset(fromList(df), 
           nintersects = NA, 
           mainbar.y.label = "Intersection Size",
           #sets.bar.color = layerAPOEcols[gsub(idx2,"WM",names(genesOnly), perl = T)], 
           sets = rev(names(df)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
); p
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBDvsCtrl_combinedGenotype.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

elements <- unique(unlist(df))
data <- unlist(lapply(df, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(df), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(df)
rownames(data) <- elements

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)

patterns <- unname(unlist(unique(unite(p$New_data, sep = ".", col = "pattern"))))
tmp <- unite(data, sep = ".", col = "pattern")
data$intersection <- tmp$pattern

outpaths <- lapply(patterns, function(pattern) {
  pathways <- data.frame(rownames(data[data$intersection == pattern,]))
  colnames(pathways) <- paste(colnames(data)[c(
    as.logical(
      as.numeric(
        unlist(
          strsplit(pattern,".", fixed = T)
        )
      )
    ),
    FALSE)], 
    collapse = ".")
  
  pathways
  addWorksheet(wb, pattern)
  writeData(wb, sheet = pattern ,pathways,rowNames=F,colNames=T)
})
names(outpaths) <- patterns


fileName <- paste0(outdir,timeStamp,"_Cano_UpsetPlot_LBDvsCtrl_combinedGenotype.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)


#### SNCA high vs low 
### deconvoluted ex 
files <- list.files(path = "DEGs/", pattern = "092324.*common.xlsx", full.names = T)
for (file in files){
  #file <- "DEGs/DESeq2_common_deconvolutedCelltypes_byLBAnnotation_APOE_Layers_marker_genes_LFC-consistent.xlsx"
  sheets <- getSheetNames(file)
  selectedSheets <- grep("Gray",sheets, value = T, invert = T)
  outname <- basename(gsub(".xlsx","",file))
  df <- lapply(selectedSheets, function(x) {
    read.xlsx(file, sheet = x, skipEmptyRows = T, skipEmptyCols = T)
  })
  names(df) <- selectedSheets
  for (name in names(df)){
    if(is.null(df[[name]])){
      df[[name]] <- NULL
    }
  }
  
  splitDF25 <- c()
  for (layer in names(df)) {
    i <- gsub("ayer","", layer)
    outName <- paste0(i,"_",gsub("high","",str_split_i(outname, "_", 3)))
    tmp <- df[[layer]]
    High <- tmp[tmp$avg_log2FC.spot > 0.25,]
    Low <- tmp[tmp$avg_log2FC.spot < -0.25,]
    splitDF25[[paste0(i,"_High")]] <- High
    splitDF25[[paste0(i,"_Low")]] <- Low
    #assign(outName, splitDF)
  }
  rm(outName, tmp, High, Low)

  
  genesOnly <- lapply(splitDF25, function(x) x$gene)
  genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
  p <- upset(fromList(genesOnly), 
             nintersects = NA, mainbar.y.label = paste0(gsub("high","",str_split_i(outname, "_", 3)), " Intersection Size"),
             #sets.bar.color = layerAPOEcols[names(genesOnly)], 
             sets = rev(names(genesOnly)), 
             keep.order = T,
             order.by = "freq", 
             point.size = 3, 
             main.bar.color = "#0c6b58", 
             text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
             #boxplot.summary = c("L1_E3","L2_E3")
  )
  numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
  numLabels <- length(p$labels)
  fileName <- paste0(outdir,outname,".pdf")
  pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
  print(p)
  dev.off()
  #rm(fileName, df)

  ### get genes
  elements <- unique(unlist(genesOnly))
  data <- unlist(lapply(genesOnly, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(genesOnly)
  rownames(data) <- elements
  
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeData(wb = wb, sheet = "Summary", x = data, rowNames = T)
  
  patterns <- unname(unlist(unique(unite(p$New_data, sep = ".", col = "pattern"))))
  tmp <- unite(data, sep = ".", col = "pattern")
  data$intersection <- tmp$pattern
  
  outGenes <- lapply(patterns, function(pattern) {
    genes <- data.frame(rownames(data[data$intersection == pattern,]))
    colnames(genes) <- paste(colnames(data)[c(
      as.logical(
        as.numeric(
          unlist(
            strsplit(pattern,".", fixed = T)
          )
        )
      ),
      FALSE)], 
      collapse = ".")
    
    genes
    addWorksheet(wb, pattern)
    writeData(wb, sheet = pattern ,genes,rowNames=F,colNames=T)
  })
  names(outGenes) <- patterns
  
  fileName <- paste0(outdir,outname,".xlsx")
  saveWorkbook(wb, file = fileName, overwrite = T)
  
}

### Layer 5 specific
layer <- "L5"
idx <- grep(layer, names(mergedGenes), value = T)
idx <- grep("con", idx, value=T, invert = T)
genesOnly <- mergedGenes[idx]
genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
plotCols <- layerAPOEcols[gsub(".*_L","L",names(genesOnly), perl = T)]

names(genesOnly) <- c("LB- Up", "LB- Down", "LB+ Up", "LB+ Down", "LBsur Up", "LBsur Down")

p <- upset(fromList(genesOnly), 
           nintersects = NA, 
           mainbar.y.label = paste0(layer, " Intersection Size"),
           sets.bar.color = plotCols, 
           sets = rev(names(genesOnly)), 
           keep.order = T,
           order.by = "freq", 
           point.size = 3, 
           main.bar.color = "#0c6b58", 
           text.scale = c(1.5,.8,1.5,.8,2,1) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
           #boxplot.summary = c("L1_E3","L2_E3")
); p
numIntersections <- nrow(unique(unite(p$New_data, sep = ".", col = "pattern")))
numLabels <- length(p$labels)
fileName <- paste0(timeStamp,"_DESeq2_common_UpsetPlot_",layer,"noCLB.pdf")
pdf(file = fileName,width = numIntersections/2, height = numLabels*(3/4)+1)
print(p)
dev.off()

wb <- createWorkbook()

genesOnly <- mergedGenes[idx]
genesOnly <- genesOnly[!unlist(lapply(genesOnly, is_empty))]
elements <- unique(unlist(genesOnly))
data <- unlist(lapply(genesOnly, function(x) {
  x <- as.vector(match(elements, x))
}))
data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(genesOnly), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- names(genesOnly)
rownames(data) <- elements
addWorksheet(wb, layer)
writeData(wb = wb, sheet = layer, x = data, rowNames = T)

fileName <- paste0(timeStamp,"_DESeq2_common_UpsetPlot_L5noCLB_genes.xlsx")
saveWorkbook(wb, file = fileName, overwrite = T)
