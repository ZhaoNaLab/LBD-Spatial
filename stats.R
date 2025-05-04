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
library(patchwork)
library(ggstatsplot, lib.loc = "./Rlib4.2.2/")
library(prismatic, lib.loc = "./Rlib4.2.2/")
library(lme4)
library(afex, lib.loc = "./Rlib4.2.2/")
#library(reshape)
library(colorRamps)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggridges)


require(GGally)
require(reshape2)
require(compiler)
require(parallel)
require(boot)
require(lattice)

source("../../scripts/helper.r")

layerCols2 <- c("Layer1" ="#8D405C", 
                "Layer23"= "#E7BDE1",
                "Layer4"= "#CF8CA4",
                "Layer5"= "#9F6E80",
                "Layer6"= "#CDADB9",
                "WM"= "#67A9D8"
                #"Unknown" = "#33A02C"
                )

lbAnnoCols <- c("LB(+)" ="green", 
                "LB Surround"= "purple",
                "LB(-)"= "grey",
                "CLB(-)"= "grey")

lbAnnoCols2 <- c("WM_LB(+)" ="green", 
                 "WM"= "purple",
                 "Unselected"= "grey30",
                 "WM_CLB(-)"= "purple")
lbAnnoCols3 <- c("WM_LB(+)" ="green", 
                 "WM_LB_Surround"= "purple",
                 "Unselected"= "grey",
                 "Layer6" = "red",
                 "WM"= "white")

lbd_images <- c("LBD_3","LBD_4", "LBD_6","LBD_8","LBD_11", "LBD_12")
lbdImageOrder <- c("LBD_4","LBD_11", "LBD_12", "LBD_3","LBD_6","LBD_8")

### redo with dropped unknowns
data2 <- readRDS("bothCohorts_spatial_withAdjExpAssays_dropUnknowns_correctedWeight_fixedImageRotation_051024.rds")

data2@meta.data$LBannotation <- factor(data2@meta.data$LBannotation, levels = c("CLB(-)",
                                                                                "LB Surround",
                                                                                "LB(+)",
                                                                                "LB(-)"))

imagescales <- getImagePointSizes(data2)

modelDF <- data2@meta.data[data2@meta.data$LBannotation != "CLB(-)",c("disease_state", "LBannotation", "sampleid", 
                                                                      "CorrectedLayers", "apoe_geno", "slide_serial")]
modelDF <- data2@meta.data[,c("disease_state", "LBannotation", "sampleid", 
                                                                      "CorrectedLayers", "apoe_geno", "slide_serial")]

table(modelDF$LBannotation, modelDF$CorrectedLayers)
#modelDF$LBannotation <- factor(modelDF$LBannotation)
#data(obk.long, package = "afex")
#modelDF$LBannotation <- ifelse(modelDF$LBannotation == "LB(+)", 1,0)

modelDF$CorrectedLayers <- factor(modelDF$CorrectedLayers, levels = names(layerCols2))
modelDF$disease_state <- factor(modelDF$disease_state, levels = c("Tri","LBD"))
modelDF$apoe_geno <- factor(modelDF$apoe_geno, levels = c("E3/3","E3/4"))

# estimate mixed ANOVA on the full design:
# aov_ez(id = "sampleid", "LBannotation", data = modelDF, 
#        between = "slide_serial", 
#        within = c("disease_state","apoe_geno"),
#        observed = c("LBannotation","CorrectedLayers"))

##### using mike's suggestions
# Save the data in two different vector
#before <- mice2$before
#after <- mice2$after
totals <- table(modelDF$sampleid)
#modelDF$sampleid <- factor(modelDF$sampleid, levels = lbd_images)
tmp <- reshape2::melt(table(modelDF$LBannotation, modelDF$CorrectedLayers, modelDF$sampleid, modelDF$apoe_geno))
tmp <- tmp[tmp$Var1 == "LB(+)",]
names(tmp) <- c("LBannotation", "Layer","Sample","APOE","Count")
apoeKey <- melt(table(data2$sampleid, data2$apoe_geno))
apoeKey <- apoeKey[apoeKey$value != 0,c(1,2)]
tmp <- lapply(apoeKey$Var1, function(i) tmp[tmp$Sample == i & tmp$APOE == apoeKey[apoeKey$Var1 == i,2],])
tmpDF <- tmp[[1]]
for (i in 2:10){
  tmpDF <- rbind(tmpDF, tmp[[i]])
}
tmp <- tmpDF
### percentage of total spots per sample
for( i in names(totals)){
  tmp[tmp$Sample == i, "Percentage"] <- tmp[tmp$Sample == i, "Count"]/totals[[i]]
}
tmp$Percentage <- ifelse(is.nan(tmp$Percentage), 0, tmp$Percentage )
### percentage of LB+ spots per sample
for( i in names(totals)) {
  tmpTotal <- sum(tmp[tmp$Sample == i, "Count"])
  tmp[tmp$Sample == i, "PercentageLB"] <- tmp[tmp$Sample == i, "Count"]/tmpTotal
}

### percentage of spots per layer per samples
for( i in names(totals)){
  layerTotals <- melt(table(modelDF[modelDF$sampleid == i, "CorrectedLayers"]))
  for (layer in unique(tmp$Layer)){
    tmp[tmp$Sample == i & tmp$Layer == layer, "LayerPercentage"] <- tmp[tmp$Sample == i & tmp$Layer == layer, "Count"]/layerTotals[layerTotals$Var1 == layer, "value"]
  }
}

tmp$PercentageLB <- ifelse(is.nan(tmp$PercentageLB), 0, tmp$PercentageLB )
# Compute t-test
### can only do 2 groups at a time...
library(tidyverse)
library(rstatix)
library(ggpubr) 

#tmp <- tmp[!tmp$Layer=="WM",]
### drop controls
controls <- c("LBD_1","LBD_5","LBD_7","LBD_9")
tmpNoC <- tmp[!tmp$Sample %in% controls,]
tmpNoC <- tmpNoC[!tmpNoC$Layer == "WM",]
for (i in c("Count","Percentage","PercentageLB","LayerPercentage")){
  outName <- paste0("LB_stats/paired_t_test_",i,"_GroupByLayer.csv")
  outName2 <- paste0("LB_stats/paired_t_test_",i,"_GroupByLayer.png")
  stat.test <- tmpNoC %>%
    group_by(Layer) %>%
    pairwise_t_test(
      as.formula(paste0(i, "~ APOE")), paired = TRUE, 
      p.adjust.method = "none", 
    )
  # Create the plot
  bxp <- ggboxplot(
    tmpNoC, x = "Layer", y = i,
    color = "APOE", palette = "jco", outlier.shape = 2, add = c("jitter"), 
  )
  # Add statistical test p-values
  stat.test <- stat.test %>% add_xy_position(x = "Layer", group = "APOE")
  stat.test$groups <- NULL
  write.csv(stat.test, file = outName)
  
  bxp <- bxp + stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    step.increase = 0.08, hide.ns = T
  )
  saveImages(outName2, bxp)
  
  outName <- paste0("LB_stats/paired_t_test_",i,"_GroupByAPOE.csv")
  outName2 <- paste0("LB_stats/paired_t_test_",i,"_GroupByAPOE.png")
  stat.test <- tmpNoC %>%
    group_by(APOE) %>%
    pairwise_t_test(
      as.formula(paste0(i, "~ Layer")), paired = TRUE, 
      p.adjust.method = "none"
    )
  # Create the plot
  bxp <- ggboxplot(
    tmpNoC, x = "Layer", y = i,
    color = "APOE", palette = "jco", outlier.shape = 2
  )
  # Add statistical test p-values
  stat.test <- stat.test %>% add_xy_position(x = "Layer", group = "APOE")
  stat.test$groups <- NULL
  write.csv(stat.test, file = outName)
  
  bxp <- bxp + stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    step.increase = 0.08, hide.ns = T
  )
  saveImages(outName2, bxp)
}

# compute the difference
d <- with(tmpNoC, 
          Count[APOE == "E3/3"] - Count[APOE == "E3/4"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.6141

################# no White matter
tmpNoC <- tmp[!tmp$Sample %in% controls,]
tmpNoC <- tmpNoC[!tmpNoC$Layer == "WM",]
tmpNoC$Layer <- factor(tmpNoC$Layer, levels = c("Layer1","Layer23","Layer4","Layer5","Layer6"))
for (i in c("Count","Percentage","PercentageLB","LayerPercentage")){
  outName <- paste0("LB_stats/paired_t_test_",i,"_OnlyLayer.csv")
  outName2 <- paste0("LB_stats/paired_t_test_",i,"_OnlyLayer.png")
  stat.test <- tmpNoC %>%
    pairwise_t_test(as.formula(paste0(i," ~ Layer")),
                    paired = TRUE, 
      p.adjust.method = "none", 
    )
  # Create the plot
  bxp <- ggboxplot(
    tmpNoC, x = "Layer", y = i,
    color = "Layer", palette = layerCols2, outlier.shape = 2, add = c("jitter"), xlab = FALSE
  )
  # Add statistical test p-values
  stat.test <- stat.test %>% add_xy_position(x = "Layer")
  stat.test$groups <- NULL
  write.csv(stat.test, file = outName)
  
  bxp <- bxp + stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    step.increase = 0.08, hide.ns = T
  )
  saveImages(outName2, bxp)
}




