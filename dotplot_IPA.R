library("tidyverse")
library("data.table")
library(openxlsx)

layerCols <- c( 
  "L1" ="#8D405C",
  "L23"= "#E7BDE1",
  "L4"= "#CF8CA4",
  "L5"= "#9F6E80",
  "L6"= "#CDADB9",
  "WM"= "#67A9D8")

apoeColors <- c("E3"="#2367AC", 
                "E4"="#B21F2C")

## The following codes will be used to create a dotplot of the pathways 
## comparing among the APOE genotypes in control and AD subjects in 
## excitatory neurons

Root.dir <- "./pathways"
input.file <- paste0(Root.dir,"/2405 DEseq LB+ spots L2-6 E3 E4-toxic combined_AQW.xlsx")
sheetnames <- getSheetNames(input.file)
BP.df.list <- read.xlsx(input.file, sheet = "Summary")
BP.df.list <- na.omit(BP.df.list)


BP.df.list$Group <- paste(BP.df.list$Layer, BP.df.list$APOE, sep = ".")
colnames(BP.df.list)[4] <- "Pval"
BP.df.list$Ratio <- ifelse(BP.df.list$APOE == "E4", BP.df.list$Ratio*-1, BP.df.list$Ratio)

BP.df.list$Ingenuity.Toxicity.Lists <- factor(BP.df.list$Ingenuity.Toxicity.Lists, 
                                                  levels = rev(unique(BP.df.list$Ingenuity.Toxicity.Lists)))

Plot.Path <- ggplot(aes(x = Group, 
                        y =Ingenuity.Toxicity.Lists),
                    data = BP.df.list) +
  #xlab("-log10 [FDR]") +
  # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_point(aes(fill = Ratio, size = Pval, group=APOE), shape = 21) +
  scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = 90),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))


setwd(Root.dir)
ggsave("061324_LBpos.Tox.APOE_correctColors.pdf", width = 5*2, height = 2.7*2, plot = Plot.Path)

### canonical
input.file <- paste0(Root.dir,"/2405 DEseq LB+ spots L2-6 E3 E4-cano combined_AQW.xlsx")
sheetnames <- getSheetNames(input.file)
BP.df.list <- read.xlsx(input.file, sheet = "Summary")
#BP.df.list[,8] <- NULL
BP.df.list <- na.omit(BP.df.list)

BP.df.list$Group <- paste(BP.df.list$Layer, BP.df.list$APOE, sep = ".")
colnames(BP.df.list)[4] <- "Pval"
BP.df.list$Ratio <- ifelse(BP.df.list$APOE == "E4", BP.df.list$Ratio*-1, BP.df.list$Ratio)

BP.df.list$Ingenuity.Canonical.Pathways <- factor(BP.df.list$Ingenuity.Canonical.Pathways, 
                                                  levels = rev(unique(BP.df.list$Ingenuity.Canonical.Pathways)))


Plot.Path <- ggplot(aes(x = Group, 
                        y =Ingenuity.Canonical.Pathways),
                    data = BP.df.list) +
  #xlab("-log10 [FDR]") +
  # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_point(aes(fill = Ratio, size = Pval, group=APOE), shape = 21) +
  scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = 90),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")); Plot.Path


setwd(Root.dir)
ggsave("061324_LBpos.Cano.APOE_correctColors.pdf", width = 5*2, height = 5*2, plot = Plot.Path)

### deconvoluted data
Root.dir <- "./pathways"
input.file <- paste0(Root.dir,"/SummaryLBannoDeconDEseq-EXE3E4canonical_fin.xlsx")
sheetnames <- getSheetNames(input.file)
BP.df.list <- lapply(sheetnames, function(sheet) read.xlsx(input.file, sheet = sheet))
names(BP.df.list) <- sheetnames
#BP.df.list <- read.xlsx(input.file, sheet = "Summary")
BP.df.list <- na.omit(BP.df.list)
for (sheet in sheetnames){
  BP.df.list[[sheet]]$sheet <- sheet
}

toPlot <- rbind(BP.df.list$`LB+ EX`, BP.df.list$`LBsur EX`, BP.df.list$`LB- EX`)
toPlot$Annotation <- gsub(" EX","", toPlot$sheet)
toPlot$Group <- paste(toPlot$APOE, toPlot$Layer, toPlot$Annotation,  sep = ".")
toPlot$Group <- factor(toPlot$Group, levels = unique(toPlot$Group))
colnames(toPlot)[4] <- "Pval"
toPlot$Ratio <- ifelse(toPlot$APOE == "E4", toPlot$Ratio*-1, toPlot$Ratio)
toPlot <- toPlot %>% group_by(Group) %>% arrange(desc(Pval))

toPlot$Ingenuity.Canonical.Pathways <- factor(toPlot$Ingenuity.Canonical.Pathways, 
                                              levels = rev(unique(toPlot$Ingenuity.Canonical.Pathways)))
#toPlot <- toPlot[toPlot$Layer == "L5",]
Plot.Path <- ggplot(aes(x = Group, 
                        y =Ingenuity.Canonical.Pathways),
                    data = toPlot) +
  #xlab("-log10 [FDR]") +
  # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_point(aes(size = Pval,group=APOE, color = Ratio)) +
  scale_color_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12, angle = 90),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) + 
  xlab("L5"); Plot.Path


setwd(Root.dir)
ggsave("082924_LBpos-LBsur-LBneg.Canon.EX-L5-APOE_correctColors_combined.pdf", width = 5*2, height = (nrow(Plot.Path$data)/12)+4, plot = Plot.Path)

 Plot.Path <- ggplot(aes(x = Group, 
                         y =Ingenuity.Canonical.Pathways),
                     data = toPlot) +
  #xlab("-log10 [FDR]") +
    #geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
   geom_point(aes(size = Pval, group=APOE, fill = Ratio), shape=21) +
   scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
   theme_bw() +
   theme(panel.grid = element_blank(),
         panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
         axis.title.y = element_blank(),
         axis.text.y = element_text(color = "black", size = 12),
         axis.text.x = element_text(color = "black", size = 12, angle = 90),
         axis.line = element_line(color = "black"),
         axis.ticks = element_line(color = "black")) ; Plot.Path
 
 
 setwd(Root.dir)
 ggsave("071324_LBpos-LBsur-LBneg.Canon.EX-L5-APOE_correctColors_nosplit.pdf", width = 5*2, height = 2.7*2, plot = Plot.Path)

 toPlot <- rbind(BP.df.list$`LB+ EX`, BP.df.list$`LBsur EX`, BP.df.list$`LB- EX`)
 toPlot$Annotation <- gsub(" EX","", toPlot$sheet)
 toPlot$Group <- paste(toPlot$APOE, toPlot$Annotation,  sep = ".")
 toPlot <- toPlot[toPlot$Layer == "L5",]
 toPlot$Group <- factor(toPlot$Group, levels = unique(toPlot$Group))
 colnames(toPlot)[4] <- "Pval"
 toPlot$Ratio <- ifelse(toPlot$APOE == "E4", toPlot$Ratio*-1, toPlot$Ratio)
 toPlot <- toPlot %>% group_by(Group) %>% arrange(desc(Pval), .by_group = T)
 
 toPlot$Ingenuity.Canonical.Pathways <- factor(toPlot$Ingenuity.Canonical.Pathways, 
                                               levels = rev(unique(toPlot$Ingenuity.Canonical.Pathways)))

 Plot.Path <- ggplot(aes(x = Group, 
                         y =Ingenuity.Canonical.Pathways),
                     data = toPlot) +
   #xlab("-log10 [FDR]") +
   # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
   geom_point(aes(size = Pval,group=APOE, fill = Ratio), shape=21) +
   scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
   theme_bw() +
   theme(panel.grid = element_blank(),
         panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
         axis.title.y = element_blank(),
         axis.text.y = element_text(color = "black", size = 12),
         axis.text.x = element_text(color = "black", size = 12, angle = 90),
         axis.line = element_line(color = "black"),
         axis.ticks = element_line(color = "black")) + 
   xlab("L5"); Plot.Path
 
 
 setwd(Root.dir)
 ggsave("082924_LBpos-LBsur-LBneg.Canon.EX-L5-APOE_correctColors_combined_ordered.pdf", width = 5*2, height = (nrow(Plot.Path$data)/12)+4, plot = Plot.Path)
 
 
 
toPlot <- rbind(BP.df.list$`LB- EX`, BP.df.list$`LB- MG`, BP.df.list$`LB- OLG`,
                BP.df.list$`LB- OPC`, BP.df.list$`LB- AS`)
toPlot$Annotation <- gsub(" .*$","", toPlot$sheet)
toPlot$Cell <- gsub("^.* ","", toPlot$sheet)
toPlot$Group <- paste(toPlot$Cell, toPlot$APOE, sep = ".")
colnames(toPlot)[4] <- "Pval"
toPlot$Ratio <- ifelse(toPlot$APOE == "E4", toPlot$Ratio*-1, toPlot$Ratio)
toPlot[grep("-alanine",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "beta-alanine Degradation"
toPlot[grep("-Adrenergic",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "alpha-Adrenergic Signaling"
toPlot[grep("HIF1",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "HIF1-alpha"
toPlot$Ingenuity.Canonical.Pathways <- factor(toPlot$Ingenuity.Canonical.Pathways, 
                                              levels = rev(unique(toPlot$Ingenuity.Canonical.Pathways)))
#toPlot <- toPlot[toPlot$Layer == "L5",]
Plot.Path <- ggplot(aes(x = Group, 
                        y =Ingenuity.Canonical.Pathways),
                    data = toPlot) +
  #xlab("-log10 [FDR]") +
  # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_point(aes(size = Pval, group=APOE, color = Ratio)) +
  scale_color_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
  #scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 90),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")); Plot.Path


setwd(Root.dir)
ggsave("082724_LBnegCelltypes.Canon.APOE_correctColors_combined.pdf", width = 5*2, height = (nrow(Plot.Path$data)/14)+2, plot = Plot.Path)


toPlot <- rbind(BP.df.list$`LB- EX`, BP.df.list$`LB- MG`, BP.df.list$`LB- OLG`,
                BP.df.list$`LB- OPC`, BP.df.list$`LB- AS`)
toPlot$Annotation <- gsub(" .*$","", toPlot$sheet)
toPlot$Cell <- gsub("^.* ","", toPlot$sheet)
toPlot$Group <- paste(toPlot$Cell, toPlot$APOE, sep = ".")
colnames(toPlot)[4] <- "Pval"
toPlot$Ratio <- ifelse(toPlot$APOE == "E4", toPlot$Ratio*-1, toPlot$Ratio)
toPlot[grep("-alanine",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "beta-alanine Degradation"
toPlot[grep("-Adrenergic",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "alpha-Adrenergic Signaling"
toPlot[grep("HIF1",toPlot$Ingenuity.Canonical.Pathways),"Ingenuity.Canonical.Pathways"] <- "HIF1-alpha"
toPlot$Ingenuity.Canonical.Pathways <- factor(toPlot$Ingenuity.Canonical.Pathways, 
                                              levels = rev(unique(toPlot$Ingenuity.Canonical.Pathways)))
for (layer in unique(toPlot$Layer)){
  toPlo <- toPlot[toPlot$Layer %in% layer,]
  #toPlot$Group <- paste(toPlot$Layer, toPlot$APOE, toPlo$ sep = ".")
  Plot.Path <- ggplot(aes(x = Group, 
                          y =Ingenuity.Canonical.Pathways),
                      data = toPlo) +
    #xlab("-log10 [FDR]") +
    # geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_point(aes(size = Pval, group=APOE, color = Ratio)) +
    scale_color_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
    #scale_fill_gradient2(low = apoeColors[2] , mid = "white", high = apoeColors[1]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_line(linetype = "dotted", color = "#606060"),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black", size = 8),
          axis.text.x = element_text(color = "black", size = 8, angle = 90),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black")) +
    xlab(paste0(layer, " LB-")); Plot.Path
  setwd(Root.dir)
  ggsave(paste0("082724_LBnegCelltypes.Canon.APOE_correctColors_",layer,".pdf"), width = 5*2, height = (nrow(Plot.Path$data)/12)+2, plot = Plot.Path)
}



