#######################
# A
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(randomcoloR)



CBC=c('Lgr5','Ascl2','Axin2','Slc12a2','Gkn3','Olfm4')
TA=c('Cdk4','Mcm5','Mcm6','Pcna','Mki67')
SPC=c('Atoh1','Dll1','Mist1','Spdef')
APC=c('Hes1','Math1','Klf4','Cdx2')
GC=c('Anxa3','Gfi1','Spdef','Muc2','Tff3','Agr2','Clca3')
PC=c('Muc2','Lyz1','Defa17','Defa22','Defa24','Defa4','Ang4')
EC=c('Arg2','Il18','Car4','Ccl25','Alpi')
EE=c('Chga','Chgb','Tac1','Tph1','Neurog3')
TC=c('Dclk1','Trpm5','Gfi1b','Il25')
Tcell=c('Cd3e','Cd3d','Cd3g','Cd8a','Prf1','Cd4','Foxp3','Ctla4')
macrophage.DC=c('Ly86','Cd74','Cd68','H2-DMb2','Ms4a4c')

ann.fm = unique(c(CBC,TA,SPC,APC,EC,GC,PC,EE,TC,Tcell,macrophage.DC))

DefaultAssay(s.integrated) = "RNA"


mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','coral1','magenta')

tiff(filename = sprintf("%s/figure/Figure5A.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
DimPlot(s.integrated, group.by = 'cellType', label=F, cols = mycol) + ggtitle('')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()

tiff(filename = sprintf("%s/figure/Figure5A.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
DimPlot(s.integrated, group.by = '', label=F, cols = mycol) + ggtitle('')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()



#######################
# B
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(randomcoloR)


mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','coral1','magenta')


cont.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "Cont") 

pge2.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "PGE2")

ABC99.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "ABC99")

cell.pop = rbind(cont.pop,pge2.pop,ABC99.pop)
cell.pop$group = as.factor(rep(c(1, 2, 3), c(nrow(cont.pop), nrow(pge2.pop), nrow(ABC99.pop))))
#cell.pop$cellType = factor(cell.pop$cellType, levels = levels(factor(cell.pop$cellType)))
cell.pop$cellType = factor(cell.pop$cellType, levels = rev(unique(cell.pop$cellType)))

ggplot(cell.pop, aes(x = group, y = percent, fill = cellType, width = 0.7)) + 
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  #geom_text(aes(label = paste0(round(cell.pop$percent, 1), "%")), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = rev(mycol)) + 
  labs(y = "% of cell population") + 
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") + 
  geom_vline(xintercept = 0.5, color = "black", size = 0.5, linetype = "solid") + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(face = "bold", color = "black", size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(), 
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +  
  scale_x_discrete(expand = c(0, 0.5)) 




#######################
# C
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(randomcoloR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOfuncR)


DefaultAssay(s.integrated) = "RNA"

tiff(filename = sprintf("%s/figure/Figure5C_PC.tiff", dir), width = 12, height = 25, units = 'cm', res = 300)
FeaturePlot(s.integrated, features = 'Lyz1', order = TRUE, split.by = 'sample', min.cutoff = 0, max.cutoff = 5, by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank())
dev.off()

tiff(filename = sprintf("%s/figure/Figure5C_GC.tiff", dir), width = 12, height = 25, units = 'cm', res = 300)
FeaturePlot(s.integrated, features = 'Muc2', order = TRUE, split.by = 'sample', min.cutoff = 0, max.cutoff = 3.5, by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank())
dev.off()



#######################
# D
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/Rdata/13-2_Trajectory-analysis.Rdata", dir)) #cds

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(monocle)
library(RColorBrewer)



mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','coral1','magenta')
scol = c('skyblue1','darkolivegreen2','tomato1')

tiff(filename = sprintf("%s/figure/supple_figure5D-1.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right", legend.text = element_text(size = 16))
dev.off()

tiff(filename = sprintf("%s/figure/supple_figure5D-2.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "State", show_branch_points = F) +
  scale_color_manual(values = scol) +
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
dev.off()

tiff(filename = sprintf("%s/figure/Figure5D-1.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "cellType", show_branch_points = F) + 
  scale_color_manual(values = mycol) +
  theme(legend.position = "right", legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3))) 
dev.off()




#######################
# E
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/Rdata/13-2_Trajectory-analysis.Rdata", dir)) #cds

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(monocle)
library(RColorBrewer)



mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','coral1','magenta')
scol = c('skyblue1','darkolivegreen2','tomato1')


df = pData(cds)
celltype.lev = levels(factor(df$cellType))
col.name = setNames(mycol[1:length(celltype.lev)], celltype.lev)
df$cellType = factor(df$cellType, levels = names(col.name))

tiff(filename = sprintf("%s/figure/Figure5E_cont.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State %in% c(3), sample == "Cont"),
       aes(x = Pseudotime, y = 2, color = cellType)) +
  geom_jitter(height = 0.3, size = 1, alpha = 0.5) +
  scale_color_manual(values = col.name, drop = FALSE) +
  scale_y_continuous(limits = c(1.8, 2.2), breaks = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        legend.position = "none",
        panel.grid.minor = element_blank()) 
dev.off()

tiff(filename = sprintf("%s/figure/Figure5E_pge2.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State %in% c(3), sample == "PGE2"),
       aes(x = Pseudotime, y = 2, color = cellType)) +
  geom_jitter(height = 0.3, size = 1, alpha = 0.5) +
  scale_color_manual(values = col.name, drop = FALSE) +
  scale_y_continuous(limits = c(1.8, 2.2), breaks = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        legend.position = "none",
        panel.grid.minor = element_blank()) 
dev.off()

tiff(filename = sprintf("%s/figure/Figure5E_inhibitor.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State %in% c(3), sample == "ABC99"),
       aes(x = Pseudotime, y = 2, color = cellType)) +
  geom_jitter(height = 0.3, size = 1, alpha = 0.5) +
  scale_color_manual(values = col.name, drop = FALSE) +
  scale_y_continuous(limits = c(1.8, 2.2), breaks = 1) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        legend.position = "none",
        panel.grid.minor = element_blank()) 
dev.off()





#######################
# J-1
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/Rdata/13-1_ABC99_vs_PGE2_DEG.Rdata", dir)) #cluster.degl

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(randomcoloR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOfuncR)
library(gridExtra)
library(SiPSiC)
library(SingleCellExperiment)



cluster.updeg = list()
i=1
for (i in 1:length(cluster.degl)){
  deg = cluster.degl[[i]]
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 0.5 & j$p_val < 0.05 & !is.na(j$p_val), ])
  cluster.updeg[[names(cluster.degl)[i]]] = updeg
}
lapply(cluster.updeg, function(i) {
  lapply(i, function(j) dim(j))
})

updeg = cluster.updeg$CBC$abc99.vs.pge2
ora = enrichGO(gene = rownames(updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
upgofilter = as.data.frame(gofilter(ora, level = 4))

upgo = head(upgofilter, 58)
upgo = upgo[!(grepl("catabolic|metabolic|regulation of", upgo$Description) &
                !grepl("negative regulation of cytokine production", upgo$Description)), ]
upgo = upgo[-c(2,3,9,12,13,17,21,22,23,24),]

upgo$pvalue = -log10(upgo$pvalue)
upgo = upgo[order(upgo$pvalue, decreasing = T), ]
upgo$Description = sapply(upgo$Description, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
ggplot(data = upgo, aes(x = reorder(Description, pvalue), y = pvalue, fill = pvalue)) +
  coord_flip() +
  #geom_hline(yintercept = -log10(0.05), color = "dimgrey", size = 0.6) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = "gold", mid = "gold", high = "gold", guide = FALSE, midpoint = 2) +
  ggtitle("ABC99 vs PGE2") +
  ylab("-log10 P-value") + xlab("") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(size = 0.1, color = "grey"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0.2, r = 0.5, b = 0.3, l = 0, unit = "cm")
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", size = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", size = 0.3) +
  geom_text(aes(label = Description, y = 0.05), hjust = 0, size = 4, color = "black")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,4)) +  
  scale_x_discrete(expand = c(0, 1))



#######################
# J-2
#######################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/Rdata/13-1_CBC_vs_others_DEG.Rdata", dir)) #markers.li


library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
#install.packages("dbplyr")
library(dbplyr)
library(colorRamp2)
library(viridis)
library(randomcoloR)
#install.packages("openai")
#remotes::install_github("Winnie09/GPTCelltype")
library(GPTCelltype)
library(openai)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOfuncR)
library(gridExtra)
library(stringr)
library(GO.db)
library(fgsea)
library(GOplot)
library(enrichplot)
library(grid)



cluster.updeg = list()
i=1
for (i in 1:length(markers.li)){
  deg = markers.li[[i]]  
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2.5 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}
lapply(cluster.updeg, function(i) {
  lapply(i, function(j) dim(j))
})


updeg = cluster.updeg$ABC99$CBC
upora = enrichGO(gene = rownames(updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
upgofilter = as.data.frame(gofilter(upora, level = 5))
upgofilter = upgofilter[!grepl('^negatvie regulation|$positive regulation|catabolic|metabolic', upgofilter$Description),]
upgofilter = upgofilter[!(grepl("catabolic|metabolic|regulation of", upgofilter$Description) &
                            !grepl("regulation of Wnt signaling pathway", upgofilter$Description)), ]
upgofilter = upgofilter[-c(2,4,5,8,9,10,13,14,16,18),]

top10 = head(upgofilter, 10)
top10$p.adjust = -log10(top10$p.adjust)
top10 = top10[order(top10$p.adjust, decreasing = TRUE), ]
top10$Description = sapply(top10$Description, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})


ggplot(data = top10, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = -log10(0.05), color = "red", size = 0.5) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = "lightpink", mid = "lightpink", high = "lightpink", guide = FALSE, midpoint = 2) +
  ggtitle('iNotum: CBC vs all others') +
  ylab("-log10FDR") + xlab("") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9), 
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(size = 0.1, color = "grey"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0.2, r = 0, b = 0.3, l = -0.5, unit = "cm"),
    legend.position = "none" 
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", size = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", size = 0.3) +
  geom_text(aes(label = Description, y = 0.05), hjust = 0, size = 4, color = "black")+
  scale_y_continuous(expand = c(0, 0)) +  
  scale_x_discrete(expand = c(0, 1))   

