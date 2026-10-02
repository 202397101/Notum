################
# Figure 2A
################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

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


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','lightslateblue','coral1')

tiff(filename = sprintf("%s/figure/Figure2B.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
DimPlot(s.integrated, group.by = 'cellType', label=F, cols = mycol) + ggtitle('')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()



################
# Figure 2B
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

#install.packages(c("Matrix", "SeuratObject"))
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


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','lightslateblue','coral1')


day4.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "4day") 

day1.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "1day")

ctrl.pop = s.integrated@meta.data %>% 
  group_by(sample, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(sample == "Ctrl")

cell.pop = rbind(ctrl.pop,day1.pop,day4.pop)
cell.pop$group = as.factor(rep(c(1, 2, 3), c(nrow(ctrl.pop), nrow(day1.pop), nrow(day4.pop))))
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



################
# Figure 2C
################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers
load(file = sprintf("%s/Rdata/11-2_Augur.Rdata", dir)) #augur.day4, augur.day1, augur.day1.day4


library('Seurat')
library('tidyverse')
library('ggplot2')
library(stringr)
library(dplyr)
library(ggrepel)
#devtools::install_github("neurorestore/Augur")
library(Augur)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
celltype.lev = levels(factor(s.integrated$cellType))
col.name = setNames(ifelse(celltype.lev %in% c("GC/PC", "CBC"), "red", "black"), celltype.lev)  # GC/PC, CBC만 빨강 강조
col.name

# Day1 vs Cont, Day4 vs Cont 의 cell type별 Augur AUC
augur.df = rbind(
  data.frame(augur.day1$AUC, comparison = "Day1 vs Cont"),
  data.frame(augur.day4$AUC, comparison = "Day4 vs Cont"))
augur.df$comparison = factor(augur.df$comparison, levels = c("Day1 vs Cont", "Day4 vs Cont"))
augur.df$cell_type = factor(augur.df$cell_type, levels = celltype.lev)
augur.df

tiff(filename = sprintf("%s/figure/Figure2C.tiff", dir), width = 9, height = 9, units = 'cm', res = 300)
ggplot(augur.df, aes(x = comparison, y = auc, group = cell_type, color = cell_type)) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed", linewidth = 0.3) +  # AUC 0.5 = 구분 못함
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 2) +
  geom_text_repel(data = augur.df %>% filter(comparison == "Day4 vs Cont"),
                  aes(label = cell_type), nudge_x = 0.3, direction = "y", hjust = 0,
                  size = 5, segment.size = 0.1, show.legend = FALSE) +
  scale_color_manual(values = col.name) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +  # 오른쪽에 라벨 공간
  labs(x = "", y = "Augur score (AUC)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.1))
dev.off()




################
# Figure 2D
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-2_FC_removeImmuneCells_scRNA.Rdata", dir)) # s.integrated, markers
load(file = sprintf("%s/Rdata/7-3_Day1_vs_Cont_mergeGC-PC_DEG.Rdata", dir)) #cluster.degl

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
library(tibble)
library(data.table)
library(msigdbr)



########## GC/PC, CBC 세포만 - cell type별 up DEG GO heatmap ##########
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
library(pheatmap)

# cell type별 GO 결과(upGO) - 위 CBC 전용 upGO와 이름이 겹치므로 별도 이름(cellGO)으로 불러옴
go.env = new.env()
load(file = sprintf("%s/Rdata/6-4_CBC_celltype_pathway_heatmap.Rdata", dir), envir = go.env) #upGO
cellGO = go.env$upGO
names(cellGO)  # cell type 이름 확인

sel.cell = c("CBC", "GC/PC")  # heatmap 열 순서
cellGO = cellGO[intersect(sel.cell, names(cellGO))]  # 이름이 다르면 여기서 빠지므로 names(cellGO) 확인
names(cellGO)

# 비교(rel.name)별로 선택 cell type의 top GO -> -log10 p 행렬
make.cellGO.mtx = function(rel.name) {
  top.go = lapply(cellGO, function(go) {
    go = data.frame(gofilter(go[[rel.name]], level = 4))
    go = go[go$pvalue < 0.05, ]
    go = go[order(go$pvalue), ]
    go = go[!grepl('negative|positive', go$Description), ]
    head(go, 7)
  })
  path = unique(unlist(lapply(top.go, function(i) i$Description)))

  mtx = matrix(1, nrow = length(path), ncol = length(cellGO), dimnames = list(path, names(cellGO)))
  for (cname in names(cellGO)) {
    go = data.frame(cellGO[[cname]][[rel.name]])
    go = go[go$Description %in% path, ]
    mtx[go$Description, cname] = -log10(go$pvalue)
  }
  mtx
}

# day1
day1.mtx = make.cellGO.mtx('1day.vs.Ctrl')
tiff(filename = sprintf("%s/figure/Figure2D_day1_GCPC_CBC_upGO.tif", dir), width = 25, height = 5 + 0.75 * nrow(day1.mtx), units = 'cm', res = 300)
pheatmap(day1.mtx,
         cluster_rows = T,
         cluster_cols = F,  # 2개 cell type은 sel.cell 순서로 고정
         color = colorRampPalette(c("white", "salmon", "red3"))(100),
         breaks = seq(0,4, length.out = 101),
         display_numbers = ifelse(day1.mtx > 1.5, "*", ""),
         number_color = "black",
         fontsize_number = 15,
         main = "Day1 vs Cont",
         cellwidth = 20, cellheight = 20,
         border_color = "black",
         fontsize_row = 15)
dev.off()

# day4
day4.mtx = make.cellGO.mtx('4day.vs.Ctrl')
tiff(filename = sprintf("%s/figure/Figure2D_day4_GCPC_CBC_upGO.tif", dir), width = 25, height = 5 + 0.75 * nrow(day4.mtx), units = 'cm', res = 300)
pheatmap(day4.mtx,
         cluster_rows = T,
         cluster_cols = F,
         color = colorRampPalette(c("white", "salmon", "red3"))(100),
         breaks = seq(0,4, length.out = 101),
         display_numbers = ifelse(day4.mtx > 1.5, "*", ""),
         number_color = "black",
         fontsize_number = 15,
         main = "Day4 vs Cont",
         cellwidth = 20, cellheight = 20,
         border_color = "black",
         fontsize_row = 15)
dev.off()







################
# Figure 2E
################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

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
#BiocManager::install("SiPSiC")
#BiocManager::install("progeny")
library(SiPSiC)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOfuncR)
library(reshape2)
library(ggpubr)


inpGS = get_anno_genes(c('GO:0031099'), database = "org.Mm.eg.db")
genes = unique(inpGS$gene)

ll = c("data.1", "data.2", "data.3")
i=1
for (i in 1:length(unique(s.integrated$sample))){
  sam = unique(s.integrated$sample)[i]
  cell.use = WhichCells(s.integrated, expression = cellType == "CBC" & sample == sam)
  DefaultAssay(s.integrated) = "RNA"
  c.mtx = GetAssayData(s.integrated, assay = "RNA", layer = ll[i])
  #c.mtx = log2(c.mtx+1)
  sce = SingleCellExperiment(assays = list(counts = as.matrix(c.mtx)))
  scores_and_idx = getPathwayScores(counts(sce), genes)
  pathway_scores = scores_and_idx$pathwayScores 
  gene_indices = scores_and_idx$index 
  s.integrated = AddMetaData(object = s.integrated, 
                             metadata = setNames(pathway_scores, cell.use), 
                             col.name = paste0("regenerationScore_", sam))
}

cbc.cells = WhichCells(s.integrated, expression = cellType == "CBC")
s.cbc = subset(s.integrated, cells = cbc.cells)
score.cols = grep("^regenerationScore_", colnames(s.cbc@meta.data), value = TRUE)

plot.df = FetchData(s.integrated, vars = c("sample", score.cols))
plot.long = melt(plot.df, id.vars = "sample", variable.name = "feature", value.name = "score", na.rm = T)
sample2 = plot.long[plot.long$sample %in% c('1day','4day'),]

grcol = c('lightslateblue','coral1')

tiff(filename = sprintf("%s/figure/Figure2F.tiff", dir), width = 4, height = 5, units = 'cm', res = 300)
ggplot(sample2, aes(x = sample, y = score, fill = sample)) +
  geom_boxplot(fill = grcol, alpha = 0.8, width = 0.8, size = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = sample), width = 0.1, size = 0.2, alpha = 0.5, stroke = 0.3) +
  scale_color_manual(values = grcol) +
  labs(x = "", y = "Regeneration Score") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major = element_line(color = "lightgrey", size = 0.1),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) 
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()



gcpc.cells = WhichCells(s.integrated, expression = cellType == "GC/PC")
s.gcpc = subset(s.integrated, cells = gcpc.cells)
score.cols = grep("^regenerationScore_", colnames(s.gcpc@meta.data), value = TRUE)

plot.df = FetchData(s.integrated, vars = c("sample", score.cols))
plot.long = melt(plot.df, id.vars = "sample", variable.name = "feature", value.name = "score", na.rm = T)
sample2 = plot.long[plot.long$sample %in% c('1day','4day'),]

grcol = c('lightslateblue','coral1')

tiff(filename = sprintf("%s/figure/Figure2F.tiff", dir), width = 4, height = 5, units = 'cm', res = 300)
ggplot(sample2, aes(x = sample, y = score, fill = sample)) +
  geom_boxplot(fill = grcol, alpha = 0.8, width = 0.8, size = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = sample), width = 0.1, size = 0.2, alpha = 0.5, stroke = 0.3) +
  scale_color_manual(values = grcol) +
  labs(x = "", y = "Regeneration Score") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major = element_line(color = "lightgrey", size = 0.1),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) 
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()

################
# Figure 2F
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/8-4_Trajectory-analysis_mergeGC-PC.Rdata", dir)) #cds

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(monocle)
library(gridExtra)
library(RColorBrewer)
library(GOfuncR)
library(clusterProfiler)
library(org.Mm.eg.db)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','lightslateblue','coral1')
scol = c('darkolivegreen2','tomato1','skyblue1','dodgerblue1','royalblue1')


tiff(filename = sprintf("%s/figure/supple_Figure2F-1.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right", legend.text = element_text(size = 16))
dev.off()

tiff(filename = sprintf("%s/figure/supple_Figure2F-1.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "State", show_branch_points = F) +
  scale_color_manual(values = scol) +
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
dev.off()


tiff(filename = sprintf("%s/figure/Figure2G-1.tiff", dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "cellType", show_branch_points = F) + 
  scale_color_manual(values = mycol) +
  theme(legend.position = "right", legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3))) 
dev.off()





################
# Figure 2G
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/8-4_Trajectory-analysis_mergeGC-PC.Rdata", dir)) #cds

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(monocle)
library(gridExtra)
library(RColorBrewer)
library(GOfuncR)
library(clusterProfiler)
library(org.Mm.eg.db)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
grcol = c('gold','lightslateblue','coral1')
scol = c('darkolivegreen2','tomato1','skyblue1','dodgerblue1','royalblue1')

stt = levels(unique(pData(cds)$sample))
stt.cc = pData(cds)[pData(cds)$sample %in% stt,]
stt.percent = stt.cc %>%
  group_by(State, cellType) %>%
  summarise(total_counts = n(), .groups = 'drop') %>%
  group_by(State) %>%
  mutate(percent = (total_counts / sum(total_counts)) * 100) %>%
  arrange(desc(percent), .by_group = TRUE)
stt.percent = expand.grid(
  State = unique(stt.percent$State),
  cellType = unique(stt.percent$cellType)
) %>%
  left_join(stt.percent, by = c("State", "cellType")) %>%
  mutate(percent = replace_na(percent, 0))


df = pData(cds)
celltype.lev = levels(factor(df$cellType))
col.name = setNames(mycol[1:length(celltype.lev)], celltype.lev)
df$cellType = factor(df$cellType, levels = names(col.name))

tiff(filename = sprintf("%s/figure/Figure2H_cont.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State == 3, sample == "Ctrl"),
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

tiff(filename = sprintf("%s/figure/Figure2H_day1.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State == 3, sample == "1day"),
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

tiff(filename = sprintf("%s/figure/Figure2H_day4.tiff", dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
ggplot(df %>% filter(State == 3, sample == "4day"),
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
