#####################################################################
# Figure 2 - scRNA-seq of PGE2-treated mouse intestinal epithelium (Ctrl, Day1, Day4)
#
# Code for: "Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration"
# GitHub : https://github.com/202397101/Notum
# Zenodo : https://doi.org/10.5281/zenodo.23102321
#
# Each panel section below can be run on its own (it loads its own input and packages).
# See README.md for the input files and the folder layout.
#####################################################################

# ---- Paths (edit before running) ----
dir     = "data"                          # input data folder (Rdata/, GEO/, MsigDB/)
out_dir = "figure"                        # output TIFF files
dir.create(out_dir, showWarnings = FALSE)



################
# Figure 2A
################
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

library(Seurat)
library(tidyverse)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")

tiff(filename = sprintf("%s/Figure2A.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
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
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

library(Seurat)
library(tidyverse)
library(ggrepel)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")

# cell type proportion per sample
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
cell.pop$cellType = factor(cell.pop$cellType, levels = rev(unique(cell.pop$cellType)))

# percentage labels (1 decimal)
cell.pop = cell.pop %>% ungroup() %>% mutate(label = sprintf("%.1f%%", percent))
as.data.frame(tidyr::pivot_wider(cell.pop, id_cols = cellType, names_from = sample, values_from = label))

# label y position = middle of each stacked segment (same order as position_stack: last level at the bottom)
cell.pop = cell.pop %>%
  group_by(group) %>%
  arrange(desc(cellType), .by_group = TRUE) %>%
  mutate(ymax = cumsum(percent), ymid = ymax - percent / 2) %>%
  ungroup()

small.cut = 3  # segments < 3% are labelled outside the bar with a leader line
bar.w = 0.5
dark.ct = c("TA", "EC", "TC")  # white text on dark segments

ggplot(cell.pop, aes(x = group, y = percent, fill = cellType, width = bar.w)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(data = subset(cell.pop, percent >= small.cut),
            aes(x = group, y = ymid, label = label, color = ifelse(cellType %in% dark.ct, "white", "black")),
            inherit.aes = FALSE, size = 2.8) +
  geom_text_repel(data = subset(cell.pop, percent < small.cut),
                  aes(x = as.numeric(group) + bar.w / 2, y = ymid, label = paste0(cellType, " ", label)), inherit.aes = FALSE,
                  size = 2.3, nudge_x = 0.1, direction = "y", hjust = 0, color = "black",
                  segment.size = 0.3, min.segment.length = 0, box.padding = 0.1, point.padding = 0) +
  scale_color_identity() +
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
  scale_x_discrete(expand = expansion(add = c(0.5, 0.8))) +
  coord_cartesian(clip = "off")



################
# Figure 2C
################
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers
load(file = sprintf("%s/Rdata/11-2_Augur.Rdata", dir)) #augur.day4, augur.day1, augur.day1.day4 (Augur::calculate_auc results)

library(Seurat)
library(tidyverse)
library(ggrepel)


celltype.lev = levels(factor(s.integrated$cellType))
col.name = setNames(ifelse(celltype.lev %in% c("GC/PC", "CBC"), "red", "black"), celltype.lev)  # highlight GC/PC and CBC
col.name

# Augur AUC per cell type: Day1 vs Cont, Day4 vs Cont
augur.df = rbind(
  data.frame(augur.day1$AUC, comparison = "Day1 vs Cont"),
  data.frame(augur.day4$AUC, comparison = "Day4 vs Cont"))
augur.df$comparison = factor(augur.df$comparison, levels = c("Day1 vs Cont", "Day4 vs Cont"))
augur.df$cell_type = factor(augur.df$cell_type, levels = celltype.lev)
augur.df

tiff(filename = sprintf("%s/Figure2C.tiff", out_dir), width = 9, height = 9, units = 'cm', res = 300)
ggplot(augur.df, aes(x = comparison, y = auc, group = cell_type, color = cell_type)) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed", linewidth = 0.3) +  # AUC 0.5 = no separation
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 2) +
  geom_text_repel(data = augur.df %>% filter(comparison == "Day4 vs Cont"),
                  aes(label = cell_type), nudge_x = 0.3, direction = "y", hjust = 0,
                  size = 5, segment.size = 0.1, show.legend = FALSE) +
  scale_color_manual(values = col.name) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
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
# GO BP (up DEGs, Day1/Day4 vs Cont) of CBC and GC/PC
go.env = new.env()
load(file = sprintf("%s/Rdata/6-4_CBC_celltype_pathway_heatmap.Rdata", dir), envir = go.env) #upGO (enrichGO result per cell type and comparison)
cellGO = go.env$upGO
names(cellGO)

library(clusterProfiler)
library(pheatmap)

sel.cell = c("CBC", "GC/PC")  # heatmap column order
cellGO = cellGO[intersect(sel.cell, names(cellGO))]
names(cellGO)

# top GO terms (level 4) of each cell type for one comparison -> -log10 p matrix (1 = not enriched)
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
tiff(filename = sprintf("%s/Figure2D_day1.tiff", out_dir), width = 25, height = 5 + 0.75 * nrow(day1.mtx), units = 'cm', res = 300)
pheatmap(day1.mtx,
         cluster_rows = T,
         cluster_cols = F,
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
tiff(filename = sprintf("%s/Figure2D_day4.tiff", out_dir), width = 25, height = 5 + 0.75 * nrow(day4.mtx), units = 'cm', res = 300)
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
# Regeneration score (SiPSiC, GO:0031099 regeneration) of CBC, Day1 vs Day4
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir)) # s.integrated, markers

library(Seurat)
library(tidyverse)
library(SiPSiC)
library(SingleCellExperiment)
library(GOfuncR)
library(reshape2)


inpGS = get_anno_genes(c('GO:0031099'), database = "org.Mm.eg.db")
genes = unique(inpGS$gene)

# SiPSiC per sample (RNA layers data.1/2/3 = samples, in the order of unique(s.integrated$sample))
ll = c("data.1", "data.2", "data.3")
for (i in 1:length(unique(s.integrated$sample))){
  sam = unique(s.integrated$sample)[i]
  cell.use = WhichCells(s.integrated, expression = cellType == "CBC" & sample == sam)
  DefaultAssay(s.integrated) = "RNA"
  c.mtx = GetAssayData(s.integrated, assay = "RNA", layer = ll[i])
  sce = SingleCellExperiment(assays = list(counts = as.matrix(c.mtx)))
  scores_and_idx = getPathwayScores(counts(sce), genes)
  pathway_scores = scores_and_idx$pathwayScores
  gene_indices = scores_and_idx$index
  s.integrated = AddMetaData(object = s.integrated,
                             metadata = setNames(pathway_scores, cell.use),
                             col.name = paste0("regenerationScore_", sam))
}

score.cols = grep("^regenerationScore_", colnames(s.integrated@meta.data), value = TRUE)

plot.df = FetchData(s.integrated, vars = c("sample", score.cols))
plot.long = melt(plot.df, id.vars = "sample", variable.name = "feature", value.name = "score", na.rm = T)
sample2 = plot.long[plot.long$sample %in% c('1day','4day'),]

grcol = c('lightslateblue','coral1')

tiff(filename = sprintf("%s/Figure2E.tiff", out_dir), width = 4, height = 5, units = 'cm', res = 300)
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
dev.off()



################
# Figure 2F
################
# Monocle2 (DDRTree) trajectory
load(file = sprintf("%s/Rdata/8-4_Trajectory-analysis_mergeGC-PC.Rdata", dir)) #cds

library(tidyverse)
library(monocle)
library(RColorBrewer)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")
scol = c('darkolivegreen2','tomato1','skyblue1','dodgerblue1','royalblue1')


tiff(filename = sprintf("%s/Figure2F_pseudotime.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right", legend.text = element_text(size = 16))
dev.off()

tiff(filename = sprintf("%s/Figure2F_State.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "State", show_branch_points = F) +
  scale_color_manual(values = scol) +
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

tiff(filename = sprintf("%s/Figure2F_cellType.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "cellType", show_branch_points = F) +
  scale_color_manual(values = mycol) +
  theme(legend.position = "right", legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()




################
# Figure 2G
################
# Pseudotime of State 3 cells in each sample
load(file = sprintf("%s/Rdata/8-4_Trajectory-analysis_mergeGC-PC.Rdata", dir)) #cds

library(tidyverse)
library(monocle)


mycol = c("greenyellow","darkgreen","orchid1","cyan","orangered1","blue","orange","red4")

df = pData(cds)
celltype.lev = levels(factor(df$cellType))
col.name = setNames(mycol[1:length(celltype.lev)], celltype.lev)
df$cellType = factor(df$cellType, levels = names(col.name))

tiff(filename = sprintf("%s/Figure2G_cont.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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

tiff(filename = sprintf("%s/Figure2G_day1.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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

tiff(filename = sprintf("%s/Figure2G_day4.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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
