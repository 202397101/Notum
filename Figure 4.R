#####################################################################
# Figure 4 - scRNA-seq of mouse intestinal epithelium treated with PGE2 and a NOTUM inhibitor (ABC99 = iNOTUM)
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



#######################
# Figure 4A
#######################
load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library(Seurat)
library(tidyverse)


DefaultAssay(s.integrated) = "RNA"

mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")

tiff(filename = sprintf("%s/Figure4A.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
DimPlot(s.integrated, group.by = 'cellType', label=F, cols = mycol) + ggtitle('')+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
dev.off()



#######################
# Figure 4B
#######################
load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library(Seurat)
library(tidyverse)
library(ggrepel)


mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")

# cell type proportion per sample (Cont, PGE2, ABC99)
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
cell.pop$cellType = factor(cell.pop$cellType, levels = rev(unique(cell.pop$cellType)))

# percentage labels (1 decimal), same as Figure 2B
cell.pop = cell.pop %>% ungroup() %>% mutate(label = sprintf("%.1f%%", percent))
as.data.frame(tidyr::pivot_wider(cell.pop, id_cols = cellType, names_from = sample, values_from = label))

# label y position = middle of each stacked segment (same order as position_stack: last level at the bottom)
cell.pop = cell.pop %>%
  group_by(group) %>%
  arrange(desc(cellType), .by_group = TRUE) %>%
  mutate(ymax = cumsum(percent), ymid = ymax - percent / 2) %>%
  ungroup()

small.cut = 4  # segments < 4% are labelled outside the bar with a leader line
bar.w = 0.45
dark.ct = c("EC", "TC")  # white text on dark segments

ggplot(cell.pop, aes(x = group, y = percent, fill = cellType, width = bar.w)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(data = subset(cell.pop, percent >= small.cut),
            aes(x = group, y = ymid, label = label, color = ifelse(cellType %in% dark.ct, "white", "black")),
            inherit.aes = FALSE, size = 3.4, fontface = "bold") +
  geom_text_repel(data = subset(cell.pop, percent < small.cut),
                  aes(x = as.numeric(group) + bar.w / 2, y = ymid, label = paste0(cellType, " ", label)), inherit.aes = FALSE,
                  size = 2.7, nudge_x = 0.05, direction = "y", hjust = 0, color = "black",
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
  scale_x_discrete(expand = expansion(add = c(0.5, 0.9))) +
  coord_cartesian(clip = "off")




#######################
# Figure 4C
#######################
load(file = sprintf("%s/Rdata/13-2_5samples_integration.Rdata", dir)) #s.integrated, markers

library(Seurat)
library(tidyverse)


DefaultAssay(s.integrated) = "RNA"

# FeaturePlot: Lyz1 (Paneth cell), Muc2 (goblet cell) per sample
tiff(filename = sprintf("%s/Figure4C_FeaturePlot_Lyz1.tiff", out_dir), width = 12, height = 25, units = 'cm', res = 300)
FeaturePlot(s.integrated, features = 'Lyz1', order = TRUE, split.by = 'sample', min.cutoff = 0, max.cutoff = 5, by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank())
dev.off()

tiff(filename = sprintf("%s/Figure4C_FeaturePlot_Muc2.tiff", out_dir), width = 12, height = 25, units = 'cm', res = 300)
FeaturePlot(s.integrated, features = 'Muc2', order = TRUE, split.by = 'sample', min.cutoff = 0, max.cutoff = 3.5, by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank())
dev.off()


########## Lyz1, Muc2 dot plot per cell type (Cont, PGE2, iNOTUM) ##########
#  - y = cell type, x = sample; left = Lyz1, right = Muc2
#  - dot size = % expressing cells, colour = mean expression (RNA assay, log-normalized)
sam.use = c("Cont", "PGE2", "ABC99")
sam.lab = c("Cont" = "Cont", "PGE2" = "PGE2", "ABC99" = "iNOTUM")
cell.lev = levels(factor(s.integrated$cellType))

dot.c.df = FetchData(s.integrated, vars = c("sample", "cellType", "Lyz1", "Muc2"), layer = "data")
dot.c.df = dot.c.df[dot.c.df$sample %in% sam.use, ]

dot.c = dot.c.df %>%
  tidyr::pivot_longer(cols = c("Lyz1", "Muc2"), names_to = "gene", values_to = "expr") %>%
  group_by(gene, cellType, sample) %>%
  summarise(n_cells = n(), mean = mean(expr), pct = mean(expr > 0) * 100, .groups = "drop") %>%
  mutate(sample = factor(sam.lab[as.character(sample)], levels = sam.lab),
         cellType = factor(cellType, levels = rev(cell.lev)),
         gene = factor(gene, levels = c("Lyz1", "Muc2"))) %>%
  as.data.frame()
dot.c

tiff(filename = sprintf("%s/Figure4C_dotplot.tiff", out_dir), width = 9, height = 7, units = "cm", res = 300)
ggplot(dot.c, aes(x = sample, y = cellType)) +
  geom_point(aes(size = pct, fill = mean), shape = 21, color = "dimgrey", stroke = 0.3) +
  facet_grid(. ~ gene) +
  scale_fill_gradient(low = "white", high = "red3", name = "Mean expr")+
  scale_size_continuous(range = c(0.3, 4.5), limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "% expressing") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold.italic"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "cm"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6),
        legend.key.size = unit(0.35, "cm"))
dev.off()




#######################
# Figure 4D
#######################
# Monocle2 (DDRTree) trajectory
load(file = sprintf("%s/Rdata/13-2_Trajectory-analysis.Rdata", dir)) #cds

library(tidyverse)
library(monocle)
library(RColorBrewer)


mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")
scol = c('skyblue1','darkolivegreen2','tomato1')

tiff(filename = sprintf("%s/Figure4D_pseudotime.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right", legend.text = element_text(size = 16))
dev.off()

tiff(filename = sprintf("%s/Figure4D_State.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "State", show_branch_points = F) +
  scale_color_manual(values = scol) +
  theme(legend.position = "right", legend.text = element_text(size = 18), legend.title = element_text(size = 18, face = "bold"), strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

tiff(filename = sprintf("%s/Figure4D_cellType.tiff", out_dir), width = 20, height = 20, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "cellType", show_branch_points = F) +
  scale_color_manual(values = mycol) +
  theme(legend.position = "right", legend.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()




#######################
# Figure 4E
#######################
# Pseudotime of State 3 cells in each sample
load(file = sprintf("%s/Rdata/13-2_Trajectory-analysis.Rdata", dir)) #cds

library(tidyverse)
library(monocle)


mycol = c("greenyellow","orchid1","cyan","orangered1","blue","orange","red4")

df = pData(cds)
celltype.lev = levels(factor(df$cellType))
col.name = setNames(mycol[1:length(celltype.lev)], celltype.lev)
df$cellType = factor(df$cellType, levels = names(col.name))

tiff(filename = sprintf("%s/Figure4E_cont.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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

tiff(filename = sprintf("%s/Figure4E_pge2.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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

tiff(filename = sprintf("%s/Figure4E_inhibitor.tiff", out_dir), width = 5.8, height = 2.5, units = 'cm', res = 300)
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
# Figure 4J-1
#######################
# GO BP ORA of CBC up DEGs, ABC99 vs PGE2
load(file = sprintf("%s/Rdata/13-1_ABC99_vs_PGE2_DEG.Rdata", dir)) #cluster.degl (Seurat FindMarkers per cell type)

library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)



cluster.updeg = list()
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

# terms shown: top 58 -> remove metabolic/catabolic/"regulation of" terms and redundant terms
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
# Figure 4J-2
#######################
# GO BP ORA of ABC99 (iNOTUM) CBC up DEGs (CBC vs all other cells)
load(file = sprintf("%s/Rdata/13-1_CBC_vs_others_DEG.Rdata", dir)) #markers.li (Seurat FindMarkers, each cell type vs others, per sample)


library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)



cluster.updeg = list()
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
# terms shown: remove metabolic/catabolic/"regulation of" terms and redundant terms -> top 10
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
