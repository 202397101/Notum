###############
# Figure 6A
###############
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

#load(file = file.path(dir_geo,"4_GSE174297_organoid","Rdata","1_seurat-object.Rdata")) #sobj
load(file = file.path(dir_geo,"4_GSE174297_organoid","Rdata","2_sobj.Rdata")) #sobj

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(homologene)
library(fgsea)
library(ggpubr)
library(monocle)
library(RColorBrewer)
library(reshape2)
library(SiPSiC)

sobj@images = list()


mycol = c('CBC' = '#7FE5D1',
          'TA' = '#87C1DD',
          'EC progenitor' = '#7F8DD9',
          'EC' = '#788785',
          'EE' = '#D6A26C',
          'GC/PC' = '#D3DF52',
          'TC' = '#8E4CE0')

sobj$cellType = factor(sobj$cellType, levels = names(mycol))
sobj$treatment = factor(sobj$treatment, levels = c('sc_control','sc_ifnr'))
grcol = c('grey','darkorange') # sc_control, sc_ifnr

tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/2_UMAP_cellType.tiff", dir_geo), width = 20, height = 10, units = 'cm', res = 300)
DimPlot(sobj, reduction = 'umap', group.by = 'cellType', split.by = 'treatment', cols = mycol, repel = F, label = F)
dev.off()




###############
# Figure 6B
###############
cell.pop = sobj@meta.data %>%
  dplyr::count(treatment, cellType, name = "counts") %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(percent = counts / sum(counts) * 100)
cell.pop$cellType = factor(cell.pop$cellType, levels = names(mycol))

tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/2_IFNr_cellPopulation_barplot.tiff", dir_geo), width = 12, height = 10, units = 'cm', res = 300)
ggplot(cell.pop, aes(x = treatment, y = percent, fill = cellType, width = 0.7)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  scale_fill_manual(values = mycol) +
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
dev.off()



###############
# Figure 6C
###############

DefaultAssay(sobj) = "RNA"
tc.markers = c('Dclk1','Trpm5')
tc.cells = WhichCells(sobj, expression = cellType == 'TC') # TC 세포만

tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/2_FeaturePlot_Dclk1-Trpm5_TC.tiff", dir_geo), width = 16, height = 14, units = 'cm', res = 300)
FeaturePlot(sobj, features = tc.markers, cells = tc.cells, order = TRUE, split.by = 'treatment', min.cutoff = 0, by.col = T, max.cutoff = 'q90') &
  theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 8),
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank())
dev.off()



###############
# Figure 6D, E
###############

dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

load(file = file.path(dir_geo,"4_GSE174297_organoid","Rdata","2_sobj.Rdata")) #sobj
load(file = file.path(dir_geo,"4_GSE174297_organoid","Rdata","3_Trajectory-analysis.Rdata")) #cds


library(dplyr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(monocle)
library(RColorBrewer)

# monocle2 plot_cell_trajectory() 가 쓰는 dplyr::select_() 가 dplyr 1.2 에서 없어짐 -> select() 로 연결 (현재 R 세션에서만 적용)
select_.fix = function(.data, ...) dplyr::select(.data, ...)
assignInNamespace("select_", select_.fix, ns = "dplyr")
mono.imp = parent.env(asNamespace("monocle"))
if (exists("select_", envir = mono.imp, inherits = FALSE)) {
  unlockBinding("select_", mono.imp)
  assign("select_", select_.fix, envir = mono.imp)
  lockBinding("select_", mono.imp)
}


n.state = length(unique(pData(cds)$State))
scol = colorRampPalette(c('skyblue1','darkolivegreen2','tomato1','gold','orchid1'))(n.state)
mycol = c('CBC' = '#7FE5D1',
          'TA' = '#87C1DD',
          'EC progenitor' = '#7F8DD9',
          'EC' = '#788785',
          'EE' = '#D6A26C',
          'GC/PC' = '#D3DF52',
          'TC' = '#8E4CE0')
sobj$cellType = factor(sobj$cellType, levels = names(mycol))
sobj$treatment = factor(sobj$treatment, levels = c('sc_control','sc_ifnr'))
grcol = c('grey','darkorange') # sc_control, sc_ifnr


tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/3_trajectory_pseudotime.tiff", dir_geo), width = 13, height = 12, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(legend.position = "right", legend.text = element_text(size = 16))
dev.off()


tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/3_trajectory_cellType_split.tiff", dir_geo), width = 24, height = 12, units = 'cm', res = 300)
plot_cell_trajectory(cds, color_by = "cellType", show_branch_points = F) +
  scale_color_manual(values = mycol) +
  theme(legend.position = "right", legend.text = element_text(size = 16), strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~treatment, nrow = 1)
dev.off()


plot_cell_trajectory(cds, color_by = "State", show_branch_points = F)


###############
# Figure 6F
###############
df = pData(cds) # 세포별 State, treatment, cellType
df$treatment = factor(df$treatment, levels = c('sc_control','sc_ifnr'))

pop.state2 = df %>%
  filter(State == 1) %>%
  dplyr::count(treatment, cellType, name = "counts") %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(percent = counts / sum(counts) * 100)
pop.state2$cellType = factor(pop.state2$cellType, levels = names(mycol))

ggplot(pop.state2, aes(x = treatment, y = percent, fill = cellType, width = 0.7)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  scale_fill_manual(values = mycol) +
  labs(y = "% of cell population (State 2)") +
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



###############
# Figure 6G
###############
sobj@images = list()
DefaultAssay(sobj) = "RNA"

gsl = gmtPathways("C:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt") # MsigDB 에서 다운받은 GMT 파일.
names(gsl)[grep('WNT', names(gsl))]
wnt = gsl[grep('PID_WNT_CANONICAL_PATHWAY', names(gsl))]
wnt = unique(unlist(wnt))
wnt = homologene(wnt, inTax = 9606, outTax = 10090)
wnt = unique(wnt$'10090')

sobj = AddModuleScore(sobj, features = list(wnt), name = "Wnt_Score")
sobj = UpdateSeuratObject(sobj)
sobj.sub = sobj[, sobj$Cd74_group != "grey_Cd74"]

sobj.cbcta = sobj[, sobj$treatment == "sc_ifnr" & sobj$Cd74_group %in% c("negative_Cd74","positive_Cd74") & sobj$cellType %in% c("CBC","TA")]
wnt.df = FetchData(sobj.cbcta, vars = c("Wnt_Score1", "Cd74_group", "cellType"))
wnt.df$cellType = droplevels(factor(wnt.df$cellType, levels = c("CBC","TA")))

########## Control vs IFNr 에서 Cd74+/- 세포 비율 ##########
# grey_Cd74 (애매한 세포) 는 위 Wnt score 그림과 같이 제외
cd74col = c('negative_Cd74' = 'grey80', 'positive_Cd74' = 'darkorange')

cd74.pop = sobj@meta.data %>%
  filter(Cd74_group != "grey_Cd74") %>%
  dplyr::count(treatment, Cd74_group, name = "counts") %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(percent = counts / sum(counts) * 100)
cd74.pop$treatment = factor(cd74.pop$treatment, levels = c('sc_control','sc_ifnr'))
cd74.pop$Cd74_group = factor(cd74.pop$Cd74_group, levels = c('positive_Cd74','negative_Cd74'))
cd74.pop

tiff(filename = sprintf("%s/figure/supple_figure6G_cell-pop.tiff", dir), width = 8, height = 6, units = 'cm', res = 300)
ggplot(cd74.pop, aes(x = treatment, y = percent, fill = Cd74_group, width = 0.7)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  #geom_text(aes(label = sprintf("%.1f%%", percent)), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = cd74col, labels = c('positive_Cd74' = 'Cd74+', 'negative_Cd74' = 'Cd74-'), name = "") +
  scale_x_discrete(labels = c('sc_control' = 'Control', 'sc_ifnr' = 'IFNr'), expand = c(0, 0.5)) +
  labs(y = "% of cells") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") +
  geom_vline(xintercept = 0.5, color = "black", size = 0.5, linetype = "solid") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        legend.text = element_text(size = 11)) +
  scale_y_continuous(expand = c(0, 0))
dev.off()


tiff(filename = sprintf("%s/4_GSE174297_organoid/figure/2_wnt_CBC-TA_pooled.tiff", dir_geo), width = 7, height = 10, units = 'cm', res = 300)
ggplot(wnt.df, aes(x = Cd74_group, y = Wnt_Score1, fill = Cd74_group)) +
  geom_violin(scale = "width", trim = T) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = c("peachpuff", "darkorange")) +
  #stat_compare_means(method = "wilcox.test", na.rm = TRUE) +
  labs(x = "", y = "Wnt Score") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 10),
        axis.line = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"))
dev.off()




