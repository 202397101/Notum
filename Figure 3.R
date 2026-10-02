################
# Figure 3A
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata", dir)) #markers.li

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
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}

cbc.1day.updeg = cluster.updeg$'1day'$CBC
cbc.1day.upora = enrichGO(gene = rownames(cbc.1day.updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
cbc.1day.upgofilter = as.data.frame(gofilter(cbc.1day.upora, level = 5))
cbc.1day.upgofilter = cbc.1day.upgofilter[cbc.1day.upgofilter$p.adjust < 0.05, ]
cbc.1day.upgofilter$p.adjust = -log10(cbc.1day.upgofilter$p.adjust)
cbc.1day.upgofilter = cbc.1day.upgofilter[order(cbc.1day.upgofilter$p.adjust, decreasing = TRUE), ]
cbc.1day.upgofilter$Description = str_to_title(cbc.1day.upgofilter$Description)

ggplot(data = cbc.1day.upgofilter, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = -log10(0.05), color = "red", size = 0.5) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = "lightpink", mid = "lightpink", high = "lightpink", guide = FALSE, midpoint = 2) +
  ggtitle('Day1: CBC vs all others') +
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



################
# Figure 3B
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata", dir)) #markers.li

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
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}

cbc.4day.updeg = cluster.updeg$'4day'$CBC
cbc.4day.upora = enrichGO(gene = rownames(cbc.4day.updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
cbc.4day.upgofilter = as.data.frame(gofilter(cbc.4day.upora, level = 5))
cbc.4day.upgofilter = cbc.4day.upgofilter[cbc.4day.upgofilter$p.adjust < 0.05, ]
cbc.4day.upgofilter$p.adjust = -log10(cbc.4day.upgofilter$p.adjust)
cbc.4day.upgofilter = cbc.4day.upgofilter[order(cbc.4day.upgofilter$p.adjust, decreasing = TRUE), ]
cbc.4day.upgofilter = head(cbc.4day.upgofilter, 15)
cbc.4day.upgofilter$Description = str_to_title(cbc.4day.upgofilter$Description)

ggplot(data = cbc.4day.upgofilter, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = -log10(0.05), color = "red", size = 0.5) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = "lightpink", mid = "lightpink", high = "lightpink", guide = FALSE, midpoint = 2) +
  ggtitle('Day4: CBC vs all others') +
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
    plot.margin = margin(t = 0.2, r = 0, b = 0.3, l = -0.5, unit = "cm")
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", size = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", size = 0.3) +
  geom_text(aes(label = Description, y = 0.05), hjust = 0, size = 4, color = "black")+
  scale_y_continuous(expand = c(0, 0)) +  
  scale_x_discrete(expand = c(0, 1))    

neg.wnt = cbc.4day.upgofilter$geneID[cbc.4day.upgofilter$ID %in% 'GO:0030178']
neg.wnt = unlist(strsplit(neg.wnt, split = "/"))
#save(neg.wnt, file = sprintf("%s/Figure3B_neg-wnt-genes.Rdata", dir))


################
# Figure 3C
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata", dir)) #markers.li

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
library(ComplexHeatmap)
library(homologene)

cluster.updeg = list()
i=1
for (i in 1:length(markers.li)){
  deg = markers.li[[i]]  
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}

cbc.cont.updeg = cluster.updeg$ctrl$CBC
cbc.1day.updeg = cluster.updeg$'1day'$CBC
cbc.4day.updeg = cluster.updeg$'4day'$CBC

cbc.updeg = list(cont = cbc.cont.updeg, day1 = cbc.1day.updeg, day4 = cbc.4day.updeg)

i=1
wnt.ORA = list()
for (i in 1:length(cbc.updeg)){
  updeg = rownames(cbc.updeg[[i]])
  ora = enrichGO(gene = updeg, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
  wnt = ora@result[grep('Wnt', ora@result$Description), ]
  wnt.ora = new("enrichResult",
                result = wnt,           
                pvalueCutoff = ora@pvalueCutoff,
                pAdjustMethod = ora@pAdjustMethod,
                organism = ora@organism,
                ontology = ora@ontology,
                gene = ora@gene,
                universe = ora@universe,
                gene2Symbol = ora@gene2Symbol,
                keytype = ora@keytype,
                readable = ora@readable)
  wnt.ORA[[names(cbc.updeg)[i]]] = wnt.ora
}

des = lapply(wnt.ORA, function(i) {
  i@result$Description
})
inter.des = Reduce(intersect, des)
inter.des = inter.des[c(5,3,9,1)]

i=1
wnt.upORA = list()
for (i in 1:length(wnt.ORA)){
  ora = wnt.ORA[[i]]
  wnt = ora@result[ora@result$Description %in% inter.des, ]
  wnt = wnt[order(wnt$Description),]
  wnt.ora = new("enrichResult",
                result = wnt,           
                pvalueCutoff = ora@pvalueCutoff,
                pAdjustMethod = ora@pAdjustMethod,
                organism = ora@organism,
                ontology = ora@ontology,
                gene = ora@gene,
                universe = ora@universe,
                gene2Symbol = ora@gene2Symbol,
                keytype = ora@keytype,
                readable = ora@readable)
  wnt.upORA[[names(wnt.ORA)[i]]] = wnt.ora
}
lapply(wnt.upORA, function(i) i@result)


merged.data = do.call(rbind, lapply(names(wnt.upORA), function(i) {
  if (!is.null(wnt.upORA[[i]])) {
    data = wnt.upORA[[i]]@result
    data$Sample = i 
    return(data)
  }
}))

merged.data$p.adjust = -log10(merged.data$p.adjust)
merged.data$Sample = factor(merged.data$Sample, levels = unique(merged.data$Sample), labels = c("Cont", "Day1", "Day4"))
range(merged.data$p.adjust)

tiff(filename = sprintf("%s/figure/Figure3C.tiff", dir), width = 13, height = 5, units = 'cm', res = 300)
ggplot(merged.data, aes(x = Sample, y = Description, size = p.adjust, color = Count)) +
  geom_point()+
  scale_color_gradient2(low = "lightgrey", mid = "red", high = "red2", midpoint = 15)+ #scale_colour_gradient2(low = "dimgrey", mid = "grey", high = "red1", midpoint = 2)+
  scale_size_continuous(
    range = c(-2, 5), 
    name = "-logFDR"
  ) +
  theme_minimal() +
  labs(x = "", y = "", color = "Count", size = "-logFDR") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
dev.off()



################
# Figure 3D
################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers

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
library(gridExtra)
library(homologene)
library(fgsea)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(SiPSiC)
library(SingleCellExperiment)


gsl = gmtPathways("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt") # MsigDB 에서 다운받은 GMT 파일.

wnt = gsl[grep('PID_WNT_CANONICAL_PATHWAY', names(gsl))]
wnt = unique(unlist(wnt))
wnt = homologene(wnt, inTax = 9606, outTax = 10090)
wnt = unique(wnt$'10090')

ll = c("data.1", "data.2", "data.3")
i=1
for (i in 1:length(unique(s.integrated$sample))){
  sam = unique(s.integrated$sample)[i]
  cell.use = WhichCells(s.integrated, expression = cellType == "CBC" & sample == sam)
  DefaultAssay(s.integrated) = "RNA"
  c.mtx = GetAssayData(s.integrated, assay = "RNA", layer = ll[i])
  #c.mtx = log2(c.mtx+1)
  sce = SingleCellExperiment(assays = list(counts = as.matrix(c.mtx)))
  scores_and_idx = getPathwayScores(counts(sce), wnt)
  pathway_scores = scores_and_idx$pathwayScores 
  gene_indices = scores_and_idx$index 
  s.integrated = AddMetaData(object = s.integrated, 
                             metadata = setNames(pathway_scores, cell.use), 
                             col.name = paste0("WNTscore_", sam))
}

cbc.cells = WhichCells(s.integrated, expression = cellType == "CBC")
s.cbc = subset(s.integrated, cells = cbc.cells)
score.cols = grep("^WNTscore_", colnames(s.cbc@meta.data), value = TRUE)

plot.df = FetchData(s.integrated, vars = c("sample", score.cols))
plot.long = melt(plot.df, id.vars = "sample", variable.name = "feature", value.name = "score", na.rm = T)
sample2 = plot.long[plot.long$sample %in% c('1day','4day'),]

grcol = c('lightslateblue','coral1')

tiff(filename = sprintf("%s/figure/Figure3D.tiff", dir), width = 4, height = 5, units = 'cm', res = 300)
ggplot(sample2, aes(x = sample, y = score, fill = sample)) +
  geom_boxplot(fill = grcol, alpha = 0.8, width = 0.8, size = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = sample), width = 0.1, size = 0.2, alpha = 0.5, stroke = 0.3) +
  scale_color_manual(values = grcol) +
  labs(x = "", y = "Wnt Score") +
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
# Figure 3H
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/GEO"

load(file = sprintf("%s/Rdata/11-1_Notum.Rdata", dir)) #mean.matrix, layer5
load(file = sprintf("%s/GEO/merge/Rdata/8_merge_tpm.Rdata", dir)) #epi1.tpm, epi2.tpm, paneth.tpm, sm6Stem.tpm, tpm.li
load(file = sprintf("%s/GEO/merge/Rdata/8_merge_DEG.Rdata", dir)) #epi1.deg, epi2.deg, paneth.deg, sm6Stem.deg, deg.li

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(ggrepel)
library(stringr)
library(GOfuncR)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(gridExtra)
library(ggvenn)



genes = colnames(mean.matrix)[mean.matrix>0]

edegl = deg.li[c(1,2,3)]
resl = lapply(edegl, function(m) {
  mm = m[match(genes, m$Genes),]
  list(logFC = mm[,colnames(mm) %in% c("log2FoldChange")], Pval = mm[,colnames(mm) %in% c("pvalue")])
})

logfc = do.call("cbind", lapply(resl, function(l) l$logFC))
pval = do.call("cbind", lapply(resl, function(l) l$Pval))
rownames(logfc) = genes
rownames(pval) = genes
logfc[is.na(logfc)] = 0
pval[is.na(pval)] = 1

sn = apply(is.na(logfc), 1, function(v) sum(v, na.rm=T)) 
use_idx = sn!=ncol(logfc)

logfc1 = logfc[use_idx,]
pval1 = pval[use_idx,]

asterisk_mat = ifelse(pval1 < 0.05, "*", "")
asterisk_mat_t = t(asterisk_mat)

color.ht = colorRamp2(c(-2,-0.02,0,0.02,2), c('blue4','blue','white','red','red4'))
lg.ht = Legend(title="Log2FC", at=c(-2,0,2), col_fun=color.ht, border='black', title_position="topcenter")
hm = Heatmap(t(logfc1), show_heatmap_legend = F, heatmap_legend_param = list(title = 'log2FC'), 
             col = color.ht, border = T, cluster_columns = T, cluster_rows = T,
             show_column_names = T, cluster_column_slices = T, column_names_rot = -90, column_names_gp = gpar(fontsize = 11),
             row_names_gp = gpar(fontsize = 10, fontface = "bold"),
             width  = ncol(t(logfc1)) * unit(5, "mm"),
             height = nrow(t(logfc1)) * unit(5, "mm"),
             rect_gp = gpar(col = "white", lwd = 1),
             cell_fun = function(j, i, x, y, width, height, fill) {
               if (asterisk_mat_t[i, j] == "*") {
                 grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold", col = "black"))
               }
             })

tiff(filename = sprintf("%s/figure/Figure4A.tiff", dir), width = 10, height = 5, units = "cm", res = 300)
draw(hm, annotation_legend_list = list(lg.ht))
dev.off()
