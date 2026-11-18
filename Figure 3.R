#####################################################################
# Figure 3 - Wnt signaling and Notum
#
# Code for: "Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration"
# GitHub : https://github.com/202397101/Notum
# Zenodo : https://doi.org/10.5281/zenodo.23102321
#
# Data   : PGE2-treated mouse epithelium scRNA-seq (this study), GSE84061 (mouse ISC/Paneth bulk RNA-seq),
#          E-MTAB-9543 (human colon scRNA-seq)
# Each panel section below can be run on its own (it loads its own input and packages).
# See README.md for the input files and the folder layout.
#####################################################################

# ---- Paths (edit before running) ----
dir     = "data"                          # input data folder (Rdata/, GEO/, MsigDB/)
dir_geo = file.path(dir, "public_scRNA")  # public scRNA-seq datasets
dir_db  = file.path(dir, "DB")            # msigdb/ (MSigDB GMT), GTEx/
out_dir = "figure"                        # output TIFF files
dir.create(out_dir, showWarnings = FALSE)



################
# Figure 3A left
################
# GO BP ORA of Day1 CBC up DEGs (CBC vs all other cells)
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata", dir)) #markers.li (Seurat FindMarkers, each cell type vs others, per sample)

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)



cluster.updeg = list()
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
# Figure 3A right
################
# GO BP ORA of Day4 CBC up DEGs (CBC vs all other cells); run after "Figure 3A left" (uses cluster.updeg)
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

# genes of "negative regulation of Wnt signaling pathway" (GO:0030178)
neg.wnt = cbc.4day.upgofilter$geneID[cbc.4day.upgofilter$ID %in% 'GO:0030178']
neg.wnt = unlist(strsplit(neg.wnt, split = "/"))



################
# Figure 3B
################
# Wnt-related GO BP terms of CBC up DEGs (vs all others) in each sample
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers
load(file = sprintf("%s/Rdata/7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata", dir)) #markers.li

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

cluster.updeg = list()
for (i in 1:length(markers.li)){
  deg = markers.li[[i]]
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}

cbc.cont.updeg = cluster.updeg$ctrl$CBC
cbc.1day.updeg = cluster.updeg$'1day'$CBC
cbc.4day.updeg = cluster.updeg$'4day'$CBC

cbc.updeg = list(cont = cbc.cont.updeg, day1 = cbc.1day.updeg, day4 = cbc.4day.updeg)

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

# Wnt terms found in all three samples -> 4 terms shown
des = lapply(wnt.ORA, function(i) {
  i@result$Description
})
inter.des = Reduce(intersect, des)
inter.des = inter.des[c(5,3,9,1)]

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

tiff(filename = sprintf("%s/Figure3B.tiff", out_dir), width = 13, height = 5, units = 'cm', res = 300)
ggplot(merged.data, aes(x = Sample, y = Description, size = p.adjust, color = Count)) +
  geom_point()+
  scale_color_gradient2(low = "lightgrey", mid = "red", high = "red2", midpoint = 15)+
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
# Figure 3C
################
# Mouse CBC Wnt score (SiPSiC, PID_WNT_CANONICAL_PATHWAY mouse orthologs), Day1 vs Day4
load(file = sprintf("%s/Rdata/1-3_merge_GC-PC.Rdata", dir))# s.integrated, markers

library(Seurat)
library(tidyverse)
library(homologene)
library(fgsea)
library(reshape2)
library(SiPSiC)
library(SingleCellExperiment)


gsl = gmtPathways(sprintf("%s/msigdb/c2.all.v2024.1.Hs.symbols.gmt", dir_db)) # MSigDB C2 v2024.1 (human symbols)

wnt = gsl[grep('PID_WNT_CANONICAL_PATHWAY', names(gsl))]
wnt = unique(unlist(wnt))
wnt = homologene(wnt, inTax = 9606, outTax = 10090)  # human -> mouse
wnt = unique(wnt$'10090')

# SiPSiC per sample (RNA layers data.1/2/3 = samples, in the order of unique(s.integrated$sample))
ll = c("data.1", "data.2", "data.3")
for (i in 1:length(unique(s.integrated$sample))){
  sam = unique(s.integrated$sample)[i]
  cell.use = WhichCells(s.integrated, expression = cellType == "CBC" & sample == sam)
  DefaultAssay(s.integrated) = "RNA"
  c.mtx = GetAssayData(s.integrated, assay = "RNA", layer = ll[i])
  sce = SingleCellExperiment(assays = list(counts = as.matrix(c.mtx)))
  scores_and_idx = getPathwayScores(counts(sce), wnt)
  pathway_scores = scores_and_idx$pathwayScores
  gene_indices = scores_and_idx$index
  s.integrated = AddMetaData(object = s.integrated,
                             metadata = setNames(pathway_scores, cell.use),
                             col.name = paste0("WNTscore_", sam))
}

score.cols = grep("^WNTscore_", colnames(s.integrated@meta.data), value = TRUE)

plot.df = FetchData(s.integrated, vars = c("sample", score.cols))
plot.long = melt(plot.df, id.vars = "sample", variable.name = "feature", value.name = "score", na.rm = T)
sample2 = plot.long[plot.long$sample %in% c('1day','4day'),]

grcol = c('lightslateblue','coral1')

tiff(filename = sprintf("%s/Figure3C.tiff", out_dir), width = 4, height = 5, units = 'cm', res = 300)
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
dev.off()



################
# Figure 3F
################
# GSE84061 ISC: Lgr5 log2(FPKM + 1) (FPKM from GEO), Young vs Old
load(file = sprintf("%s/GEO/GSE84061/Rdata/1_GSE84061_FPKM_TPM.Rdata", dir)) #tpm, fpkm

library(ggplot2)

isc.fpkm = fpkm[, c(grep("YoungStem", colnames(fpkm)), grep("AgedStem", colnames(fpkm)))]
colnames(isc.fpkm)  # YoungStem1-3, AgedStem1-3

lgr5.fpkm = data.frame(sample = colnames(isc.fpkm), FPKM = as.numeric(isc.fpkm[match("Lgr5", rownames(isc.fpkm)), ]))
lgr5.fpkm$log2FPKM = log2(lgr5.fpkm$FPKM + 1)
lgr5.fpkm$group = factor(ifelse(grepl("Aged", lgr5.fpkm$sample), "Old", "Young"), levels = c("Young", "Old"))
lgr5.fpkm

# P-value: DESeq2 Wald test (Aged stem vs Young stem, raw counts)
load(file = sprintf("%s/GEO/merge/Rdata/8_merge_ISC.Rdata", dir)) #isc.tpm, isc.deg
lgr5.fpkm.p = isc.deg$pvalue[match("Lgr5", isc.deg$Genes)]
isc.deg[match("Lgr5", isc.deg$Genes), c("Genes", "log2FoldChange", "stat", "pvalue", "padj")]
lgr5.fpkm.p

grcol = c('Young' = 'grey80', 'Old' = 'red3')

tiff(filename = sprintf("%s/Figure3F.tiff", out_dir), width = 4, height = 5.5, units = "cm", res = 300)
ggplot(lgr5.fpkm, aes(x = group, y = log2FPKM, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.8, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.2, shape = 21, fill = "white", stroke = 0.4) +
  scale_fill_manual(values = grcol) +
  scale_y_continuous(limits = c(9, NA), expand = expansion(mult = c(0, 0.35))) +
  labs(x = "", y = "log2(FPKM + 1)", title = "Lgr5 (ISC)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 9),
        plot.title = element_text(size = 10, hjust = 0.5, face = "italic"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
dev.off()



################
# Figure 3G
################
# Human colon CBC (LGR5+ stem cell) Wnt score (SiPSiC, PID_WNT_CANONICAL_PATHWAY), 25-30 vs 65-70 (E-MTAB-9543)
load(file = sprintf("%s/2_E-MTAB-9543/Rdata/2-3_integration.Rdata", dir_geo)) #s.integrated


library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)
library(SiPSiC)


DefaultAssay(s.integrated) = "RNA"
new = JoinLayers(s.integrated, overwrite = T)

gsl = gmtPathways(sprintf("%s/msigdb/c2.all.v2024.1.Hs.symbols.gmt", dir_db)) # MSigDB C2 v2024.1 (human symbols)
wnt = gsl[grep('PID_WNT_CANONICAL_PATHWAY', names(gsl))]
wnt = unique(unlist(wnt))
wnt = intersect(wnt, rownames(new))

# CBC (LGR5+ stem cell), young (25-30) vs old (65-70)
age.use = c("25-30", "65-70")
s.grp = subset(new, cells = colnames(new)[new$cellType %in% "LGR5+ stem cell" & new$age %in% age.use])

# SiPSiC per sample (orig.ident)
wnt.score = setNames(rep(NA_real_, ncol(s.grp)), colnames(s.grp))
for (sam in unique(s.grp$orig.ident)) {
  cell.use = WhichCells(s.grp, expression = orig.ident == sam)
  c.mtx = GetAssayData(s.grp, assay = "RNA", layer = "counts")[, cell.use, drop = FALSE]
  scores_and_idx = getPathwayScores(as.matrix(c.mtx), wnt)
  wnt.score[cell.use] = as.numeric(scores_and_idx$pathwayScores)
}
s.grp$WNT_SiPSiC = wnt.score

# remove outlier cells (1.5 x IQR per age group)
wnt.cell.df = s.grp@meta.data %>%
  dplyr::select(orig.ident, age, WNT_SiPSiC) %>%
  group_by(age) %>%
  filter(between(WNT_SiPSiC, quantile(WNT_SiPSiC, 0.25) - 1.5 * IQR(WNT_SiPSiC),
                 quantile(WNT_SiPSiC, 0.75) + 1.5 * IQR(WNT_SiPSiC))) %>%
  ungroup() %>% as.data.frame()
wnt.cell.df$age = factor(as.character(wnt.cell.df$age), levels = age.use)
table(wnt.cell.df$age)  # cells after outlier removal

grcol = c('darksalmon','saddlebrown')

tiff(filename = sprintf("%s/Figure3G.tiff", out_dir), width = 5, height = 6, units = 'cm', res = 300)
ggplot(wnt.cell.df, aes(x = age, y = WNT_SiPSiC, fill = age)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.7, outlier.shape = NA, linewidth = 0.2) +
  scale_fill_manual(values = grcol) +
  labs(x = "", y = "WNT_SiPSiC", title = "CBC") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 8),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(size = 9, hjust = 0.5))
dev.off()
