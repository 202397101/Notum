#####################################################################
# Figure 5 - Aged crypts and organoids, IFN-γ and Wnt
#
# Code for: "Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration"
# GitHub : https://github.com/202397101/Notum
# Zenodo : https://doi.org/10.5281/zenodo.23102321
#
# Data   : GSE169368 (aged/young mouse crypts and organoids, bulk RNA-seq), GSE174297 (IFN-γ-treated organoid scRNA-seq),
#          GTEx v8 (human small intestine)
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
# Figure 5D
################
# MSigDB ORA of DEGs in old organoids vs old crypts (GSE169368, DESeq2)
load(file = sprintf("%s/GEO/GSE169368_organoid-crypt/5_DEG_old_organoid-vs-crypt.Rdata", dir)) #deg.old, updeg, downdeg

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)

# universe: all genes tested by DESeq2 (padj not NA)
universe = deg.old$Genes[!is.na(deg.old$padj)]
deg.set = list(Up = updeg$Genes, Down = downdeg$Genes)

c2 = msigdbr(species = "Mus musculus", collection = "C2")
msig.t2g = list(
  PID          = c2 %>% filter(gs_subcollection == "CP:PID") %>% dplyr::select(gs_name, gene_symbol) %>% distinct(),
  REACTOME     = c2 %>% filter(gs_subcollection == "CP:REACTOME") %>% dplyr::select(gs_name, gene_symbol) %>% distinct(),
  WIKIPATHWAYS = c2 %>% filter(gs_subcollection == "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
)
sapply(msig.t2g, function(t) length(unique(t$gs_name)))

msig.ora = list()
for (db in names(msig.t2g)) {
  for (dr in names(deg.set)) {
    ora = enricher(gene = deg.set[[dr]], universe = universe, TERM2GENE = msig.t2g[[db]],
                   pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
    ora@result$Description = gsub("_", " ", sub("^(PID|REACTOME|WP)_", "", ora@result$ID))
    msig.ora[[paste(db, dr, sep = ".")]] = ora
  }
}
sapply(msig.ora, function(o) sum(o@result$p.adjust < 0.05))  # number of significant terms

head(msig.ora$PID.Up@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)
head(msig.ora$PID.Down@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)
head(msig.ora$REACTOME.Up@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)
head(msig.ora$REACTOME.Down@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)
head(msig.ora$WIKIPATHWAYS.Up@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)
head(msig.ora$WIKIPATHWAYS.Down@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15)


########## REACTOME (down DEGs) barplot: top 5, -log10 FDR ##########
top.n = 5
rea.down = msig.ora$REACTOME.Down@result
rea.down$p.adjust = -log10(rea.down$p.adjust)

rea.down.bp = ggplot(data = head(rea.down, top.n), aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  geom_bar(stat = "identity", width = 0.8, fill = "dodgerblue") +
  ylab("-log10(FDR)") + xlab("") +
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
  geom_text(aes(label = Description, y = 0.01), hjust = 0, size = 4, color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 1))

rea.down.bp



################
# Figure 5E
################
# SecAct (aged vs young, GSE169368): crypt vs organoid activity z-score (first quadrant)
#  - x = organoid SecAct z-score, y = crypt SecAct z-score
#  - proteins with z > 0 in both; FDR < 0.05 (crypt or organoid) coloured by side of the diagonal (y = x)
#  - others = grey (not significant)
load(file = sprintf("%s/GEO/GSE169368_organoid-crypt/6_SecAct.Rdata", dir)) #diff.l, act.l, res.l (SecAct results for crypts and organoids)

library(ggplot2)

act.df = merge(res.l$Crypts[, c("protein", "zscore", "padj")],
               res.l$Organoids[, c("protein", "zscore", "padj")],
               by = "protein", suffixes = c(".crypt", ".org"))

q1.df = act.df[act.df$zscore.org > 0 & act.df$zscore.crypt > 0, ]
q1.df$sig = (q1.df$padj.crypt < 0.05) | (q1.df$padj.org < 0.05)

q1.df$side = ifelse(!q1.df$sig, "Not significant",
                    ifelse(q1.df$zscore.crypt > q1.df$zscore.org, "Crypt > Organoid", "Organoid > Crypt"))
q1.df$side = factor(q1.df$side, levels = c("Not significant", "Crypt > Organoid", "Organoid > Crypt"))
q1.df = q1.df[order(q1.df$side), ]  # coloured points on top
table(q1.df$side)
q1.df[q1.df$protein == "IFNG", ]

side.col = c("Not significant" = "grey", "Crypt > Organoid" = "lightcoral", "Organoid > Crypt" = "skyblue")

tiff(filename = sprintf("%s/Figure5E.tiff", out_dir), width = 7, height = 8, units = "cm", res = 300)
ggplot(q1.df, aes(x = zscore.org, y = zscore.crypt)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(aes(colour = side), size = 2, alpha = 0.6) +
  scale_colour_manual(values = side.col, breaks = c("Crypt > Organoid", "Organoid > Crypt"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Organoid SecAct z-score", y = "Crypt SecAct z-score",
       title = "Secreted protein activity\n(Aged vs. Young)") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        legend.position = "none")
dev.off()




################
# Figure 5F
################
# Cell type proportion of organoids, control vs IFN-γ (GSE174297)
load(file = file.path(dir_geo,"4_GSE174297_organoid","Rdata","2_sobj.Rdata")) #sobj

library(Seurat)
library(tidyverse)
library(ggrepel)

sobj@images = list()


# same colours as Figure 2 (EC progenitor = APC colour)
mycol = c('CBC' = 'greenyellow',
          'TA' = 'darkgreen',
          'EC progenitor' = 'cyan',
          'EC' = 'blue',
          'EE' = 'orange',
          'GC/PC' = 'orangered1',
          'TC' = 'red4')

sobj$cellType = factor(sobj$cellType, levels = names(mycol))
sobj$treatment = factor(sobj$treatment, levels = c('sc_control','sc_ifnr'))


cell.pop = sobj@meta.data %>%
  dplyr::count(treatment, cellType, name = "counts") %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(percent = counts / sum(counts) * 100)
cell.pop$cellType = factor(cell.pop$cellType, levels = names(mycol))

# percentage labels (1 decimal), same as Figure 2B and 4B
# label y position = middle of each stacked segment (same order as position_stack: last level at the bottom)
cell.pop = cell.pop %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = sprintf("%.1f%%", percent)) %>%
  dplyr::group_by(treatment) %>%
  dplyr::arrange(dplyr::desc(cellType), .by_group = TRUE) %>%
  dplyr::mutate(ymax = cumsum(percent), ymid = ymax - percent / 2) %>%
  dplyr::ungroup()

small.cut = 7  # segments < 7% are labelled outside the bar with a leader line
bar.w = 0.9    # bar width
bar.dist = 1.8 # distance between bar centres
dark.ct = c("TA", "EC", "TC")  # white text on dark segments

# x position: control = 1, ifnr = 1 + bar.dist
cell.pop$x.pos = 1 + (as.numeric(cell.pop$treatment) - 1) * bar.dist

tiff(filename = sprintf("%s/Figure5F.tiff", out_dir), width = 12, height = 8, units = 'cm', res = 300)
ggplot(cell.pop, aes(x = x.pos, y = percent, fill = cellType)) +
  geom_bar(stat = "identity", color = "black", size = 0.1, width = bar.w) +
  geom_text(data = subset(cell.pop, percent >= small.cut),
            aes(x = x.pos, y = ymid, label = label, color = ifelse(cellType %in% dark.ct, "white", "black")),
            inherit.aes = FALSE, size = 3.4, fontface = "bold") +
  geom_text_repel(data = subset(cell.pop, percent < small.cut),
                  aes(x = x.pos + bar.w / 2, y = ymid, label = paste0(cellType, " ", label)), inherit.aes = FALSE,
                  size = 2.7, nudge_x = 0.06, direction = "y", hjust = 0, color = "black",
                  segment.size = 0.3, min.segment.length = 0, box.padding = 0.15, point.padding = 0) +
  scale_color_identity() +
  scale_fill_manual(values = mycol) +
  labs(y = "% of cell population") +
  geom_hline(yintercept = 0, color = "black", size = 0.5, linetype = "solid") +
  geom_vline(xintercept = 1 - bar.w / 2 - 0.3, color = "black", size = 0.5, linetype = "solid") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = 1 + (seq_along(levels(cell.pop$treatment)) - 1) * bar.dist, labels = levels(cell.pop$treatment),
                     limits = c(1 - bar.w / 2 - 0.3, 1 + bar.dist + bar.w / 2 + 0.9), expand = c(0, 0)) +
  coord_cartesian(clip = "off")
dev.off()




################
# Figure 5G
################
# Organoid CBC: PGE2 metabolism and Wnt target genes, control vs IFN-γ (GSE174297)
load(file = file.path(dir_geo, "4_GSE174297_organoid", "Rdata", "2_sobj.Rdata")) #sobj

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

sobj@images = list()
DefaultAssay(sobj) = "RNA"

meta = sobj@meta.data
cbc = meta$cellType %in% "CBC"
cell.ctrl = rownames(meta)[cbc & meta$treatment == "sc_control"]
cell.ifn  = rownames(meta)[cbc & meta$treatment == "sc_ifnr"]
cell.neg  = rownames(meta)[cbc & meta$treatment == "sc_ifnr" & meta$Cd74_group == "negative_Cd74"]
cell.pos  = rownames(meta)[cbc & meta$treatment == "sc_ifnr" & meta$Cd74_group == "positive_Cd74"]
c(ctrl = length(cell.ctrl), ifn = length(cell.ifn), cd74.neg = length(cell.neg), cd74.pos = length(cell.pos))  # 764, 318, 217, 37

cmp.levels = c("CBC", "IFN-γ CBC")

# gene groups; genes expressed in < 1% of CBC (Ptgs2, Ptger1-3, Ptges, Lef1, Nkd1, Notum, Tnfrsf19) are excluded
gene.sets = list("PGE2 metabolism" = c("Ptgs1", "Ptges2", "Ptges3", "Abcc4", "Slco2a1", "Hpgd"),  # synthesis -> transport -> degradation
                 "Wnt targets"     = c("Ascl2", "Sp5", "Sox9", "Myc", "Lgr5", "Olfm4", "Smoc2"))
genes = unlist(gene.sets, use.names = FALSE)

# dot plot table (gene x group: mean expression, % expressing, z-score, Wilcoxon p)
expr = FetchData(sobj, vars = genes)
grp.cells = list("Control" = cell.ctrl, "IFN-γ" = cell.ifn, "Cd74-" = cell.neg, "Cd74+" = cell.pos)
dot.df = do.call(rbind, lapply(names(grp.cells), function(g) {
  e = expr[grp.cells[[g]], , drop = FALSE]
  data.frame(gene = genes, group = g, mean = colMeans(e), pct = colMeans(e > 0) * 100)
}))
dot.df$gene.set = factor(rep(rep(names(gene.sets), lengths(gene.sets)), length(grp.cells)), levels = names(gene.sets))
dot.df$comparison = factor(ifelse(dot.df$group %in% c("Control", "IFN-γ"), cmp.levels[1], cmp.levels[2]), levels = cmp.levels)
dot.df = dot.df %>%
  dplyr::group_by(comparison, gene) %>%  # z-score between the two groups of each comparison
  dplyr::mutate(z = (mean - mean(mean)) / sd(mean)) %>%
  dplyr::ungroup() %>%
  as.data.frame()
dot.df$z[is.na(dot.df$z)] = 0  # equal means -> z = 0 (white)

p.ifn  = sapply(genes, function(x) suppressWarnings(wilcox.test(expr[cell.ifn, x], expr[cell.ctrl, x])$p.value))
p.cd74 = sapply(genes, function(x) suppressWarnings(wilcox.test(expr[cell.pos, x], expr[cell.neg, x])$p.value))
dot.df$p = ifelse(dot.df$group == "IFN-γ", p.ifn[dot.df$gene], ifelse(dot.df$group == "Cd74+", p.cd74[dot.df$gene], NA))
dot.df$star = ifelse(is.na(dot.df$p), "", ifelse(dot.df$p < 0.01, "**", ifelse(dot.df$p < 0.05, "*", "")))
dot.df$gene = factor(dot.df$gene, levels = genes)
dot.df


# CBC: Control vs IFN-γ
dot.ifn = dot.df[dot.df$comparison == cmp.levels[1], ]
dot.ifn$group = factor(dot.ifn$group, levels = c("IFN-γ", "Control"))  # top to bottom: Control, IFN-γ

p.dot.ifn = ggplot(dot.ifn, aes(x = gene, y = group)) +
  geom_point(aes(size = pmin(pct, 60), fill = z), shape = 21, color = "dimgrey", stroke = 0.4) +  # dots > 60% drawn at the maximum size
  geom_text(aes(label = star), nudge_y = -0.5, size = 3.5) +  # stars below the IFN-γ dots
  facet_grid(. ~ gene.set, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "magenta", midpoint = 0, name = "Mean expr\n(z-score)",
                       limits = c(-0.8, 0.8), breaks = c(-0.8, -0.4, 0, 0.4, 0.8), oob = scales::squish) +
  scale_size_continuous(range = c(0.5, 5), limits = c(0, 60), breaks = c(0, 20, 40, 60), name = "% expressing") +
  scale_y_discrete(expand = expansion(add = c(0.8, 0.35))) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "", title = "CBC: Control vs IFN-γ") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 10),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))

# gene groups stacked: PGE2 metabolism (top), Wnt targets (bottom); empty space on the right keeps dot spacing equal
n.gene = lengths(gene.sets)
p.ifn.top = p.dot.ifn %+% dot.ifn[dot.ifn$gene.set == names(gene.sets)[1], ]
p.ifn.bot = p.dot.ifn %+% dot.ifn[dot.ifn$gene.set == names(gene.sets)[2], ] + labs(title = NULL)

tiff(filename = sprintf("%s/Figure5G.tiff", out_dir),
     width = 3 + 0.8 * max(n.gene), height = 9, units = "cm", res = 300)
((p.ifn.top + plot_spacer() + plot_layout(widths = c(n.gene[1], max(n.gene) - n.gene[1]))) / p.ifn.bot) &
  theme(legend.position = "none")
dev.off()

# legend (colour -0.8 to 0.8, size 0-60%)
tiff(filename = sprintf("%s/Figure5G_legend.tiff", out_dir), width = 4.2, height = 3, units = "cm", res = 300)
ggpubr::as_ggplot(ggpubr::get_legend(p.dot.ifn +
                                       guides(fill = guide_colorbar(barheight = unit(1.6, "cm"), barwidth = unit(0.25, "cm"))) +
                                       theme(legend.box = "horizontal",
                                             legend.title = element_text(size = 7), legend.text = element_text(size = 6))))
dev.off()




################
# Figure 5H
################
# GTEx v8 small intestine: Reactome IFN-γ signaling vs PID WNT canonical pathway (ssGSEA)
load(file = sprintf("%s/GTEx/GTEx_v8_SmallIntestine_tpm.Rdata", dir_db)) #si.tx.tpm (transcript), si.gene.tpm (gene symbol), samp (sample annotation)

library(dplyr)
library(ggplot2)
library(GSVA)
library(fgsea)
library(ppcor)

log.tpm = log2(si.gene.tpm + 1)

#  - IFN-γ set = REACTOME_INTERFERON_GAMMA_SIGNALING (top 5 of the Figure 5D REACTOME ORA), human symbols (MSigDB C2)
#  - ssGSEA normalizes over all gene sets -> the three sets (IFN-γ, WNT_PID, WNT_targets) are scored together
gsl = gmtPathways(sprintf("%s/msigdb/c2.all.v2024.1.Hs.symbols.gmt", dir_db))
gs.rea = list(
  IFNG_REACTOME = gsl[["REACTOME_INTERFERON_GAMMA_SIGNALING"]],
  WNT_PID       = gsl[["PID_WNT_CANONICAL_PATHWAY"]],
  WNT_targets   = c("AXIN2", "LGR5", "ASCL2", "SP5", "OLFM4", "SMOC2", "SOX9", "MYC",
                    "NKD1", "CD44", "RNF43", "ZNRF3", "CCND1", "EPHB2")
)
gs.rea = lapply(gs.rea, intersect, rownames(log.tpm))
lengths(gs.rea)

sc.rea = gsva(ssgseaParam(log.tpm, gs.rea, minSize = 5), verbose = FALSE)
cor.rea = data.frame(samp[colnames(sc.rea), c("SAMPID", "AGE", "age.num")], t(sc.rea), check.names = FALSE)

# Spearman correlation and partial correlation adjusted for age
sp.rea = cor.test(cor.rea$IFNG_REACTOME, cor.rea$WNT_PID, method = "spearman", exact = FALSE)
pc.rea = pcor.test(cor.rea$IFNG_REACTOME, cor.rea$WNT_PID, cor.rea$age.num, method = "spearman")
c(rho = unname(sp.rea$estimate), p = sp.rea$p.value, partial_rho_age = pc.rea$estimate, partial_p_age = pc.rea$p.value)

tiff(filename = sprintf("%s/Figure5H.tiff", out_dir), width = 7, height = 7, units = "cm", res = 300)
ggplot(cor.rea, aes(x = IFNG_REACTOME, y = WNT_PID)) +
  geom_point(shape = 16, colour = "darkorange", size = 1, alpha = 0.85) +
  geom_smooth(method = "lm", se = T, colour = "red3", linewidth = 0.6, linetype = "solid") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "Reactome IFN-γ signaling (ssGSEA)", y = "WNT_PID (ssGSEA)", title = "GTEx small intestine") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 8),
        panel.grid = element_blank(),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))
dev.off()
