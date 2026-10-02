###############
# Figure1A
###############
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

load(file = sprintf("%s/2_E-MTAB-9543/Rdata/2-3_integration.Rdata", dir_geo)) #s.integrated


library(Seurat)
library(tidyverse)
library(Matrix) 
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOfuncR)
library(fgsea)
library(forcats)
library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggtext)
library(factoextra)


DefaultAssay(s.integrated) = "integrated"
uniqueIDs = unique(s.integrated$individual)
age.map = c("417C"='67.5', "411C"='27.5', "A33 (414C)"='52.5') 
s.integrated@meta.data$numeric.age = age.map[as.character(s.integrated$individual)]
expr.mtx = t(as.matrix(
  GetAssayData(s.integrated, layer = "data")
))
ages = as.numeric(s.integrated$numeric.age)

valid.cells = which(!is.na(ages))
expr.mtx = expr.mtx[valid.cells, ]
ages = ages[valid.cells]

age.group = as.character(s.integrated$age[valid.cells])
df = as.data.frame(expr.mtx)
df$age_group = age.group

avg.expr = df %>%
  group_by(age_group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
expr.mat = t(as.matrix(avg.expr[,-1])) 
colnames(expr.mat) = avg.expr$age_group

scaled.expr = t(scale(t(expr.mat)))
scaled.expr = scaled.expr[complete.cases(scaled.expr), ]  

ep = fviz_nbclust(scaled.expr, FUN = hcut, method = "wss", k.max = 10, linecolor = "black") + 
  ggtitle("Elbow method for hierarchical clustering") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        title = element_blank()) 
clu3 = ep$data[2,]

ep + geom_line(size = 0.1, color = "black") +              
  geom_point(size = 1, color = "black") +                
  geom_point(data = clu3, aes(x = clusters, y = y),        
             color = "red", size = 2.5)


scaled.expr = t(scale(t(expr.mat)))
scaled.expr = scaled.expr[complete.cases(scaled.expr), ]  
gene.dist = dist(scaled.expr) 
hc = hclust(gene.dist, method = "ward.D2") 
k = 2  
clusters = cutree(hc, k = k)
clu.col = c("gold", "red")
names(clu.col) = as.character(1:k)

plot.df = scaled.expr %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(cluster = factor(clusters[gene])) %>%
  pivot_longer(cols = -c(gene, cluster),
               names_to = "age_group",
               values_to = "scaled_expr")
plot.df$age_group = factor(plot.df$age_group, levels = c("25-30", "50-55", "65-70"))
plot.df$cluster = factor(plot.df$cluster, labels = paste("Cluster", levels(plot.df$cluster)))

gene.counts = plot.df %>%
  group_by(cluster) %>%
  summarize(n_gene = n_distinct(gene)) %>%
  mutate(cluster_label = paste0(cluster, " (n=", n_gene, ")"))
plot.df$cluster = factor(
  plot.df$cluster,
  levels = gene.counts$cluster,
  labels = gene.counts$cluster_label)

age_levels = levels(factor(plot.df$age_group))
x_start = 1
x_end = length(age_levels)

plot.df$cluster_num = as.character(clusters[plot.df$gene])

tiff(filename = sprintf("%s/figure/Figure1A.tiff", dir), width = 22/3*k, height = 6, units = 'cm', res = 300)
ggplot(plot.df, aes(x = age_group, y = scaled_expr, group = gene)) +
  geom_line(alpha = 0.4, color = "grey", size = 0.3) +
  stat_summary(aes(group = 1, color = cluster_num), fun = mean, geom = "line", size = 1.2) + 
  geom_segment(aes(x = x_start, xend = x_end, y = -2, yend = -2),
               inherit.aes = FALSE, color = "black", size = 0.5) + 
  facet_wrap(~cluster, scales = "fixed", nrow = 1) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(margin = margin(r = -25), color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(size=0.5, color = "black"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)) +
  scale_color_manual(values = clu.col)+
  labs(title = "Gene Clusters by Age-dependent Expression Pattern", x = "Age Group", y = "Z-scaled Expression")
dev.off()


########## Pathways ##########
cluster.genes = split(names(clusters), clusters)
lapply(cluster.genes, function(i) {
  i[grepl('PTGS|PTGES', i)]
})

oral = list()
i=1
for (i in 1:length(cluster.genes)) {
  genes = cluster.genes[[i]]
  msig.wp = read.gmt("E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/MsigDB/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")
  ora = data.frame(enricher(gene =genes, TERM2GENE = msig.wp, pAdjustMethod = "BH",pvalueCutoff = 0.05))
  oral[[names(cluster.genes)[i]]] = ora
}

lapply(oral, function(i) head(i, 10))

top5 = lapply(oral, function(df) {
  df = head(df, 5)
  df$ID = gsub("WP_", "", df$ID)
  df$ID = gsub("_", " ", df$ID)
  df$ID = str_to_sentence(df$ID)
  df$pvalue = -log10(df$pvalue)
  return(df)
})

pg.gene = top5$'2'[4,]$geneID
pg.gene = unlist(strsplit(pg.gene, split = "/"))

bar.plots = lapply(names(top5), function(cl) {
  ggplot(data = top5[[cl]], aes(x = reorder(ID, pvalue), y = pvalue)) +
    coord_flip() +
    #geom_hline(yintercept = -log10(0.01), color = "red", size = 0.5) +
    geom_bar(stat = "identity", width = 0.8, fill = clu.col[cl]) +
    ylab("-log10 P-value") + xlab("") +
    ggtitle(paste('Cluster', cl))+
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
    geom_text(aes(label = ID, y = 0.01), hjust = 0, size = 4, color = "black")+
    scale_y_continuous(expand = c(0, 0)) +  
    scale_x_discrete(expand = c(0, 1))
})
names(bar.plots) = names(top5)
bar.plots


###############
# Figure1B
###############
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

load(file = sprintf("%s/2_E-MTAB-9543/Rdata/2-3_integration.Rdata", dir_geo)) #s.integrated


library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)


DefaultAssay(s.integrated) = "integrated"
age.map = c("417C"='67.5', "411C"='27.5', "A33 (414C)"='52.5') 
s.integrated@meta.data$numeric.age = age.map[as.character(s.integrated$individual)]
expr.mtx = t(as.matrix(
  GetAssayData(s.integrated, layer = "data")
))
ages = as.numeric(s.integrated$numeric.age)

valid.cells = which(!is.na(ages))
expr.mtx = expr.mtx[valid.cells, ]

df = as.data.frame(expr.mtx)
df$age_group = as.character(s.integrated$age[valid.cells])

avg.expr = df %>%
  group_by(age_group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
expr.mat = t(as.matrix(avg.expr[,-1])) 
colnames(expr.mat) = avg.expr$age_group
scaled.expr = t(scale(t(expr.mat)))
scaled.expr = scaled.expr[complete.cases(scaled.expr), ]  
gene.dist = dist(scaled.expr)
hc = hclust(gene.dist, method = "ward.D2")
k = 2  
clusters = cutree(hc, k = k)
clu.col = c("gold", "red")[1:k]
names(clu.col) = as.character(1:k)

cluster.genes = split(names(clusters), clusters)
msig.wp = read.gmt("E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/MsigDB/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")
oral = list()
for (i in 1:length(cluster.genes)) {
  ora = data.frame(enricher(gene = cluster.genes[[i]], TERM2GENE = msig.wp, pAdjustMethod = "BH", pvalueCutoff = 0.05))
  oral[[names(cluster.genes)[i]]] = ora
}

pg.gene = head(oral$'2', 5)[4,]$geneID  
pg.gene = unlist(strsplit(pg.gene, split = "/"))


pg.use = intersect(pg.gene, rownames(scaled.expr))
pg.use

pg.df = scaled.expr[pg.use, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "age_group", values_to = "scaled_expr")
pg.df$age_group = factor(pg.df$age_group, levels = c("25-30", "50-55", "65-70"))
pg.df$cluster_num = as.character(clusters[pg.df$gene])

tiff(filename = sprintf("%s/figure/Figure1B.tiff", dir), width = 10, height = 5, units = 'cm', res = 300)
ggplot(pg.df, aes(x = age_group, y = scaled_expr, group = gene)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3, linetype = "dashed") +  
  geom_line(alpha = 0.6, color = "lightsalmon", size = 0.2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2, color = clu.col["2"]) + 
  #geom_text_repel(data = pg.df %>% filter(age_group == "65-70"),
  #               aes(label = gene), nudge_x = 0.3, direction = "y", hjust = 0,
  #                size = 2.5, fontface = "italic", segment.size = 0.1, max.overlaps = Inf) +
  scale_x_discrete(expand = expansion(add = c(0.2, 0.8))) +  # 오른쪽에 유전자 이름 공간
  coord_cartesian(ylim = c(-2, 2)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(size = 0.5, color = "black"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5)) +
  labs(title = paste0("Prostaglandin synthesis and regulation (n=", length(pg.use), ")"),
       x = "Age Group", y = "Z-scaled Expression")
dev.off()



###############
# Figure1C
###############
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

load(file = sprintf("%s/2_E-MTAB-9543/Rdata/2-3_integration.Rdata", dir_geo)) #s.integrated


library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOfuncR)
library(fgsea)
library(forcats)
library(RColorBrewer)
library(singscore)
library(GSEABase)
library(ggpubr)



DefaultAssay(s.integrated) = "integrated"
uniqueIDs = unique(s.integrated$individual)
age.map = c("417C"='67.5', "411C"='27.5', "A33 (414C)"='52.5') 
s.integrated@meta.data$numeric.age = age.map[as.character(s.integrated$individual)]
expr.mtx = t(as.matrix(
  GetAssayData(s.integrated, layer = "data")
))
ages = as.numeric(s.integrated$numeric.age)

valid.cells = which(!is.na(ages))
expr.mtx = expr.mtx[valid.cells, ]
ages = ages[valid.cells]

age.group = as.character(s.integrated$age[valid.cells])
df = as.data.frame(expr.mtx)
df$age_group = age.group

avg.expr = df %>%
  group_by(age_group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
expr.mat = t(as.matrix(avg.expr[,-1]))  
colnames(expr.mat) = avg.expr$age_group
scaled.expr = t(scale(t(expr.mat)))
scaled.expr = scaled.expr[complete.cases(scaled.expr), ]  
gene.dist = dist(scaled.expr)
hc = hclust(gene.dist, method = "ward.D2")
k = 2 
clusters = cutree(hc, k = k)

cluster.genes = split(names(clusters), clusters)
oral = list()
for (i in 1:length(cluster.genes)) {
  genes = cluster.genes[[i]]
  msig.wp = read.gmt("C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/MsigDB/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")
  ora = data.frame(enricher(gene = genes, TERM2GENE = msig.wp, pAdjustMethod = "BH", pvalueCutoff = 0.05))
  oral[[names(cluster.genes)[i]]] = ora
}

top5 = lapply(oral, function(df) {
  df = head(df, 5)
  df$ID = gsub("WP_", "", df$ID)
  df$ID = gsub("_", " ", df$ID)
  df$ID = str_to_sentence(df$ID)
  df$pvalue = -log10(df$pvalue)
  return(df)
})

pg.gene = top5$'2'[4,]$geneID  
pg.gene = unlist(strsplit(pg.gene, split = "/"))

DefaultAssay(s.integrated) = 'RNA'
new = JoinLayers(s.integrated, overwrite = T)


c.mtx = GetAssayData(new, assay = "RNA", layer = "counts")
pb.counts = sapply(split(colnames(c.mtx), new$orig.ident), function(cells) Matrix::rowSums(c.mtx[, cells, drop = FALSE]))
rownames(pb.counts) = rownames(c.mtx)

pb.cpm = sweep(pb.counts, 2, colSums(pb.counts), "/") * 1e6
sample.age = new@meta.data %>% distinct(orig.ident, age)

pg.gene = intersect(pg.gene, rownames(pb.cpm))
c(n_gene = length(pg.gene))

rankData = rankGenes(pb.cpm)
singscore.df = simpleScore(rankData, upSet = GeneSet(pg.gene))
singscore.df$orig.ident = rownames(singscore.df)
singscore.df = left_join(singscore.df, sample.age, by = "orig.ident")

singscore.df = singscore.df[singscore.df$age %in% c("25-30","50-55", "65-70"), ]
singscore.df$age = factor(singscore.df$age, levels = c("25-30","50-55", "65-70"))
table(singscore.df$age)



library(SiPSiC)
cell.use = colnames(new)[new$age %in% c("25-30","50-55", "65-70")]
cell.cpm = sweep(c.mtx[, cell.use], 2, Matrix::colSums(c.mtx[, cell.use]), "/") * 1e6  
scores_and_idx = getPathwayScores(as.matrix(cell.cpm[pg.gene, , drop = FALSE]), pg.gene)
sipsic.df = data.frame(score = as.numeric(scores_and_idx$pathwayScores), 
                       age = factor(new$age[match(cell.use, colnames(new))], levels = c("25-30","50-55", "65-70")),
                       orig.ident = new$orig.ident[match(cell.use, colnames(new))])
table(sipsic.df$age)



sp.pb = cor.test(as.numeric(singscore.df$age), singscore.df$TotalScore, method = "spearman", alternative = "greater", exact = FALSE)


sp.cell = cor.test(as.numeric(sipsic.df$age), sipsic.df$score, method = "spearman", alternative = "greater", exact = FALSE)
sp.cell

fmt.p = function(p) ifelse(p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(p, 2)))
stat.label = paste0("Pseudobulk: ρ = ", sprintf("%.2f", sp.pb$estimate), ", ", fmt.p(sp.pb$p.value), "\n",
                    "Cell (SiPSiC): ρ = ", sprintf("%.2f", sp.cell$estimate), ", ", fmt.p(sp.cell$p.value))

sipsic.df = sipsic.df %>%
  group_by(age) %>%
  filter(between(score, quantile(score, 0.25) - 1.5 * IQR(score),
                        quantile(score, 0.75) + 1.5 * IQR(score))) %>%
  ungroup() %>%
  as.data.frame()
table(sipsic.df$age)  # 제거 후 cell 수


sip.rng = quantile(sipsic.df$score, c(0.05, 0.95), na.rm = TRUE)
ss.rng = range(singscore.df$TotalScore)
b = unname(diff(ss.rng) / diff(sip.rng))
a = unname(ss.rng[1] - b * sip.rng[1])
sipsic.df$score.scaled = a + b * sipsic.df$score

grcol = c('darksalmon',"salmon",'saddlebrown')

tiff(filename = sprintf("%s/figure/Figure1C.tiff", dir), width = 7, height = 5, units = 'cm', res = 300)
ggplot(singscore.df, aes(x = age, y = TotalScore)) +
  geom_violin(data = sipsic.df, aes(y = score.scaled, fill = age), width = 0.8, alpha = 0.6, scale = "width", linewidth = 0.3) +  # cell (SiPSiC)
  geom_jitter(shape = "★", width = 0.1, size = 4,color = "red") +  # pseudobulk 샘플 (singscore)
  stat_summary(aes(group = 1), fun = median, geom = "line",
               color = "red", linewidth = 0.4, linetype = "dashed") +  # 샘플 그룹별 중앙값 추세선
  scale_fill_manual(values = grcol) +
  annotate("text", x = 0.6, y = Inf, label = stat.label, hjust = 0, vjust = 1.5, size = 2.2) +
  scale_y_continuous(name = "Prostaglandin score (singscore, pseudobulk)",
                     expand = expansion(mult = c(0.05, 0.35)),  # rho/p 글씨(2줄) 들어갈 공간
                     sec.axis = sec_axis(~ (. - a) / b, name = "Prostaglandin score (SiPSiC, cell)")) +
  labs(x = "") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"),
        axis.title.y.left = element_text(color = "red"), axis.text.y.left = element_text(color = "red"),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.1),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()




###############
# Figure1D
###############
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"
dir_geo = "E:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/public_scRNA"

load(file = sprintf("%s/2_E-MTAB-9543/Rdata/2-4_DEG_417C_vs_411C.Rdata", dir_geo)) #deg_417C_vs_411C


library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(ggrepel)
library(stringr)
library(igraph)
library(RCy3)
library(tidyr)
library(tibble)
library(data.table)
library(igraph)
library(RCy3)
library(fgsea)
library(igraph)
library(RCy3)
library(homologene)
library(BiocParallel)
library(ggvenn)
library(EnhancedVolcano)


sigGene = deg_417C_vs_411C[deg_417C_vs_411C$p_val < 0.01 & !is.na(deg_417C_vs_411C$p_val),]
sigGene = sigGene[order(sigGene$avg_log2FC, decreasing = T),]

gsl = read.gmt("C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2/MsigDB/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")
gsl = split(x = gsl$gene, f = gsl$term)
gset = gsl[grep('WP_PROSTAGLANDIN_SYNTHESIS_AND_REGULATION', names(gsl))]
gset = unique(unlist(gset))

updeg = deg_417C_vs_411C[deg_417C_vs_411C$avg_log2FC>1 & deg_417C_vs_411C$p_val<0.01 & !is.na(deg_417C_vs_411C$p_val),] 
inter = intersect(rownames(updeg), gset)

pgeGene = sigGene[rownames(sigGene) %in%  gset,]

plot.df = data.frame(gene = rownames(pgeGene),
                     logfc = pgeGene$avg_log2FC,
                     updeg = ifelse(rownames(pgeGene) %in% inter, "upDEG", "background"))
plot.df$rank = rank(-plot.df$logfc, ties.method = "first") 
topGenes = subset(plot.df, logfc >= 1)

tiff(filename = sprintf("%s/figure/Figure1B.tiff", dir), width = 8, height = 10, units = 'cm', res = 300)
ggplot(plot.df, aes(x = rank, y = logfc)) +
  geom_point(data = subset(plot.df, logfc < 1), shape = 21, col = "darkgrey", fill = "darkgrey", alpha = 0.8, size = 2.5, stroke = 0.1) +
  geom_point(data = subset(plot.df, logfc >= 1), shape = 21, col = "red", fill = "red3", aes(size = ifelse(gene %in% topGenes$gene, 3.5, 2)), stroke = 0.2, alpha = 0.8) +
  scale_size_identity() +
  geom_text_repel(data = topGenes, aes(label = gene), size = 4, color = "red3", segment.size = 0.05, force = 3, nudge_x = 3, fontface = "italic") +
  scale_x_continuous(breaks = plot.df$rank, labels = plot.df$gene) +  
  scale_y_continuous(breaks = seq(-5, 4, 1), limits = c(-5,3.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=15),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_line(color = "grey", linewidth = 0.1),
        panel.border = element_blank(), 
        axis.line = element_line(color="black", size=0.5), 
        axis.line.x.top = element_blank(),              
        axis.line.y.right = element_blank()) +
  labs(x = "", y = "log2FC")
dev.off()



tiff(filename = sprintf("%s/figure/Figure1B.tiff", dir), width = 10, height = 10, units = 'cm', res = 300)
ggplot(plot.df, aes(x = rank, y = logfc)) +
  geom_point(data = subset(plot.df, logfc < 1),
             shape = 21, col = "darkgrey", fill = "darkgrey",
             alpha = 0.8, size = 2.5, stroke = 0.1) +
  geom_point(data = subset(plot.df, logfc >= 1),
             shape = 21, col = "red", fill = "red3",
             aes(size = ifelse(gene %in% topGenes$gene, 3.5, 2)),
             stroke = 0.2, alpha = 0.8) +
  scale_size_identity() +
  geom_text_repel(
    data = topGenes,
    aes(label = gene,
        size  = ifelse(gene == "PTGES", 5, 4),
        color = ifelse(gene == "PTGES", "red3", "black")),
    segment.size = 0.05, force = 3, nudge_x = 4, show.legend = FALSE
  ) +
  scale_color_identity() +
  scale_x_continuous(breaks = plot.df$rank, labels = plot.df$gene) +
  scale_y_continuous(breaks = seq(-5, 4, 1), limits = c(-5, 3.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=15),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_line(color = "grey", linewidth = 0.1),
        panel.border = element_blank(),
        axis.line = element_line(color="black", size=0.5),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank()) +
  labs(x = "", y = "log2FC")
dev.off()



###############
# Figure1E
###############
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/GEO/merge/Rdata/8_merge_DEG.Rdata", dir)) #epi1.deg, epi2.deg, paneth.deg, sm6Stem.deg, deg.li

source("C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rscripts/11-2_GSEAplot.R")


library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(ggrepel)
library(stringr)
library(igraph)
library(RCy3)
library(tidyr)
library(tibble)
library(data.table)
library(igraph)
library(RCy3)
library(fgsea)
library(igraph)
library(RCy3)
library(lsa)


# Epi1
idx = !is.na(epi1.deg$stat)
use = epi1.deg[idx,]
rk = use$stat
names(rk) = use$Genes
rk = sort(rk, decreasing = T)

pg = get_anno_genes('GO:0006693', database = "org.Mm.eg.db")
pg.gset = split(x = pg$gene, f = pg$go_id)

gs = as.data.frame(msigdbr(species = "Mus musculus", category = "C2")) 
pg.gs = unique(gs$gs_name[grep('PROSTAGLANDIN', gs$gs_name)])
pg.gs = pg.gs[c(5,6)]

pg.gset.wp = gs[gs$gs_name %in% pg.gs, c("gs_name", "gene_symbol", "gs_exact_source")] 
pg.gset.wp = split(x=pg.gset.wp$gene_symbol, f=pg.gset.wp$gs_name)

gset = pg.gset.wp[1]
fname = sprintf("%s/figure/FigureS1.B_epi1.tif", dir)
labs = list(mt="Prostaglandin signaling", redgroup.lab="Aged", bluegroup.lab="Young", mlab="GSE213723")
gseaPlot(fname, rk, gset, labs, xmax = 28000) 
#xmax = length(rk)





########################
# Figure 1F
########################
# kegg pathview
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/GEO/merge/Rdata/8_merge_DEG.Rdata", dir)) #epi1.deg, epi2.deg, paneth.deg, sm6Stem.deg, deg.li

library(pathview)

deg.entrez = bitr(epi1.deg$Genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db")
deg.entrez[grep('Ptgs|Ptges', deg.entrez$SYMBOL),]

logfc = epi1.deg$log2FoldChange[match(deg.entrez$SYMBOL, epi1.deg$Genes)]
names(logfc) = epi1.deg$Genes[match(deg.entrez$SYMBOL, epi1.deg$Genes)]
deg.entrez$log2FC = logfc 
entrez.log2fc = deg.entrez$log2FC
names(entrez.log2fc) = deg.entrez$ENTREZID

deg.entrez[grep('19224|19225', deg.entrez$ENTREZID),]

setwd(sprintf("%s/pathview", dir))
pathview(gene.data = entrez.log2fc, pathway.id = "mmu00590", species = "mmu") 



