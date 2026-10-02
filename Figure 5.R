################
# Figure 5E,F
################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/aging_PGE2"

load(file = sprintf("%s/GEO/GSE169368_organoid-crypt/2_RCM_TPM.Rdata", dir)) #annot, sinfo, ginfo, tpm, rcm
load(file = sprintf("%s/GEO/GSE169368_organoid-crypt/2_DEG.Rdata", dir)) #degl, degl1

library(ggplot2)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOfuncR)
library(enrichplot)
library(ggvenn)
library(gridExtra)
library(homologene)
library(fgsea)
library(ComplexHeatmap)
library(colorRamp2)
library(pheatmap)



gol=list()
i=2
for (i in 1:length(degl1)){
  deg = degl1[[i]]
  updeg = deg[deg$log2FoldChange > 0 & deg$padj < 0.05 & !is.na(deg$padj),]
  downdeg = deg[deg$log2FoldChange < 0 & deg$padj < 0.05 & !is.na(deg$padj),]
  
  upora = enrichGO(gene = updeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
  downora = enrichGO(gene = downdeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
  
  gol[[names(degl1)[i]]] = list(up=upora, down=downora)
}

lapply(gol, function(i) {
  lapply(i, function(j) head(j@result, 20))
})

go.bpl=list()
i=1
for(i in 1:length(gol)){
  go = gol[[i]]
  up = go$up
  down = go$down
  upf = data.frame(gofilter(up, level = 5))
  upf$p.adjust = -log10(upf$p.adjust)
  dnf = data.frame(gofilter(down, level = 5))
  dnf$p.adjust = -log10(dnf$p.adjust)
  
  up.bp = ggplot(data = head(upf, 7), aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
    coord_flip() +
    geom_bar(stat = "identity", width = 0.8, fill = "orange") +
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
    geom_text(aes(label = Description, y = 0.01), hjust = 0, size = 4, color = "black")+
    scale_y_continuous(expand = c(0, 0)) +  
    scale_x_discrete(expand = c(0, 1))
  
  down.bp = ggplot(data = head(dnf, 7), aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
    coord_flip() +
    geom_bar(stat = "identity", width = 0.8, fill = "deepskyblue") +
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
    geom_text(aes(label = Description, y = 0.01), hjust = 0, size = 4, color = "black")+
    scale_y_continuous(expand = c(0, 0)) +  
    scale_x_discrete(expand = c(0, 1))
  
  go.bpl[[names(gol)[i]]] = list(up=up.bp, down=down.bp)
}

grid.arrange(go.bpl[[1]][[1]], go.bpl[[1]][[2]], ncol = 2)
grid.arrange(go.bpl[[2]][[1]], go.bpl[[2]][[2]], ncol = 2) 
