
# Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration
---

**Abstract**  
Aging is associated with increased intestinal prostaglandin E2 (PGE2) production and impaired resolution of inflammation, creating a chronic inflammatory environment that compromises epithelial homeostasis by driving intestinal stem cells (ISCs) toward secretory differentiation. Here, we analyzed single-cell RNA-sequencing datasets from young and aged human and mouse intestine and identified age-associated alterations in eicosanoid-response programs. We then investigated the epithelial consequences of acute and sustained PGE2 exposure in vivo. Acute PGE2 exposure elicited a regenerative response, whereas sustained exposure induced secretory-lineage skewing without depletion of the ISC compartment and was accompanied by attenuated Wnt signaling. Sustained PGE2 exposure markedly increased expression of Notum, which encodes an extracellular Wnt deacylase, in Paneth cells. Pharmacological inhibition of NOTUM mitigated secretory-lineage skewing and restored epithelial regenerative capacity under sustained PGE2 exposure. Notably, sustained PGE2 treatment alone did not reproduce these phenotypes in intestinal organoids, indicating a requirement for inflammatory niche-derived signals. IFNγ, a cytokine associated with chronic intestinal inflammation, induced Notum expression and recapitulated key epithelial features of sustained PGE2 exposure in organoids, including reduced Wnt responsiveness and secretory-lineage skewing. Together, our findings identify NOTUM as an epithelial niche effector linking chronic inflammatory signaling to reduced Paneth cell–ISC Wnt activity, secretory-lineage skewing, and impaired epithelial regeneration.

---

**Overview of repository**  
Overview

We combined single-cell and bulk transcriptome analyses of human and mouse intestine to show that prostaglandin programs increase in the aged intestine, and that sustained PGE2 exposure induces secretory-lineage skewing through Paneth cell-derived NOTUM. The analyses in this repository cover:

Age-associated prostaglandin programs in human small intestine (scRNA-seq) and mouse intestinal epithelium (bulk RNA-seq)
scRNA-seq of mouse intestinal epithelium after PGE2 treatment (control, Day 1, Day 4)
scRNA-seq of mouse intestinal epithelium after PGE2 with or without NOTUM inhibitor
Wnt signaling network analysis and selection of Notum among Wnt antagonists
Comparison of aged crypts and crypt-derived organoids, and transcriptomic age estimation in IFNγ-treated organoids

## Code description

**Figure 1. Prostaglandin programs in the aged intestine**
- Fig. 1A: Hierarchical clustering of age-associated genes in human small intestine (E-MTAB-9543)
- Fig. 1B: Expression heatmap of prostaglandin synthesis and regulation genes (WP98) across age groups
- Fig. 1C: Prostaglandin pathway scoring at sample (singscore) and single-cell (SiPSiC) level
- Fig. 1D: Differential expression of prostaglandin genes in EPCAM+ epithelial cells (65–70 vs 25–30 years)
- Fig. 1E: GSEA of prostaglandin metabolic process in aged vs young mouse intestinal epithelial cells (GSE213723)
- Fig. 1F: Mapping of log2 fold changes onto the KEGG arachidonic acid metabolism pathway

**Figure 2. scRNA-seq of PGE2-treated intestinal epithelium**
- Fig. 2A: QC, integration, clustering, and cell-type annotation (UMAP)
- Fig. 2B: Cell-type proportions by condition
- Fig. 2C: Cell-type prioritization by responsiveness to PGE2 (Augur)
- Fig. 2D: GO enrichment of upregulated genes in CBCs and GC/PCs
- Fig. 2E: Regeneration gene set scoring in CBCs
- Fig. 2F–H: Pseudotime trajectory inference (Monocle 2) and cell-type ordering along pseudotime

**Figure 3. Wnt signaling and Notum**
- Fig. 3A–B: GO enrichment of genes upregulated in CBCs at Day 1 and Day 4
- Fig. 3C: Enrichment of Wnt signaling GO terms and child terms in CBCs
- Fig. 3D–E: Wnt and BMP signaling module scores
- Fig. 3H: Log2 fold changes of Wnt upstream regulators in aged vs young bulk RNA-seq datasets (GSE84061, GSE213723, GSE190082)

**Figure 4. scRNA-seq of PGE2 + NOTUM inhibitor-treated intestinal epithelium**
- Fig. 4A: Clustering and cell-type annotation (UMAP)
- Fig. 4B: Cell-type proportions by condition
- Fig. 4C: Expression of Muc2 and Lyz1
- Fig. 4D–E: Pseudotime trajectory inference (Monocle 2)
- Fig. 4J: GO enrichment of upregulated genes in CBCs

**Figure 5. Aged crypts and organoids**
- Fig. 5E–F: GO enrichment of DEGs in aged vs young crypts and crypt-derived organoids (GSE169368)
