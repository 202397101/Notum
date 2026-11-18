# Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration


**Abstract**  
Aging is associated with increased intestinal prostaglandin E2 (PGE2) production and impaired resolution of inflammation, creating a chronic inflammatory environment that compromises epithelial homeostasis by driving intestinal stem cells (ISCs) toward secretory differentiation. Here, we analyzed single-cell RNA-sequencing datasets from young and aged human and mouse intestine and identified age-associated alterations in eicosanoid-response programs. We then investigated the epithelial consequences of acute and sustained PGE2 exposure in vivo. Acute PGE2 exposure elicited a regenerative response, whereas sustained exposure induced secretory-lineage skewing without depletion of the ISC compartment and was accompanied by attenuated Wnt signaling. Sustained PGE2 exposure markedly increased expression of Notum, which encodes an extracellular Wnt deacylase, in Paneth cells. Pharmacological inhibition of NOTUM mitigated secretory-lineage skewing and restored epithelial regenerative capacity under sustained PGE2 exposure. Notably, sustained PGE2 treatment alone did not reproduce these phenotypes in intestinal organoids, indicating a requirement for inflammatory niche-derived signals. IFNγ, a cytokine associated with chronic intestinal inflammation, induced Notum expression and recapitulated key epithelial features of sustained PGE2 exposure in organoids, including reduced Wnt responsiveness and secretory-lineage skewing. Together, our findings identify NOTUM as an epithelial niche effector linking chronic inflammatory signaling to reduced Paneth cell–ISC Wnt activity, secretory-lineage skewing, and impaired epithelial regeneration.

---

**Overview of repository**

We combined single-cell and bulk transcriptome analyses of human and mouse intestine to show that prostaglandin programs increase in the aged intestine, and that sustained PGE2 exposure induces secretory-lineage skewing through Paneth cell-derived NOTUM. The analyses in this repository cover:

- Age-associated prostaglandin programs in human intestine (scRNA-seq) and mouse intestinal epithelium (bulk RNA-seq)
- scRNA-seq of mouse intestinal epithelium after PGE2 treatment (control, Day 1, Day 4)
- Wnt signaling and Notum in intestinal stem cells
- scRNA-seq of mouse intestinal epithelium after PGE2 with or without NOTUM inhibitor
- Comparison of aged crypts and crypt-derived organoids, IFNγ-treated organoids, and IFNγ–Wnt association in human small intestine (GTEx)

---
R code used to generate the main figures (Figures 1–5) of the manuscript.

- GitHub: https://github.com/202397101/Notum
- Archived version (Zenodo): https://doi.org/10.5281/zenodo.23102321

## Repository contents

| File | Figure | Description |
|---|---|---|
| `Figure 1.R` | Figure 1 | Prostaglandin programs in the aged intestine: age-dependent gene clusters and WikiPathways ORA (1A), prostaglandin synthesis genes (1B), prostaglandin score, singscore pseudobulk + SiPSiC per cell (1C), log2FC of prostaglandin genes (1D), GSEA (1E), KEGG pathview of mmu00590 (1F) |
| `Figure 2.R` | Figure 2 | scRNA-seq of PGE2-treated epithelium (Ctrl, Day 1, Day 4): UMAP (2A), cell type proportion (2B), Augur cell type prioritization (2C), GO heatmap of CBC and GC/PC up-regulated genes (2D), regeneration score (2E), Monocle2 trajectory (2F), pseudotime of State 3 cells (2G) |
| `Figure 3.R` | Figure 3 | Wnt signaling and Notum: GO ORA of Day 1 / Day 4 CBC (3A), Wnt-related GO terms (3B), mouse CBC Wnt score (3C), *Lgr5* in young vs aged ISCs (3F), human CBC Wnt score, young vs old (3G) |
| `Figure 4.R` | Figure 4 | scRNA-seq with PGE2 plus the NOTUM inhibitor ABC99 (iNOTUM): UMAP (4A), cell type proportion (4B), *Lyz1* / *Muc2* expression (4C), Monocle2 trajectory (4D), pseudotime of State 3 cells (4E), GO ORA of CBC (4J) |
| `Figure 5.R` | Figure 5 | Aged crypts and organoids: MSigDB ORA of old organoid vs old crypt DEGs (5D), SecAct secreted protein activity (5E), organoid cell type proportion, control vs IFN-γ (5F), PGE2 metabolism and Wnt target genes in organoid CBC (5G), IFN-γ vs WNT ssGSEA correlation in GTEx small intestine (5H) |
| `gseaPlot.R` | Figure 1E | GSEA enrichment plot function, sourced by `Figure 1.R` |

Panels not listed above were not generated with this code.

## Public datasets

| Dataset | Description | Used in |
|---|---|---|
| [E-MTAB-9543](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9543) | Human intestine scRNA-seq across age groups | Figure 1A–D, 3G |
| [GSE213723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213723) | Mouse intestinal epithelium bulk RNA-seq, aged vs young | Figure 1E–F |
| [GSE84061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84061) | Mouse ISC and Paneth cell bulk RNA-seq, aged vs young | Figure 3F |
| [GSE169368](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169368) | Aged and young mouse crypts and crypt-derived organoids, bulk RNA-seq | Figure 5D–E |
| [GSE174297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174297) | Mouse intestinal organoid scRNA-seq, control vs IFN-γ | Figure 5F–G |
| [GTEx v8](https://gtexportal.org) | Human small intestine (terminal ileum) bulk RNA-seq | Figure 5H |
| [MSigDB v2024.1](https://www.gsea-msigdb.org) | `c2.all.v2024.1.Hs.symbols.gmt`, `c2.cp.wikipathways.v2024.1.Hs.symbols.gmt` | Figures 1, 3, 5 |

## Input files and folder layout

The scripts start from processed R objects: Seurat objects, DEG tables and Monocle2 objects. Set `dir` at the top of each script to the input data folder. Scripts write TIFF files to `figure/`. Run them from the repository root, because `Figure 1.R` sources `gseaPlot.R`.

```
data/
├── Rdata/
│   ├── 1-3_merge_GC-PC.Rdata                         # Fig 2, 3: Seurat object, PGE2 Ctrl / Day1 / Day4
│   ├── 6-4_CBC_celltype_pathway_heatmap.Rdata        # Fig 2D: GO ORA per cell type
│   ├── 7-2_sample_cellType_vs_others_mergeGC-PC_DEG.Rdata  # Fig 3A-B: cell type vs others DEGs per sample
│   ├── 8-4_Trajectory-analysis_mergeGC-PC.Rdata      # Fig 2F-G: Monocle2 object
│   ├── 11-2_Augur.Rdata                              # Fig 2C: Augur results
│   ├── 13-1_ABC99_vs_PGE2_DEG.Rdata                  # Fig 4J
│   ├── 13-1_CBC_vs_others_DEG.Rdata                  # Fig 4J
│   ├── 13-2_5samples_integration.Rdata               # Fig 4A-C: Seurat object, Cont / PGE2 / ABC99
│   └── 13-2_Trajectory-analysis.Rdata                # Fig 4D-E: Monocle2 object
├── GEO/
│   ├── merge/Rdata/8_merge_DEG.Rdata                 # Fig 1E-F: DESeq2 results, aged vs young
│   ├── merge/Rdata/8_merge_ISC.Rdata                 # Fig 3F: DESeq2 results, aged vs young ISC
│   ├── GSE84061/Rdata/1_GSE84061_FPKM_TPM.Rdata      # Fig 3F: FPKM from GEO
│   └── GSE169368_organoid-crypt/
│       ├── 5_DEG_old_organoid-vs-crypt.Rdata         # Fig 5D
│       └── 6_SecAct.Rdata                            # Fig 5E
├── MsigDB/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt  # Fig 1
├── public_scRNA/
│   ├── 2_E-MTAB-9543/Rdata/2-3_integration.Rdata     # Fig 1A-C, 3G: Seurat object
│   ├── 2_E-MTAB-9543/Rdata/2-4_DEG_417C_vs_411C.Rdata  # Fig 1D
│   └── 4_GSE174297_organoid/Rdata/2_sobj.Rdata       # Fig 5F-G: Seurat object
└── DB/
    ├── msigdb/c2.all.v2024.1.Hs.symbols.gmt          # Fig 3C, 3G, 5H
    └── GTEx/GTEx_v8_SmallIntestine_tpm.Rdata         # Fig 5H
```

Each `# Figure ...` section loads its own input and packages, so you can run sections one at a time.

## Software

The code was tested with R 4.6.1 and these package versions:

| Package | Version | Package | Version |
|---|---|---|---|
| Seurat | 5.5.1 | SeuratObject | 5.4.0 |
| tidyverse | 2.0.0 | ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 | patchwork | 1.3.2 |
| Matrix | 1.7.5 | SingleCellExperiment | 1.34.0 |
| clusterProfiler | 4.20.2 | org.Mm.eg.db | 3.23.0 |
| org.Hs.eg.db | 3.23.1 | GOfuncR | 1.31.0 |
| fgsea | 1.38.0 | BiocParallel | 1.46.0 |
| msigdbr | 26.1.1 | factoextra | 2.2.0 |
| singscore | 1.32.0 | GSEABase | 1.74.0 |
| SiPSiC | 1.12.0 | pathview | 1.52.0 |
| homologene | 1.4.68.19.3.27 | reshape2 | 1.4.5 |
| pheatmap | 1.0.13 | monocle | 2.40.0 |
| RColorBrewer | 1.1.3 | ggpubr | 1.0.0 |
| GSVA | 2.6.15 | ppcor | 1.1 |

## Citation

If you use this code, please cite the article and the archived code (https://doi.org/10.5281/zenodo.23102321).

## Contact

Kyuwon Son (Pusan National University)
