#####################################################################
# gseaPlot(): GSEA enrichment plot (running ES curve + hit bars + ranked-metric colour bar)
# Used by: Figure 1.R (Figure 1E)
#
# Code for: "Paneth cell-derived Notum to drive lineage skewing and impair intestinal regeneration"
# GitHub : https://github.com/202397101/Notum
# Zenodo : https://doi.org/10.5281/zenodo.23102321
#####################################################################
library(fgsea)
library(BiocParallel)

# fname    : output TIFF file
# rk       : named numeric vector of ranked statistics (names = gene symbols)
# gset     : list with ONE gene set
# labs     : list(mt = title, redgroup.lab = left label, bluegroup.lab = right label, mlab = centre label)
# xmax     : x-axis maximum (number of ranked genes)
# hcol     : colours of the ranked-metric bar (negative, 0, positive)
# pv.x     : right end of the P/NES text (fraction of xmax; smaller = further left)
# pv.cex   : size of the P/NES text
# xmax.lab : TRUE = last x-axis tick shows the number of ranked genes (length(rk))
gseaPlot = function(fname, rk, gset, labs, xmax, hcol = c("blue", "white", "red"), pv.x = 0.99, pv.cex = 1, xmax.lab = FALSE) {
	mt = labs$mt
	redgroup.lab = labs$redgroup.lab
	bluegroup.lab = labs$bluegroup.lab
	class.name = labs$mlab
	a = plotEnrichment(gset[[1]], rk)
	fres = a$data
	colnames(fres) = c("x", "y")
	es.range <- c(min(fres$y), max(fres$y))
	gsea.rnk = data.frame(names(rk), rk)
	colnames(gsea.rnk) = c("id", "metric")
	gsea.rnk = gsea.rnk[order(-gsea.rnk$metric), ]
	metric.range = c(min(gsea.rnk$metric), max(gsea.rnk$metric))
	ss = fgsea(pathways = gset, stats = rk, nperm = 10000, nproc = 0, BPPARAM = SerialParam())
	nes = round(ss$NES, digits = 2)
	# P < 0.0001 -> "P < 0.0001"; P < 0.01 -> 2 significant digits; otherwise 2 decimals
	pv = ifelse(ss$pval < 0.0001, "P < 0.0001", paste("P =", ifelse(ss$pval < 0.01, format(signif(ss$pval, digits = 2), scientific = FALSE), round(ss$pval, digits = 2))))
	tiff(filename = fname, width = 8, height = 7, units = 'cm', res = 300, type = "cairo")

	gsea.layout = layout(matrix(c(1, 2, 3)), heights = c(1.5, 0.2, 0.6))
	# panel 1: running enrichment score
	par(mar = c(0, 3.5, 0.1, 0.8))
	ln1 = round(seq(es.range[1], es.range[2], 0.1), digits = 1)
	plot(fres$x, fres$y, type = "l", col = "green", lwd = 2.5, xaxt = "n", xlab = "", ylab = ""
		, xaxs = "i", yaxt = "n", yaxs = "i", main = "", xpd = T
		, ylim = c(es.range[1] - 0.05, es.range[2] + 0.05), xlim = c(0, xmax)
		, panel.first = { abline(h = ln1, col = "gray95", lty = 2)
		    abline(h = 0, col = "black", lty = 2)
		    abline(h = es.range[ifelse(nes < 0, 1, 2)], col = "red", lty = 2)})
	ln2 = round(seq(es.range[1], es.range[2] + 0.05, 0.2), digits = 1)
	axis(2, at = ln2, labels = rep("", length(ln2)), cex.axis = 1, pos = 0, tck = -0.03)
	axis(2, at = ln2, labels = ln2, las = 2, cex.axis = 1, pos = -length(gsea.rnk$metric) * 0.0001, tick = FALSE)
	axis(2, at = es.range[1] + (es.range[2] - es.range[1]) / 2, labels = "Enrichment score", tick = FALSE, cex.axis = 1.2, pos = -length(gsea.rnk$metric) * 0.068)
	mtext(mt, line = 0.2, cex = 0.8)
	plot.coordinates <- par("usr")
	text(xmax * pv.x, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.15)
		, paste0(pv, "\nNES: ", nes), adj = c(1, 1), cex = pv.cex)
	# panel 2: positions of gene-set members
	par(mar = c(0, 3.5, 0.1, 0.8))
	plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n", ylab = "", xlim = c(1, xmax))
	abline(v = fres$x, lwd = 0.8)
	# panel 3: ranked-metric colour bar
	mx = min(abs(quantile(gsea.rnk$metric, c(0.01, 0.99))))
	metric.range = c(-mx, mx)
	rank.colors = gsea.rnk$metric
	rank.colors[rank.colors >= mx] = mx
	rank.colors[rank.colors <= -mx] = -mx
	rank.colors = (rank.colors - metric.range[1]) / (metric.range[2] - metric.range[1])
	rank.colors = ceiling(rank.colors * 255 + 1)
	rank.colors = colorRampPalette(hcol)(256)[rank.colors]
	rank.colors = rle(rank.colors)

	par(mar = c(3.3, 3.5, 0.5, 0.8))
	bp = barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n"
		, xlim = c(0, xmax), width = 0.2)
	box()
	text(length(gsea.rnk$metric) / 2, 0.1, labels = ifelse(!missing(class.name), class.name, ""), cex = 1)
	text(length(gsea.rnk$metric) * 0.03, 0.1, redgroup.lab, adj = c(0, NA), cex = 1.5)
	text(length(gsea.rnk$metric) * 0.97, 0.1, bluegroup.lab, adj = c(1, NA), cex = 1.5)
	ln3 = seq(0, length(rk), 5000)
	if (xmax.lab) ln3 = c(ln3[ln3 <= length(rk) - 2500], length(rk))  # last tick = number of ranked genes
	axis(1, at = ln3, labels = rep("", length(ln3)), cex.axis = 0.7, pos = -0.05, tck = -0.08, xpd = T)
	ln3.adj = ifelse(ln3 == max(ln3) & xmax.lab, 1, 0.5)  # right-align the last tick label
	for (k in seq_along(ln3)) axis(1, at = ln3[k], labels = ln3[k], cex.axis = 1, pos = 0.05, tick = FALSE, xpd = T, hadj = ln3.adj[k])
	mtext("Gene rank", side = 1, line = 2.1, at = xmax / 2, cex = 0.67)
	dev.off()
}
