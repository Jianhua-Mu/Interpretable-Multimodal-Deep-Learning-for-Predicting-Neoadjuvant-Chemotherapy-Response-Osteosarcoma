# Chemotherapy Response Analysis: Bulk RNA-seq, Immune Infiltration, and Single-cell
# Author: [your name]
# Date: [date]

# Load required packages
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(IOBR)
library(survival)
library(survminer)
library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)


# 1. Differential expression analysis (resistant vs sensitive)

# Assume: expr matrix (genes x samples), group vector indicating resistant/sensitive
group <- factor(c(rep("Sensitive", n_sensitive), rep("Resistant", n_resistant)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# voom transformation for RNA-seq count data
dge <- DGEList(counts = expr)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(Resistant - Sensitive, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# Define significance
deg$change <- ifelse(deg$adj.P.Val < 0.05 & deg$logFC > 1, "Up",
                     ifelse(deg$adj.P.Val < 0.05 & deg$logFC < -1, "Down", "Stable"))
deg <- rownames_to_column(deg, var = "Gene")


# 2. Volcano plot

ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = change), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Stable" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(title = "Chemotherapy Resistant vs Sensitive")

# Label top genes
top_genes <- deg %>% filter(change != "Stable") %>% top_n(10, abs(logFC))
ggplot(deg, aes(logFC, -log10(adj.P.Val), color = change)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Stable" = "grey80")) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 20) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 adjusted P-value")


# 3. Lipid metabolism gene expression (ACOT7, PPT1)

lipid_genes <- c("ACOT7", "PPT1")
expr_lipid <- as.data.frame(t(expr[lipid_genes, ]))
expr_lipid$group <- group
expr_lipid <- pivot_longer(expr_lipid, -group, names_to = "gene", values_to = "expression")

ggplot(expr_lipid, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(aes(group = group), label = "p.format") +
  scale_fill_manual(values = c("Sensitive" = "#2b9672", "Resistant" = "#ce5f17")) +
  theme_bw() +
  labs(y = "Normalized expression", x = "")


# 4. Functional enrichment for upregulated genes in resistant group

gene_up <- deg$Gene[deg$change == "Up"]
entrez_up <- bitr(gene_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(entrez_up$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05)
# Focus on lipid-related terms
lipid_go <- ego@result[grep("lipid|acyl-CoA|palmitoyl|CoA hydrolase", ego@result$Description, ignore.case = TRUE), ]
barplot(ego, showCategory = 10) + ggtitle("GO Biological Process: Resistant vs Sensitive Up")

ekegg <- enrichKEGG(entrez_up$ENTREZID, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
dotplot(ekegg, showCategory = 15) + ggtitle("KEGG pathways: Resistant Up")


# 5. Immune infiltration with xCell

# Input: log2(TPM+1) expression matrix, rows = genes, columns = samples
tme <- deconvo_tme(eset = expr, method = "xcell", arrays = FALSE)
tme$group <- group
# Compare macrophage infiltration
ggplot(tme, aes(x = group, y = Macrophage, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2) +
  stat_compare_means(label = "p.format") +
  scale_fill_manual(values = c("Sensitive" = "#2b9672", "Resistant" = "#ce5f17")) +
  theme_bw() +
  labs(y = "Macrophage infiltration (xCell)", title = "TME Macrophage")

# Compare M2 macrophage signature score (GSVA)
m2_signature <- list(M2 = c("CD163", "IL10", "TGFB1", "MSR1", "CCL18", "MRC1"))
gsva_res <- gsva(as.matrix(expr), m2_signature, method = "ssgsea", kcdf = "Gaussian")
m2_score <- as.numeric(gsva_res)
names(m2_score) <- colnames(expr)
tme$M2_score <- m2_score

ggplot(tme, aes(x = group, y = M2_score, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white") +
  stat_compare_means(label = "p.format") +
  scale_fill_manual(values = c("Sensitive" = "#2b9672", "Resistant" = "#ce5f17")) +
  theme_bw() +
  labs(y = "M2 signature score (ssGSEA)", title = "M2 Macrophage Enrichment")


# 6. Survival analysis (if OS data available)

# Assume surv_data contains columns: OS.time, OS, group (High/Low M2 score)
surv_data <- tme
surv_data$M2_group <- ifelse(surv_data$M2_score > median(surv_data$M2_score), "High", "Low")
fit <- survfit(Surv(OS.time, OS) ~ M2_group, data = surv_data)
ggsurvplot(fit, data = surv_data, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, palette = c("#2b9672", "#ce5f17"),
           legend.labs = c("Low M2", "High M2"),
           xlab = "Time (days)", ylab = "Overall survival probability")


# 7. Single-cell analysis (Seurat workflow) – single-center, no batch correction

rm(list = ls())

pbmc.data <- Read10X(data.dir = "C:/Users/28978/Desktop/xena/single_cell")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)


# Quality control
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- PercentageFeatureSet(pbmc, pattern = '^RP[SL]', col.name = 'percent.RP')

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "orig.ident")
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

pbmc <- NormalizeData(pbmc)

# Remove ribosomal genes
ribo_genes <- grep("^RPL|^RPS", rownames(pbmc), value = TRUE)
pbmc <- pbmc[!rownames(pbmc) %in% ribo_genes, ]


gene_expr_counts <- Matrix::rowSums(pbmc@assays$RNA@counts > 0)
genes_above_1 <- names(gene_expr_counts[gene_expr_counts > 0])


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Dimensional reduction and clustering (using PCA only, no batch correction)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")  # view by sample

pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.03)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE)

# Cell type annotation (using canonical markers)
# Example marker list; adjust according to tissue
sc <- RenameIdents(sc,
                   "0" = "Macrophage",
                   "1" = "T cell",
                   "2" = "Fibroblast",
                   "3" = "Endothelial",
                   "4" = "B cell"
)
sc$cell_type <- Idents(sc)
DimPlot(sc, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend()

# Add M2 score module
sc <- AddModuleScore(sc, features = list(m2_signature$M2), name = "M2_score")
VlnPlot(sc, features = "M2_score1", group.by = "cell_type", pt.size = 0) +
  stat_compare_means(label = "p.signif") +
  labs(y = "M2 signature score")

# If response metadata available (e.g., sc$response <- "Resistant"/"Sensitive")
VlnPlot(sc, features = "M2_score1", group.by = "response", pt.size = 0) +
  stat_compare_means(label = "p.format") +
  scale_fill_manual(values = c("Sensitive" = "#2b9672", "Resistant" = "#ce5f17")) +
  labs(y = "M2 score per cell")

# DotPlot for key markers
markers <- c("CD14", "CD68", "CD163", "IL10", "TGFB1", "ACOT7", "PPT1")
DotPlot(sc, features = markers, group.by = "cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 8. Heatmap of immune infiltration across methods

# Combine multiple deconvolution results (example using IOBR)
tme_list <- list(
  cibersort = deconvo_tme(eset = expr, method = "cibersort", arrays = FALSE),
  epic = deconvo_tme(eset = expr, method = "epic", arrays = FALSE),
  mcp = deconvo_tme(eset = expr, method = "mcpcounter"),
  xcell = deconvo_tme(eset = expr, method = "xcell", arrays = FALSE),
  estimate = deconvo_tme(eset = expr, method = "estimate")
)

# Merge all and annotate by response group
combined <- Reduce(function(x, y) inner_join(x, y, by = "ID"), tme_list)
rownames(combined) <- combined$ID
combined <- combined[, -1]
combined <- t(scale(t(combined)))  # Z-score scaling

# Annotation
ann_col <- data.frame(Response = group, row.names = colnames(expr))
pheatmap(combined, annotation_col = ann_col,
         show_colnames = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
         main = "Immune cell infiltration (multi-algorithm)")

# Visualization with plot1cell
library(plot1cell)
iri.integrated <- pbmc

circ_data <- prepare_circlize_data(iri.integrated, scale = 0.8)
set.seed(1234)

color_celltype <- unique(c(c4a("moonrise3"), c4a('classic10light', 9),
                           "#e5192c", "#3a77b7", "#3cac4c", "#813c93", "#f36c24",
                           "#37b8c3", "#a54922", "#6b7627", "#28996b",
                           "#965b6a", "#e9148f", "#595b5e",
                           "#80d08a", "#d29099", "#f2e010", "#DC143C", "#0000FF",
                           "#20B2AA", "#FFA500", "#9370DB",
                           "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
                           "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3",
                           "#800080", "#A0522D", "#D2B48C",
                           "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD",
                           "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4",
                           "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347", "#6A5ACD",
                           "#9932CC", "#8B008B", "#8B4513", "#DEB887"))

plot_circlize(circ_data, do.label = T, pt.size = 0.5, 
              col.use = color_celltype[1:9],
              bg.color = 'white', kde2d.n = 200, 
              repel = F, label.cex = 0.6)

add_track(circ_data, group = "state", 
          colors = color_celltype[1:2], track_num = 2)
add_track(circ_data, group = "orig.ident",
          colors = color_celltype[1:6], track_num = 3)
add_track(circ_data, group = "RNA_snn_res.0.5",
          colors = color_celltype[1:29], track_num = 4)
add_track(circ_data, group = "celltype",
          colors = color_celltype[1:9], track_num = 5)

# Additional single-cell visualizations (adjust object names as needed)
complex_featureplot(macr, 
                    features = single_cell_function_geneset[[9]], 
                    group = "state", select = c("Re", "Pr"),
                    order = F)

complex_upset_plot(iri.integrated, celltype = "NewPT2", 
                   group = "Group", min_size = 10, logfc = 0.25)

plot_cell_fraction(iri.integrated, 
                   celltypes = c("PTS1", "PTS2", "PTS3", "NewPT1", "NewPT2"), 
                   groupby = "Group", 
                   show_replicate = T, rep_colname = "orig.ident")