suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
option_list <- list(
    make_option(c('-m', '--mode'), type = 'character', help = 'path of celltype'),
    make_option(c('-r', '--rds'), type = 'character', help = 'path of rds'),
    make_option(c('-c', '--celltype'), type = 'character', help = 'path of celltype'),
    make_option(c('-o', '--path_out'), type = 'character', help = 'path of celltype')
)

opt_parser <- OptionParser(option_list = option_list);
opts <- parse_args(opt_parser);

rds <- opts$rds
celltype <- opts$celltype
path_out <- opts$path_out
mode <- opts$mode

# FindMarkers <- function(object, cluster, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
#                                  only.pos = FALSE, return.thresh = 1e-2,
#                                  slot = "data") {
#   matrix = GetAssayData(object, slot = slot)
#   features = rownames(matrix)
#   y <- cluster
#   y <- factor(y)
#   test.res = wilcoxauc(matrix, y)

#   # Calculate logFC
#   if (slot != "scale.data") {
#     if (slot == "data") {
#       X <- expm1(matrix)
#     }
#     group_sums <- sumGroups(X, y, 1)
#     group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
#     cs <- colSums(group_sums)
#     gs <- as.numeric(table(y))
#     lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
#       log(group_means[, g] + 1) - log(((cs - group_sums[g,]) / (length(y) - gs[g])) + 1)
#     }))
#   } else {
#     group_sums = sumGroups(X, y, 1)
#     group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
#     cs <- colSums(group_sums)
#     gs <- as.numeric(table(y))
#     lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
#       group_means[, g] - ((cs - group_sums[g,]) / (length(y) - gs[g]))
#     }))
#   }

#   test.res$avg_logFC <- as.vector(lfc)
#   res <- test.res[, c("pval", "avg_logFC", "pct_in", "pct_out", "padj", "group", "feature")]
#   res[, c("pct_in", "pct_out")] = round(res[, c("pct_in", "pct_out")] / 100, digits = 3)
#   colnames(res) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
#   res <- res %>% dplyr::filter(.data$p_val < return.thresh &
#                          abs(.data$avg_logFC) > logfc.threshold &
#                          (.data$pct.1 > min.pct |
#                          .data$pct.2 > min.pct) &
#                          .data$gene %in% features)
#   if (only.pos) {
#     res <- res %>% dplyr::filter(.data$avg_logFC > 0)
#   }
#   res <- res %>% dplyr::arrange(.data$cluster, .data$p_val, desc(.data$avg_logFC))
#   return(res)
# }

SeuratObj <- readRDS(rds)
celltype <- read.table(celltype, header = TRUE, sep = '\t')[, 2]
if (mode == 'single') {
  SeuratObj$pred <- celltype
} else {
  SeuratObj$pred <- 'pred'
  for (i in 1:length(celltype)) {
    SeuratObj$pred[SeuratObj$seurat_clusters == as.character(i - 1)] = celltype[i]
  }
}

# cluster.genes <- FindMarkers(object = SeuratObj, cluster = celltype)
# cluster.genes <- cluster.genes[cluster.genes$p_val_adj < 1e-05,]
# write.table(cluster.genes, file.path(path_out, paste0(SeuratObj@project.name, "_DiffGenes.tsv")), quote = F, sep = "\t")
p <- DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2, repel = TRUE, group.by = 'pred')
ggsave(file.path(path_out, paste0(SeuratObj@project.name, "_pred.png")), p, width = 7, height = 5)