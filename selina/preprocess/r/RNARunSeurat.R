RNARunSeurat <- function(inputMat, type = "matrix", project = "query", orig.ident = NULL,
                         min.c = 10, min.g = 200, mito = FALSE, mito.cutoff = 0.2, variable.genes = 2000, assembly = "GRCh38", cluster.res = 0.6, npcs = 50, outdir = ".", mode = "single") {
  if (type == "matrix") { SeuratObj <- CreateSeuratObject(inputMat, project = project, min.cells = min.c, min.features = min.g) }
  if (type == "object") { SeuratObj <- inputMat }

  #=========Mitochondria and Spike-in========  
  if (mito) {
    message("Check the mitochondria and spike-in percentage ...")
    if (assembly == "GRCh38") {
      mito.genes <- grep("^MT-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
      ercc.genes <- grep("^ERCC", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
    }
    else {
      mito.genes <- grep("^mt-", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
      ercc.genes <- grep("^ercc", rownames(GetAssayData(object = SeuratObj)), value = TRUE)
    }
    percent.mito <- Matrix::colSums(GetAssayData(object = SeuratObj)[mito.genes,]) / Matrix::colSums(GetAssayData(object = SeuratObj))
    percent.ercc <- Matrix::colSums(GetAssayData(object = SeuratObj)[ercc.genes,]) / Matrix::colSums(GetAssayData(object = SeuratObj))
    SeuratObj$percent.mito <- percent.mito
    SeuratObj$percent.ercc <- percent.ercc
    p1 <- VlnPlot(SeuratObj, c("percent.mito", "percent.ercc"), ncol = 2)
    ggsave(file.path(outdir, paste0(SeuratObj@project.name, ".spikein.png")), p1, width = 6, height = 4.5)

    SeuratObj <- subset(x = SeuratObj, subset = percent.mito < mito.cutoff)
    SeuratObj <- subset(x = SeuratObj, subset = percent.ercc < 0.05)
    vars.to.regress <- c("nCount_RNA", "percent.mito", "percent.ercc")
  }
  else {
    vars.to.regress <- "nCount_RNA"
  }

  #=========Variable genes========
  message("Normalization and identify variable genes ...")
  SeuratObj <- NormalizeData(object = SeuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObj <- FindVariableFeatures(object = SeuratObj, selection.method = "vst", nfeatures = variable.genes)
  SeuratObj <- ScaleData(object = SeuratObj, vars.to.regress = vars.to.regress)

  #=========PCA===========
  message("PCA analysis ...")
  SeuratObj <- RunPCA(object = SeuratObj, features = VariableFeatures(SeuratObj), npcs = npcs)
  # p2 = ElbowPlot(object = SeuratObj, ndims = SeuratObj@commands$RunPCA.RNA@params$npcs)
  # ggsave(file.path(outdir, paste0(SeuratObj@project.name, "_PCElbowPlot.png")), p2, width = 10, height = 4)

  #=========UMAP===========
  message("UMAP analysis ...")
  pc.contribution <- SeuratObj@reductions$pca@stdev / sum(SeuratObj@reductions$pca@stdev) * 100
  pc.contribution.cum <- cumsum(pc.contribution)
  pc.first <- which(pc.contribution.cum > 75)[1]
  dims.use <- 1:pc.first
  SeuratObj <- RunUMAP(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindNeighbors(object = SeuratObj, reduction = "pca", dims = dims.use)
  SeuratObj <- FindClusters(object = SeuratObj, resolution = cluster.res)
  # p3 <- DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2)
  # ggsave(file.path(outdir, paste0(SeuratObj@project.name, "_cluster.png")), p3, width = 5, height = 4)
  single_exprmat <- as.matrix(SeuratObj@assays$RNA@counts)
  genes <- as.data.frame(rownames(single_exprmat))
  colnames(genes) <- 'Gene'
  single_exprmat <- cbind(genes, single_exprmat)
  cluster_exprmat = as.matrix(AverageExpression(SeuratObj, assays = 'RNA', slot = 'counts')$RNA)
  genes <- as.data.frame(rownames(cluster_exprmat))
  colnames(genes) <- 'Gene'
  cluster_exprmat <- cbind(genes, cluster_exprmat)
  if (mode == 'single') {
    fwrite(single_exprmat, file.path(outdir, paste0(SeuratObj@project.name, "_single_expr.txt")), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  } else if (mode == 'cluster') {
    fwrite(cluster_exprmat, file.path(outdir, paste0(SeuratObj@project.name, "_cluster_expr.txt")), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  } else {
    fwrite(single_exprmat, file.path(outdir, paste0(SeuratObj@project.name, "_single_expr.txt")), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    fwrite(cluster_exprmat, file.path(outdir, paste0(SeuratObj@project.name, "_cluster_expr.txt")), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  }
  return(SeuratObj)
}