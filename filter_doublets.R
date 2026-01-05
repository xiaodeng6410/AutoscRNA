
####新的代码
print("metadata中不要出现sample和SAMPLE等变量名")
### filter_doublets.R
no_double <- function(SeuratObject) {
    # 获取所有样本ID
    samples <- unique(SeuratObject$orig.ident)
    # 创建空的metadata向量
    doublet_class <- rep(NA, ncol(SeuratObject))
    doublet_score <- rep(NA, ncol(SeuratObject))
    names(doublet_class) <- colnames(SeuratObject)
    names(doublet_score) <- colnames(SeuratObject)
    doublet_cxds_score <- rep(NA, ncol(SeuratObject))
    doublet_weighted <- rep(NA, ncol(SeuratObject))
    names(doublet_cxds_score) <- colnames(SeuratObject)
    names(doublet_weighted) <- colnames(SeuratObject)
    # 存储每个样本过滤后的细胞名称
    cells_to_keep <- c()
    # 对每个样本分别处理
    for (SAMPLE in samples) {
        message("Processing sample: ", SAMPLE )
        # 提取当前样本的细胞
        sample_cells <- colnames(SeuratObject)[SeuratObject$orig.ident == SAMPLE ]
        mat <- GetAssayData(SeuratObject, assay = "RNA", slot = "counts")[, sample_cells]
        # 创建SingleCellExperiment对象
        sce <- SingleCellExperiment(assays = list(counts = mat))
        # 计算该样本的doublet rate
        cellnumber <- ncol(sce)
        doubrate <- (cellnumber * 0.008) / 1000
        message("  Cells: ", cellnumber, ", Doublet rate: ", round(doubrate, 4))
        # 运行scDblFinder
        sce <- scDblFinder(sce, dbr = doubrate)
        # 找出singlet细胞
        singlet_cells <- colnames(sce)[sce$scDblFinder.class == "singlet"]
        cells_to_keep <- c(cells_to_keep, singlet_cells)
        # 保存doublet信息
        doublet_class[sample_cells] <- sce$scDblFinder.class
        doublet_score[sample_cells] <- sce$scDblFinder.score
        doublet_cxds_score[sample_cells] <- sce$scDblFinder.cxds_score
        doublet_weighted[sample_cells] <- sce$scDblFinder.weighted      
    }
    # 添加metadata
    SeuratObject <- AddMetaData(SeuratObject, metadata = doublet_class, col.name = "scDblFinder.class")
    SeuratObject <- AddMetaData(SeuratObject, metadata = doublet_score, col.name = "scDblFinder.score")
    SeuratObject <- AddMetaData(SeuratObject, metadata = doublet_cxds_score, col.name = "scDblFinder.cxds_score")
    SeuratObject <- AddMetaData(SeuratObject, metadata = doublet_weighted, col.name = "scDblFinder.weighted")
    # 只保留所有样本中的singlet细胞
    SeuratObject <- subset(SeuratObject, cells = cells_to_keep)
    message("Total cells kept: ", length(cells_to_keep))
    return(SeuratObject)
}
