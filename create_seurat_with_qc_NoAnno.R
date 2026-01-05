
create_seurat_with_qc <- function(sampleID, 
                                   localdir,
                                   min.cells = 3,
                                   min.features = 200,
                                   tiles = 1,
                                   cores = 20,
                                   verbose = TRUE){
  require(Seurat)
  require(DropletQC)  # 用于 nuclear_fraction_tags
  if(verbose) cat("开始处理样本:", sampleID, "\n")
  # ========== 构建数据路径 ==========
  if(verbose) cat("步骤1: 构建数据路径...\n")
  bamdir <- paste0(localdir, "/cellranger_out/", sampleID, 
                   "/outs/possorted_genome_bam.bam")
  datadir <- paste0(localdir, "/cellbender_out/", sampleID, "/10x_mtx")                 
  barcodedir <- paste0(datadir,"/barcodes.tsv.gz")
  metadir <- paste0(datadir, "/metadata.csv")
  # 检查文件是否存在
  if(!file.exists(bamdir)) stop("BAM文件不存在: ", bamdir)
  if(!file.exists(barcodedir)) stop("Barcode文件不存在: ", barcodedir)
  if(!dir.exists(datadir)) stop("数据目录不存在: ", datadir)
  if(!file.exists(metadir)) stop("Metadata文件不存在: ", metadir)
  # ========== Step 1: 计算核基因比例 ==========
  if(verbose) cat("步骤2: 计算核基因比例 (nuclear fraction)...\n")
  meta_nucli <- nuclear_fraction_tags(
    bam = bamdir,
    barcodes = barcodedir,
    tiles = tiles,
    cores = cores,
    verbose = FALSE
  )
  if(verbose) cat("  - 完成核基因比例计算，共", nrow(meta_nucli), "个细胞\n")
  # ========== Step 2: 读取数据并创建Seurat对象 ==========
  if(verbose) cat("步骤3: 读取数据并创建Seurat对象...\n")
  Seurat.data <- Read10X(data.dir = datadir)
  seurat_obj <- CreateSeuratObject(
    counts = Seurat.data,
    project = sampleID,
    min.cells = min.cells,
    min.features = min.features
  )
  if(verbose) cat("  - Seurat对象创建完成，包含", ncol(seurat_obj), "个细胞\n")
  # 添加核基因比例信息
  barcode <- rownames(seurat_obj@meta.data)
  meta_nucli <- meta_nucli[barcode, , drop = FALSE]
  seurat_obj <- AddMetaData(seurat_obj, meta_nucli)
  if(verbose) cat("  - 已添加核基因比例信息\n")
  # ========== Step 3: 计算EmptyDrops ==========
  if(verbose) cat("步骤4: 识别空液滴 (EmptyDrops)...\n")
  gbm.nf.umi <- data.frame(seurat_obj@meta.data[, c("nuclear_fraction", "nCount_RNA")])
  meta_drops <- identify_empty_drops(nf_umi = gbm.nf.umi)
  if(verbose) cat("  - EmptyDrops识别完成\n")
  # ========== Step 4: 添加CellBender质控信息 ==========
  if(verbose) cat("步骤5: 添加CellBender质控信息...\n")
  meta_cellbender <- read.csv(metadir, header = TRUE, row.names = 1)
  meta_cellbender <- meta_cellbender[barcode, , drop = FALSE]
  meta_cellbender$barcode <- rownames(meta_cellbender)
  # 合并EmptyDrops和CellBender信息
  meta_drops <- merge(meta_drops, meta_cellbender, by = "row.names", all.x = TRUE)
  rownames(meta_drops) <- meta_drops$Row.names
  meta_drops$Row.names <- NULL
  seurat_obj <- AddMetaData(seurat_obj, meta_drops)
  if(verbose) cat("  - 已添加CellBender质控信息\n")
  # ========== Step 5: 计算线粒体和血红蛋白基因比例 ==========
  if(verbose) cat("步骤6: 计算线粒体和血红蛋白基因比例...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(p)]")
  if(verbose) cat("  - 已计算线粒体基因比例 (percent.mt)\n")
  if(verbose) cat("  - 已计算血红蛋白基因比例 (percent.hb)\n")
  # ========== 添加样本信息到metadata ==========
  seurat_obj$sample <- sampleID
  # ========== 重命名一 ==========
  seurat_obj<- RenameCells(seurat_obj, new.names = paste0(sampleid_list[i], "_", colnames(seurat_obj)))
  # ========== 输出质控统计 ==========
  if(verbose) {
    cat("\n=== 质控统计摘要 ===\n")
    cat("样本ID:", sampleID, "\n")
    cat("细胞数:", ncol(seurat_obj), "\n")
    cat("基因数:", nrow(seurat_obj), "\n")
    cat("\nMetadata列名:\n")
    print(colnames(seurat_obj@meta.data))
    cat("\n线粒体基因比例统计:\n")
    print(summary(seurat_obj$percent.mt))
    cat("\n血红蛋白基因比例统计:\n")
    print(summary(seurat_obj$percent.hb))
    cat("\n核基因比例统计:\n")
    print(summary(seurat_obj$nuclear_fraction))
  }
  return(seurat_obj)
}
