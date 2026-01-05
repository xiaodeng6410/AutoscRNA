
###主函数
process_seurat_integrate <- function(seurat_obj, 
                                     norm_methods = c("LogNormalize", "SCT"),
                                     integrate_methods = c("CCA", "RPCA", "Harmony", "scVI"),
                                     output_dir = "./Results/integrated_data/",
                                     verbose = TRUE) {
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # 循环标准化方法
  for (norm_method in norm_methods) {
    if (verbose) message("Processing normalization method: ", norm_method)
    # 标准化
    normalized_obj <- normalize_seurat(seurat_obj, method = norm_method)
    # 循环整合方法
    for (int_method in integrate_methods) {
      # 条件判断：如果是 SCT，就跳过 scVI
      if (norm_method == "SCT" && int_method == "scVI") {
        if (verbose) message("  Skipping scVI integration for SCT-normalized data")
        next
      }
      if (verbose) message("  Processing integration method: ", int_method)
      # 整合
      result  <- integrate_data(normalized_obj, 
                                       norm_method_verbose = norm_method, 
                                       integrate_method_verbose = int_method)
      integrated_obj <- result$seurat_obj_integrated
      int_reduction <- result$reduction_name                                 
      # 聚类
      clustered_obj <- cluster_seurat(integrated_obj, 
                                        reduction = int_reduction,
                                        dims = 1:30,
                                        resolution = 0.5,
                                        verbose = verbose)
      # 构建保存文件名
      save_file <- paste0(output_dir, "total_", norm_method, "_", int_method, "_clustered.rdata")
      if (verbose) message("    Saving result to: ", save_file)
      
      # 保存
      save(clustered_obj, file = save_file)
    }
  }
  if (verbose) message("All combinations processed.")
}

###构建一个降维函数
cluster_seurat <- function(seurat_obj, reduction = "pca", dims = 1:30, resolution = 0.5, verbose = TRUE) {
  if (verbose) message("Clustering using ", reduction, " with dims ", paste(dims, collapse = ","), "...")
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Reduction method ", reduction, " not found in Seurat object.")
  }
  if (reduction == "pca") {
    cluster_name<- paste0("pca", "_cluster")
  }else {
    cluster_name <- paste0(paste(strsplit(reduction, split = "\\.")[[1]], collapse = "_"), "_clusters")
  }
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = dims, verbose = verbose)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = verbose, cluster.name = cluster_name)
  # UMAP 可视化
  umap_name <- paste0("umap.", reduction)
  seurat_obj <- RunUMAP(seurat_obj, reduction.name = umap_name, reduction = reduction, dims = dims, verbose = verbose)
  if (verbose) message("Clustering completed.")
  return(seurat_obj)
}
  
###构建一个标准化函数

# 标准化函数
normalize_seurat <- function(seurat_obj, method = c("LogNormalize", "SCT"),verbose = TRUE) {
  method <- match.arg(method)  # 确保输入合法
  
  if (method == "LogNormalize") {
    if (verbose) message("Running LogNormalize...")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", verbose = verbose)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000, verbose = verbose)
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = verbose)
    seurat_obj <- RunPCA(seurat_obj, verbose = verbose)
    
  } else if (method == "SCT") {
    if (verbose) message("Running SCTransform...")
    options(future.globals.maxSize = 80 * 1024^3)  # 80 GB  
    seurat_obj <- SCTransform(seurat_obj, verbose = verbose)
    seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = verbose)
  }
  
  if (verbose) message("Normalization completed: ", method)
  return(seurat_obj)
}


####构建一个整合函数 

integrate_data <- function(seurat_obj_log, 
                           norm_method_verbose = c("LogNormalize", "SCT"),
                           integrate_method_verbose = c("CCA", "RPCA", "Harmony", "scVI"), 
                           verbose = TRUE) {                        
  # 确保输入合法
  options(future.globals.maxSize = 500 * 1024^3)  # 500 GB  
  norm_method <- match.arg(norm_method_verbose)
  int_method <- match.arg(integrate_method_verbose)
  reduction_name <- NULL
  seurat_obj_integrated <- NULL
  if (verbose) message("Normalizing with ", norm_method, " and integrating with ", int_method)
  if (int_method == "RPCA" && norm_method %in% c("LogNormalize", "SCT")){
      if (norm_method == "SCT" ) {
        reduction_name <- "sct.rpca"
         seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = RPCAIntegration,
                                                  normalization.method = "SCT", 
                                                  new.reduction = reduction_name, 
                                                  verbose = verbose 
                                                  ) ###RPCA SCT 整合
     }else {
        reduction_name <- 'integrated.rpca'
        seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = RPCAIntegration,
                                                  new.reduction = reduction_name, 
                                                  verbose = verbose 
                                                  ) ###RPCA log 整合
    }
    }else if( int_method  == "CCA" && norm_method %in% c("LogNormalize", "SCT")) {
        if (norm_method == "SCT" ) {
            reduction_name <- "sct.cca"
            seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = CCAIntegration,
                                     normalization.method = "SCT",
                                     new.reduction = reduction_name , 
                                     verbose = verbose
                                     ) ###CCA SCT 整合
         }else {
            reduction_name <- 'integrated.cca'
            seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = CCAIntegration,
                                     new.reduction = reduction_name , 
                                     verbose = verbose
                                     )###CCA log 整合
          }
    }else if( int_method == "Harmony" && norm_method %in% c("LogNormalize", "SCT")) {
        if (norm_method == "SCT" ) {
            reduction_name <- "sct.harmony"
            seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = HarmonyIntegration,
                                     normalization.method = "SCT",
                                     new.reduction = reduction_name , 
                                     verbose = verbose
                                     ) ###Harmony SCT 整合
         }else {
            reduction_name <- "integrated.harmony"
            seurat_obj_integrated <- IntegrateLayers(seurat_obj_log, method = HarmonyIntegration,
                                     new.reduction = reduction_name , 
                                     verbose = verbose
                                     )###Harmony log 整合
          }
    }else if(int_method  == "scVI" && norm_method %in% c("LogNormalize")) {
        reduction_name<- "integrated.scvi"
        seurat_obj_integrated <-IntegrateLayers(seurat_obj_log , method = scVIIntegration,
                                                new.reduction = reduction_name,
                                                conda_env = "/home/dengys/anaconda3/envs/scvi", 
                                                verbose = verbose
                                                ) ###scVI  整合
        }else {
                 stop("Unsupported combination: ", norm_method, " + ", int_method)
    }
  # 返回对象
   return(list(
    seurat_obj_integrated = seurat_obj_integrated,
    reduction_name = reduction_name
    ))
}
