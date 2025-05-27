# 设置 Python 环境
setwd("D:/ADS2")
# 设置 Python 环境
library(reticulate)
virtualenv_create("r-reticulate")
use_virtualenv("r-reticulate", required = TRUE)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("anndata", quietly = TRUE)) {
  install.packages("anndata")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# 加载 R 包
  library(Seurat)
  library(anndata)
  library(Matrix)
  library(ggplot2)
  
  # 设置数据集名称
  data <- "pancreas"
  
  # 设置参数
  nHVGs <- 4000
  verbose <- TRUE
  ncpu <- 8
  dimr <- 150
  dimqm <- 20  # 减小维度以降低 NA 风险
  
  # 数据路径
  DATA_DIR <- paste0("./", data, ".h5ad")
  
  # 检查数据集文件是否存在
  if (!file.exists(DATA_DIR)) {
    stop("Dataset file not found: ", DATA_DIR)
  }
  
  # 根据数据集名称设置参考和查询
  if (data == 'pancreas') {
    reference <- c("inDrop1", "inDrop2", "inDrop3", "inDrop4", "fluidigmc1", "smartseq2", "smarter")
    query <- c("celseq", "celseq2")
    batch_key <- "study"
    ct_key <- "cell_type"
  } else {
    stop("Unsupported dataset: ", data)
  }
  
  # 定义 project_only_seurat 函数
  project_only_seurat <- function(adata_file, batch_col, query_names, label_col=NULL, dim=50) {
    print("Reading data.")
    ad <- read_h5ad(adata_file)
    batch_names <- as.character(unique(ad$obs[[batch_col]]))
    is_query <- batch_names %in% query_names
    batch_names <- c(sort(batch_names[!is_query]), sort(batch_names[is_query]))
    batch_order <- order(factor(ad$obs[[batch_col]], levels = batch_names))
    ad <- ad[batch_order, ]
    se <- CreateSeuratObject(counts=t(as(ad$X, "CsparseMatrix")), assay = "RNA")
    VariableFeatures(se) <- rownames(se)
    se[[batch_col]] <- ad$obs[[batch_col]]
    do_transfer <- !is.null(label_col)
    if (do_transfer) {
      se[[label_col]] <- ad$obs[[label_col]]
      n_ref <- sum(!(se[[]][[batch_col]] %in% query_names))
    }
    datas <- SplitObject(se, split.by = batch_col)
    for (i in 1:length(datas)) {
      datas[[i]] <- NormalizeData(datas[[i]], verbose = FALSE, assay = "RNA")
    }
    q_mask <- names(datas) %in% query_names
    refs <- datas[!q_mask]
    queries <- datas[q_mask]
    print("Integrating reference:")
    print(names(refs))
    ref_start <- Sys.time()
    anchors_refs <- FindIntegrationAnchors(object.list = refs, dims = 1:dim)
    ref <- IntegrateData(anchorset = anchors_refs, dims = 1:dim)
    ref_end <- Sys.time()
    ref_time <- difftime(ref_end, ref_start, units = "secs")
    cat("Reference Time (seconds):", ref_time, "\n")
    ref <- ScaleData(ref, assay = "RNA")
    ref <- RunPCA(ref, npcs = dim)
    ref <- FindNeighbors(ref, reduction = "pca", dims = 1:dim, graph.name = "snn")
    ref <- RunSPCA(ref, npcs = dim, graph = "snn")
    ref <- FindNeighbors(ref, reduction = "spca", dims = 1:dim, graph.name = "spca.nn",
                         k.param = 50, cache.index = TRUE, return.neighbor = TRUE, l2.norm = TRUE)
    latent <- Embeddings(ref, reduction = "spca")
    if (do_transfer) {
      pred_labels <- rep(NA, n_ref)
      pred_scores <- rep(NA, n_ref)
    }
    for (i in 1:length(queries)) {
      query <- queries[i]
      print("Mapping query to reference:")
      print(names(query))
      query <- query[[1]]
      query_start <- Sys.time()
      anchors_query <- FindTransferAnchors(reference = ref, query = query, reference.reduction = "spca",
                                           reference.neighbors = "spca.nn", dims = 1:dim, k.anchor = 10)
      if (do_transfer) {
        query <- TransferData(anchorset = anchors_query, reference = ref, query = query, refdata = list(label = label_col))
        pred_labels <- c(pred_labels, query$predicted.label)
        pred_scores <- c(pred_scores, query$predicted.label.score)
        print("Unique pred_labels:", unique(query$predicted.label))
        print("Any NA in pred_labels:", any(is.na(query$predicted.label)))
      }
      query <- IntegrateEmbeddings(anchorset = anchors_query, reference = ref, query = query, reductions = "pcaproject",
                                   dims = 1:dim, new.reduction.name = "qrmapping", reuse.weights.matrix = FALSE)
      query_end <- Sys.time()
      query_time <- difftime(query_end, query_start, units = "secs")
      cat("Query Time (seconds):", query_time, "\n")
      latent <- rbind(latent, Embeddings(query, reduction = "qrmapping"))
    }
    print("Dimensions of latent:", dim(latent))
    print("Any NA in latent:", any(is.na(latent)))
    ad$obsm[["X_seurat"]] <- latent[ad$obs_names,]
    if (do_transfer) {
      ad$obs$pred_label <- pred_labels
      ad$obs$pred_score <- pred_scores
    }
    ad$write_h5ad(adata_file)
  }
  
  # 执行分析
  project_only_seurat(DATA_DIR, batch_key, query, ct_key, dim = dimqm)
  
  # 执行分析
  project_only_seurat(DATA_DIR, batch_key, query, ct_key, dim = dimqm)
  
  # 读取数据
  ad <- read_h5ad(DATA_DIR)
  print(paste("Post-processing obs columns:", paste(names(ad$obs), collapse=", ")))

  
  # 构建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = t(as(ad$X, "CsparseMatrix")), assay = "RNA")
  seurat_obj[["study"]] <- ad$obs[["study"]]
  seurat_obj[["cell_type"]] <- ad$obs[["cell_type"]]
  if ("pred_label" %in% names(ad$obs)) {
    seurat_obj[["pred_label"]] <- ad$obs[["pred_label"]]
  } else {
    warning("pred_label not found in ad$obs")
  }
  
  # 检查 obsm 里是否有 X_seurat，否则自动降维
  if ("X_seurat" %in% names(ad$obsm)) {
    x_seurat <- ad$obsm[["X_seurat"]]
    print(paste("X_seurat dimensions:", paste(dim(x_seurat), collapse = " x ")))
    if (ncol(x_seurat) > 0) {
      seurat_obj[["spca"]] <- CreateDimReducObject(embeddings = x_seurat, key = "SPCA_", assay = "RNA")
      red_name <- "spca"
      dims_use <- 1:min(20, ncol(seurat_obj[["spca"]]))
      message("Using existing X_seurat for降维与UMAP。")
    } else {
      stop("X_seurat has zero columns")
    }
  } else {
    message("X_seurat not found in ad$obsm，自动用PCA降维。")
    # 标准降维流程
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, npcs = 20)
    red_name <- "pca"
    dims_use <- 1:20
  }
  
  # UMAP & 可视化
  seurat_obj <- RunUMAP(seurat_obj, reduction = red_name, dims = dims_use)
  p1 <- DimPlot(seurat_obj, group.by = "cell_type", label = TRUE) + ggtitle("Cell Types")
  p2 <- DimPlot(seurat_obj, group.by = "study", label = TRUE) + ggtitle("Batches")
  print(p1)
  print(p2)
  ggsave("umap_celltype.png", p1, width = 10, height = 8)
  ggsave("umap_batch.png", p2, width = 10, height = 8)
  print(names(ad$obsm))
  latent <- Embeddings(seurat_obj, reduction = "pca") # 或 "pca"
  ad$obsm[["X_seurat"]] <- latent
  
  # 检查meta/obs信息
  print(names(ad$obs))
  
  # 确保X_seurat和meta行顺序一致
  stopifnot(rownames(ad$obsm[["X_seurat"]]) == rownames(ad$obs))
  
  # 保存为新h5ad（建议不要覆盖原文件，便于追溯）
  ad$write_h5ad(paste0(data, "seurat_integrated.h5ad"))
  