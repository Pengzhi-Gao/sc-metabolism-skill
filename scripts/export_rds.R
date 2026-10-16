args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("usage: export_rds.R input.rds expression.csv metadata.csv layer max_elements")
object <- readRDS(args[[1]])
layer <- args[[4]]

if (inherits(object, "Seurat")) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject package is missing")
  assay <- if ("RNA" %in% names(object@assays)) "RNA" else SeuratObject::DefaultAssay(object)
  matches <- SeuratObject::Layers(object[[assay]], search = layer)
  if (length(matches) != 1L || matches[[1]] != layer) stop("RDS assay requires one exact layer name")
  expr <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  meta <- object@meta.data
} else if (inherits(object, "SingleCellExperiment") || inherits(object, "SummarizedExperiment")) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("SummarizedExperiment package is missing")
  if (!layer %in% SummarizedExperiment::assayNames(object)) stop("requested layer/assay not present")
  expr <- SummarizedExperiment::assay(object, layer)
  meta <- as.data.frame(SummarizedExperiment::colData(object))
} else {
  stop("RDS must contain a Seurat or SingleCellExperiment object")
}
if (is.null(rownames(expr)) || is.null(colnames(expr))) stop("RDS expression matrix lacks gene/cell IDs")
if (!all(colnames(expr) %in% rownames(meta))) stop("RDS metadata lacks expression cell IDs")
if (nrow(expr) * ncol(expr) > as.numeric(args[[5]])) stop("RDS matrix exceeds max_dense_elements")
meta <- meta[colnames(expr), , drop = FALSE]
utils::write.csv(as.matrix(expr), args[[2]], quote = TRUE)
utils::write.csv(meta, args[[3]], quote = TRUE)
