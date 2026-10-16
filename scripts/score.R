args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop("usage: score.R expression.csv output.csv method genesets.gmt database cores min_genes r_seed")
input <- args[[1]]
output <- args[[2]]
method <- args[[3]]
gmt <- args[[4]]
database <- args[[5]]
cores <- as.integer(args[[6]])
min_genes <- as.integer(args[[7]])
r_seed <- as.integer(args[[8]])
if (is.na(r_seed) || r_seed < 1L) stop("r_seed must be a positive integer")
set.seed(r_seed)

if (grepl("\\.mtx$", input, ignore.case = TRUE)) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("sparse GSVA input needs Matrix")
  expr <- methods::as(Matrix::readMM(input), "dgCMatrix")
  genes <- utils::read.csv(paste0(input, ".genes.csv"), header = FALSE,
                           colClasses = "character", check.names = FALSE)[[1]]
  cells <- utils::read.csv(paste0(input, ".cells.csv"), header = FALSE,
                           colClasses = "character", check.names = FALSE)[[1]]
  if (length(genes) != nrow(expr) || length(cells) != ncol(expr)) {
    stop("sparse matrix dimensions do not match gene/cell IDs")
  }
  dimnames(expr) <- list(genes, cells)
} else {
  expr <- as.matrix(utils::read.csv(input, row.names = 1, check.names = FALSE))
  storage.mode(expr) <- "numeric"
}
if (anyNA(expr) || any(expr < 0) || nrow(expr) == 0L || ncol(expr) == 0L) {
  stop("expression matrix must contain nonnegative finite values")
}

read_gmt <- function(file, observed, minimum) {
  lines <- readLines(file, warn = FALSE)
  sets <- lapply(strsplit(lines, "\t", fixed = TRUE), function(x) if (length(x) >= 3L) x[-c(1L, 2L)] else character())
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), function(x) x[[1]], character(1))
  lookup <- stats::setNames(observed, toupper(observed))
  sets <- lapply(sets, function(x) unique(unname(lookup[toupper(x)])))
  sets <- lapply(sets, function(x) x[!is.na(x) & nzchar(x)])
  sets[lengths(sets) >= minimum]
}

if (method == "scmetabolism") {
  if (!requireNamespace("scMetabolism", quietly = TRUE)) stop("scMetabolism R package is missing")
  if (!requireNamespace("AUCell", quietly = TRUE)) stop("AUCell R package is missing")
  score <- scMetabolism::sc.metabolism(
    countexp = expr, method = "AUCell", imputation = FALSE,
    ncores = cores, metabolism.type = database
  )
  score <- as.matrix(score)
  # Apply the configured coverage threshold to the installed package's sets.
  native_gmt <- system.file("data", if (database == "KEGG") "KEGG_metabolism_nc.gmt" else "REACTOME_metabolism.gmt",
                            package = "scMetabolism")
  fields <- strsplit(readLines(native_gmt, warn = FALSE), "\t", fixed = TRUE)
  eligible <- vapply(fields, function(x) length(intersect(x[-c(1L, 2L)], rownames(expr))) >= min_genes, logical(1))
  set_names <- vapply(fields, function(x) x[[1]], character(1))
  score <- score[rownames(score) %in% set_names[eligible], , drop = FALSE]
  if (!nrow(score)) stop("no scMetabolism gene set passed min_genes")
} else if (method == "gsva") {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA R package is missing")
  sets <- read_gmt(gmt, rownames(expr), min_genes)
  if (!length(sets)) stop("no metabolic gene set passed min_genes")
  param <- GSVA::gsvaParam(expr, sets, kcdf = "auto", minSize = min_genes,
                           sparse = inherits(expr, "sparseMatrix"))
  score <- GSVA::gsva(param, verbose = FALSE)
} else if (method == "vision") {
  if (!requireNamespace("VISION", quietly = TRUE)) stop("VISION R package is missing")
  if (ncol(expr) < 30L) stop("VISION requires at least 30 cells for stable neighborhood analysis")
  sets <- read_gmt(gmt, rownames(expr), min_genes)
  if (!length(sets)) stop("no metabolic gene set passed min_genes")
  mapped_gmt <- tempfile(fileext = ".gmt")
  writeLines(vapply(names(sets), function(name) paste(c(name, "metabolism", sets[[name]]), collapse = "\t"), character(1)), mapped_gmt)
  depth <- colSums(expr)
  if (any(depth <= 0)) stop("VISION input contains empty cells")
  scaled <- t(t(expr) / depth) * stats::median(depth)
  options(mc.cores = cores)
  vision <- VISION::Vision(scaled, signatures = mapped_gmt)
  vision <- VISION::analyze(vision)
  score <- t(as.matrix(VISION::getSignatureScores(vision)))
  unlink(mapped_gmt)
} else {
  stop(paste("unsupported method", method))
}

if (!is.matrix(score)) score <- as.matrix(score)
if (identical(colnames(score), make.names(colnames(expr), unique = TRUE))) {
  colnames(score) <- colnames(expr)
}
if (!identical(colnames(score), colnames(expr))) {
  if (identical(rownames(score), colnames(expr))) {
    score <- t(score)
  } else {
    stop("backend score cell IDs do not match input")
  }
}
if (anyNA(score) || any(!is.finite(score))) stop("backend returned nonfinite scores")
utils::write.csv(score, output, row.names = TRUE, quote = TRUE)
cat(method, "complete:", nrow(score), "pathways x", ncol(score), "cells\n")
