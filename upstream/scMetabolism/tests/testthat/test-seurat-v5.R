test_that("counts are read from a single Seurat v5 layer", {
  skip_if_not_installed("Seurat", minimum_version = "5.0.0")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 0, 3, 4, 0, 5),
      nrow = 2,
      dimnames = list(c("gene1", "gene2"), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_equal(scMetabolism:::.get_seurat_counts(obj), counts)
})

test_that("split Seurat v5 counts layers are joined without changing the object", {
  skip_if_not_installed("Seurat", minimum_version = "5.0.0")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 0, 3, 4, 0, 5),
      nrow = 2,
      dimnames = list(c("gene1", "gene2"), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj[["RNA"]] <- split(obj[["RNA"]], f = c("sample1", "sample1", "sample2", "sample2"))
  layers_before <- SeuratObject::Layers(obj[["RNA"]])

  expect_equal(scMetabolism:::.get_seurat_counts(obj), counts)
  expect_equal(SeuratObject::Layers(obj[["RNA"]]), layers_before)
})

test_that("counts are still read from a v3-style assay", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3),
      nrow = 2,
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    ),
    sparse = TRUE
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  suppressWarnings(obj[["RNA"]] <- Seurat::CreateAssayObject(counts = counts))

  expect_s4_class(obj[["RNA"]], "Assay")
  expect_equal(as.matrix(scMetabolism:::.get_seurat_counts(obj)), as.matrix(counts))
})

test_that("metabolism scores form a valid Seurat v5 assay", {
  skip_if_not_installed("Seurat", minimum_version = "5.0.0")

  counts <- matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
  )
  scores <- data.frame(
    cell1 = c(0.1, 0.2),
    cell2 = c(0.3, 0.4),
    row.names = c("pathway1", "pathway2")
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))

  obj <- scMetabolism:::.set_metabolism_scores(obj, scores)

  expect_true(validObject(obj))
  expect_equal(as.matrix(obj@assays$METABOLISM$score), as.matrix(scores))
})

test_that("DimPlot accepts lower-case Seurat v5 reduction keys", {
  skip_if_not_installed("Seurat", minimum_version = "5.0.0")

  counts <- matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
  )
  scores <- data.frame(
    cell1 = 0.1,
    cell2 = 0.3,
    row.names = "pathway1"
  )
  embeddings <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("cell1", "cell2"), c("umap_1", "umap_2"))
  )
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  obj[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = "umap_",
    assay = "RNA"
  )
  obj <- scMetabolism:::.set_metabolism_scores(obj, scores)

  plot <- expect_output(
    DimPlot.metabolism(
      obj,
      pathway = "pathway1",
      dimention.reduction.run = FALSE
    ),
    "Please Cite"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data$metabolism_score, c(0.1, 0.3))
})
