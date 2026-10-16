#' scMetabolism
#'
#' scMetabolism
#' @param obj
#' @keywords scMetabolism
#' @examples
#' sc.metabolism.Seurat()
#' @export sc.metabolism.Seurat

.get_seurat_counts <- function(obj, assay = "RNA") {
  assay_obj <- obj[[assay]]

  if (is.null(assay_obj)) {
    stop(sprintf("Assay '%s' was not found in the Seurat object.", assay), call. = FALSE)
  }

  # Seurat v5 stores expression matrices in layers.  A merged or integrated
  # object can have several counts layers (for example, counts.sample1 and
  # counts.sample2), which GetAssayData() cannot read as a single matrix.
  if (inherits(assay_obj, "Assay5")) {
    layers <- getExportedValue("SeuratObject", "Layers")(
      assay_obj,
      search = "^counts($|\\.)"
    )

    if (length(layers) == 0L) {
      stop(sprintf("Assay '%s' does not contain a counts layer.", assay), call. = FALSE)
    }

    if (length(layers) > 1L) {
      assay_obj <- getExportedValue("SeuratObject", "JoinLayers")(
        assay_obj,
        layers = "counts",
        new = "counts"
      )
      layer <- "counts"
    } else {
      layer <- layers[[1L]]
    }

    return(getExportedValue("SeuratObject", "LayerData")(
      assay_obj,
      layer = layer
    ))
  }

  # A Seurat v5 installation can still contain a v3-style Assay, but recent
  # SeuratObject releases have made the old slot argument defunct.
  seurat_object_namespace <- asNamespace("SeuratObject")
  if (exists("LayerData", envir = seurat_object_namespace, inherits = FALSE)) {
    return(getExportedValue("SeuratObject", "LayerData")(
      assay_obj,
      layer = "counts"
    ))
  }

  # SeuratObject v4 does not provide LayerData(), so use its slot API.
  Seurat::GetAssayData(obj, assay = assay, slot = "counts")
}

.set_metabolism_scores <- function(obj, signature_exp) {
  seurat_object_namespace <- asNamespace("SeuratObject")

  if (exists("CreateAssay5Object", envir = seurat_object_namespace,
             inherits = FALSE)) {
    create_assay5 <- getExportedValue("SeuratObject", "CreateAssay5Object")
    suppressWarnings({
      metabolism_assay <- create_assay5(data = as.matrix(signature_exp))
      metabolism_assay$score <- as.matrix(signature_exp)
      metabolism_assay$data <- NULL
      obj[["METABOLISM"]] <- metabolism_assay
    })
  } else {
    # Preserve the package's original representation with Seurat v4.
    obj@assays$METABOLISM$score <- signature_exp
  }

  obj
}

.get_metabolism_scores <- function(obj) {
  scores <- obj@assays$METABOLISM$score
  data.frame(as.matrix(scores), check.names = FALSE)
}


sc.metabolism.Seurat <- function(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG") {

  countexp <- .get_seurat_counts(obj, assay = "RNA")

  countexp<-data.frame(as.matrix(countexp))

  #signatures_KEGG_metab <- "./data/KEGG_metabolism_nc.gmt"
  #signatures_REACTOME_metab <- "./data/REACTOME_metabolism.gmt"

  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")


  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}

  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {

    cat("Start imputation...\n")

    #Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588
    #Github: https://github.com/KlugerLab/ALRA

    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")


    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }

  #signature method
  cat("Start quantify the metabolism activity...\n")

  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)

    options(mc.cores = ncores)

    vis <- analyze(vis)

    signature_exp<-data.frame(t(vis@SigScores))
  }

  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }

  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }

  #GSVA
  if (method == "GSVA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")

  .set_metabolism_scores(obj, signature_exp)
}
