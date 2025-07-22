#' Example Seurat objet for a single 384-well plate
#'
#' A real Seurat object (stored in `data/mac.rda`) used in tests and vignettes.
#'
#' @format A Seurat object
#' @source From the PMMSq033 plate
#' @docType data
#' @keywords datasets
"mac"

#' Cell-level metadata for the `mac` dataset
#'
#' The `@meta.data` slot from the `mac` Seurat object, exported as a standalone
#' `data.frame` in `data/metadata.rda`.
#'
#' @format A `data.frame` with columns:
#'   - `orig.ident`: character, plate ID
#'   - `Treatment_1`: character, treatment group
#'   - `Cell_type`:    character, cell type
#'   - `nFeature_RNA`, `nCount_RNA`: numeric, per-cell metrics
#' @docType data
#' @keywords datasets
"metadata"
