#' Compute median metrics by control groups
#'
#'@import dplyr
#'@import tibble


get_median_metrics <- function (metadata = NULL,
                                metrics = "nFeature_RNA",
                                group_by = "Treatment_1",
                                controls = c("^DMSO|^Media|^Staurosporine|^Untreated"),
                                Model_type = NULL){
  # Check if metadata is a data frame
  if (!is.data.frame(metadata)) {
    stop("Metadata must be a data frame.")
  }
  # Check if metrics is a character vector
  if (!is.character(metrics) | !(metrics %in% c("nFeature_RNA", "nCount_RNA"))) {
    stop("Metrics must be a character vector. It is either nFeature_RNA or nCount_RNA.")
  }
  # Check if group_by is NULL
  if (is.null(group_by)) { group_by <- "Treatment_1" }
  # Check if group_by is in the columns in metadata
  if (is.null(metadata[[group_by]])) {
    stop(paste("Group_by", group_by, "is not a column in metadata."))
  }
  # Check if controls is a character vector
  if (!is.character(controls)) {
    stop("Controls must be a character vector.")
  }
  # Check model type if it is one of three types
  if (is.null(Model_type) | !Model_type %in% c("2D_non_adherent", "2D_adherent", "3D")){
    stop("Model_type must be one of '2D_non_adherent', '2D_adherent', or '3D'.")
  }


  # For entire plate
  # Overall
  overall_median <- data.frame(sample_type = "Overall",
                               median_val = median(metadata[[metrics]], na.rm = TRUE))
  colnames(overall_median) <- c("sample_type", paste0("median_", metrics))
  # Filter metadata for controls
  control_metadata <- metadata[grepl(paste(controls, collapse = "|"), metadata[[group_by]]), ]
  control_median <- control_metadata %>%
    group_by(.data[[group_by]]) %>%
    summarise(median_gene_numbers = median(.data[[metrics]]))
  colnames(control_median) <- c("sample_type",  paste0("median_", metrics))

  median_metrics <- rbind(overall_median, control_median)

  # Standardise names of controls
  median_metrics$sample_type <- str_replace(median_metrics$sample_type,
                                            "DMSO_normalization|dmso","DMSO")
  median_metrics$sample_type <- str_replace(median_metrics$sample_type,
                                            "Staurosporine|STAUR|Stauro|STS|Stouro|ML_STS","STAURO")
  median_metrics$sample_type <- str_replace(median_metrics$sample_type,
                                            "Media_Only|media","Media")
  median_metrics$sample_type <- str_replace(median_metrics$sample_type,
                                            "untreated","Untreated")

  # rename sample_type
  median_metrics$sample_type <- paste0(median_metrics$sample_type, "_", metrics)
  # Transpose the data frame with sample_type as column names
  median_metrics <- median_metrics %>%
    column_to_rownames(var = "sample_type") %>%
    t() %>%
    as.data.frame()
  # Add columns
  median_metrics$Model_type <- Model_type
  median_metrics$Sample <- unique(metadata$orig.ident)
  # get Sample column as 1st column
  rownames(median_metrics) <- NULL
  median_metrics <- median_metrics %>%
    select(Sample, everything())

  # For each cell type if Cell_type column contains multiple cell types
  if (length(unique(metadata$Cell_type))==1) {
    median_metrics$Cell_type <- unique(metadata$Cell_type)
  } else {
    median_metrics$Cell_type <- paste(unique(metadata$Cell_type), collapse = ", ")
  }
  return (median_metrics)
}
