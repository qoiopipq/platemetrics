#' Compute median metrics by control groups
#' @param metadata The metadata data frame with nCount_RNA, nFeature_RNA, orig.ident columns, and metadata information.
#' This metadata should be able to get from Seurat object.
#' @param metrics The metrics to compute median for. It can be either "nFeature_RNA" or "nCount_RNA".
#' @param group_by The column name in metadata to group by. Default is "Treatment_1".
#' @param controls A character vector of control group names to filter by.
#' Default is c("^DMSO|^Media|^Staurosporine|^Untreated").
#' @param Model_type The model type of the plate. It can be either "2D_non_adherent", "2D_adherent", or "3D".
#'@import dplyr
#'@import tibble
#'


get_plate_median_metrics <- function (metadata = NULL,
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
  # Check model type if it is in any of these three types
  if (is.null(Model_type) | !Model_type %in% c("2D_non_adherent", "2D_adherent", "3D")){
    stop("Model_type must be at least one of '2D_non_adherent', '2D_adherent', or '3D'.")
  }
  # Check orig.ident column
  if(is.null(metadata$orig.ident)){
    stop("Metadata must contain a column named 'orig.ident'.")
  }
  # For entire plate
  # Overall
  overall_median <- data.frame(sample_type = "Overall",
                               median_val = median(metadata[[metrics]], na.rm = TRUE))

  colnames(overall_median) <- c("sample_type",  paste0("median_", metrics))
  # Rename the columns based on metrics
  if (metrics == "nFeature_RNA"){
    colnames(overall_median) <- c("sample_type",  "median_genes")
  } else {
    colnames(overall_median) <- c("sample_type",  "median_UMIs")
  }
  # Filter metadata for controls
  control_metadata <- metadata[grepl(paste(controls, collapse = "|"), metadata[[group_by]]), ]
  control_median <- control_metadata %>%
    group_by(.data[[group_by]]) %>%
    summarise(median_gene_numbers = median(.data[[metrics]]))
  # Rename the columns based on metrics
  if (metrics == "nFeature_RNA"){
    colnames(control_median) <- c("sample_type",  "median_genes")
  } else {
    colnames(control_median) <- c("sample_type",  "median_UMIs")
  }
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
  median_metrics$sample_type <- str_replace(median_metrics$sample_type,
                                            "mock","Mock")
  # rename sample_type with either median_genes or median_UMIs
  if (colnames(median_metrics)[2]%in%c("median_genes", "median_UMIs")) {
    colname_suffix <- colnames(median_metrics)[2]
  } else {
    stop("Column name for median metrics is not recognized. It should be either 'median_genes' or 'median_UMIs'.")
  }
  median_metrics$sample_type <- paste0(median_metrics$sample_type, "_", colname_suffix)
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
    select(Sample, Model_type,everything())
  return (median_metrics)
}



#' Get median metrics for all cell types in a plate
#' @param metadata The metadata data frame with nCount_RNA, nFeature_RNA, orig.ident columns, and metadata information.
#' This metadata should be able to get from Seurat object.
#' @param metrics The metrics to compute median for. It can be either "nFeature_RNA" or "nCount_RNA".
#' @param group_by The column name in metadata to group by. Default is "Treatment_1".
#' @param controls A character vector of control group names to filter by.
#' Default is c("^DMSO|^Media|^Staurosporine|^Untreated").
#' @param Model_type The model type of the plate. It can be either "2D_non_adherent", "2D_adherent", or "3D".
#' @import dplyr
#' @import tibble

get_celltypes_medain_metrics <-  function (metadata = NULL,
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
  # Check model type if it is in any of these three types
  # Check model type can be a vector with length >= 1
  if (is.null(Model_type) | !all(Model_type %in% c("2D_non_adherent", "2D_adherent", "3D"))){
    stop("Model_type must be at least one of '2D_non_adherent', '2D_adherent', or '3D'.")
  }
  # Check orig.ident column
  if(is.null(metadata$orig.ident)){
    stop("Metadata must contain a column named 'orig.ident'.")
  }
  # filter cell type if it is empty
  cell_type_metadata <- metadata %>%
    filter(!is.na(Cell_type) | Cell_type != "" | Cell_type != "EMPTY")
  if(length(unique(cell_type_metadata$Cell_type)) == 1){
    cell_type_metadata <- get_plate_median_metrics(metadata = metadata,
                                                   metrics = metrics,
                                                   group_by = group_by,
                                                   controls = controls,
                                                   Model_type = Model_type)
    cell_type_metadata$Cell_type <- unique(metadata$Cell_type)
    # Order columns
    cell_type_metadata<- cell_type_metadata%>%
      select(Sample, Cell_type, Model_type, everything())
    return(cell_type_metadata)
  }else{
    # Overall metrics for each cell type
    cell_type_overall <- cell_type_metadata %>%
      group_by(Cell_type) %>%
      summarise(median_val = median(.data[[metrics]], na.rm = TRUE)) %>%
      ungroup()
    # Add Overall as Treatment_1
    cell_type_overall$Treatment_1 <- "Overall"
    # Get median metrics for each cell type by control groups
    # Filter metadata for controls
    cell_type_control_metadata <- cell_type_metadata[grepl(paste(controls, collapse = "|"), metadata[[group_by]]), ]
    cell_type_control_median <- cell_type_control_metadata %>%
      group_by(.data[[group_by]], .data$Cell_type)%>%
      summarise(median_val = median(.data[[metrics]]), .groups = "keep")
    ct_median_metrics <- rbind(cell_type_overall, cell_type_control_median)
    ct_median_metrics$Sample <- unique(metadata$orig.ident)
    # Standardise names of controls
    ct_median_metrics$Treatment_1 <- str_replace(ct_median_metrics$Treatment_1,
                                              "DMSO_normalization|dmso","DMSO")
    ct_median_metrics$Treatment_1 <- str_replace(ct_median_metrics$Treatment_1,
                                              "Staurosporine|STAUR|Stauro|STS|Stouro|ML_STS","STAURO")
    ct_median_metrics$Treatment_1 <- str_replace(ct_median_metrics$Treatment_1,
                                              "Media_Only|media","Media")
    ct_median_metrics$Treatment_1 <- str_replace(ct_median_metrics$Treatment_1,
                                              "untreated","Untreated")
    ct_median_metrics$Treatment_1 <- str_replace(ct_median_metrics$Treatment_1,
                                              "mock","Mock")
    ct_median_metrics$sample_type <- ct_median_metrics$Treatment_1
    ct_median_metrics <- ct_median_metrics %>% select(-Treatment_1)
    # Rename the columns based on metrics
    if (metrics == "nFeature_RNA"){
      ct_median_metrics$sample_type <- paste0(ct_median_metrics$sample_type, "_median_genes")
    } else {
      ct_median_metrics$sample_type <- paste0(ct_median_metrics$sample_type, "_median_UMIs")
    }
    # Transpose
    ct_median_metricsw <- ct_median_metrics %>% pivot_wider(names_from = "sample_type",
                                                            values_from = "median_val")
    # Add Model_type
    ct_median_metricsw$Model_type <- Model_type
    # Order columns
    ct_median_metricsw <- ct_median_metricsw %>%
      select(Sample, Cell_type, Model_type, everything())
    ct_median_metricsw <- as.data.frame(ct_median_metricsw)
  return(ct_median_metricsw)
}
}

