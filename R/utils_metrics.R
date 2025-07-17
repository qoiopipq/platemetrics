#' Compute median metrics by control groups
#'
#'@import dplyr


get_median_metrics <- function (metadata,
                                metrics = "nFeature_RNA",
                                group_by = "Treatment_1",
                                controls = c("^DMSO|^Media|^Staurosporine|^Untreated")){
  # Check if metadata is a data frame
  if (!is.data.frame(metadata)) {
    stop("Metadata must be a data frame.")
  }
  # Check if metrics is a character vector
  if (!is.character(metrics)) {
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

  # Overall
  overall_median <- data.frame(sample_type = "Overall",
                               median_val = median(metadata[[metrics]], na.rm = TRUE))
  colnames(overall_median) <- c("sample_type", paste0("median_", metrics))
  # Filter metadata for controls
  control_metadata <- metadata[grepl(paste(controls, collapse = "|"), metadata[[group_by]]), ]
  control_median <- control_metadata %>%
    group_by(!!sym(group_by)) %>%
    summarise(median_gene_numbers = median(!!sym(metrics)))
  colnames(control_median) <- c("sample_type",  paste0("median_", metrics))

  median_metrics <- rbind(overall_median, control_median)
  return (median_metrics)
}
