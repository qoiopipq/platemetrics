#' Read in CSV files from each plate or
#'
#' @param path The directory path where the CSV files are located.
#' @param pattern The pattern to match files against, default is "*median_genes.csv".
#' It can be only one of "\*median_genes.csv" or "\*median_counts.csv".
#' @param date The date to filter files by modification time, default is the current system time.
#' @param plate Logical value indicating whether to read plate-level data (TRUE) or cell type-level data (FALSE).
#' @return A data frame containing the concatenated data from the CSV files.
#' @import dplyr
#' @import stringr
#' @import readr
# path <- "/Users/liuxin/Documents/MACseq/Projects"
# pattern <- "*median_genes.csv"
# date <- as.POSIXct("2025-03-07 00:00:00")
# test_read_metrics <- read_metrics(path = path,
#                          pattern = pattern,
#                          date = date,
#                          plate = TRUE)


read_metrics <- function(path,
                         pattern = "*median_genes.csv",
                         date = Sys.time(),
                         plate = TRUE){
  # Check path exists
  if (!dir.exists(path)) {
    stop("The specified path does not exist.")
  }
  # Check pattern is a character string
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("The pattern must be a single character string.")
  }
  # Check date is in correct format
  if (!inherits(as.POSIXct(date), "POSIXct")) {
    stop("The date must be in a valid POSIXct format.")
  }
  # Check plate is either TRUE or FALSE
  if (!is.logical(plate) || length(plate) != 1) {
    stop("The plate argument must be a single logical value (TRUE or FALSE).")
  }
  # Get recent files based on the pattern and date
  recent_dirs<-list.files(path = path,recursive = TRUE,full.names=T,pattern = pattern)
  file_info <- file.info(recent_dirs)
  cutoff <- as.POSIXct(date)
  recent_files <- rownames(file_info)[file_info$mtime > cutoff]
  if (length(recent_files) == 0) {
    stop("No files found matching the specified pattern and date.")
  }
  # Read and concatenate CSV files
  # pattern = "*median_genes.csv"
  if (pattern == "*median_genes.csv") {
    if (plate) {
      # plate level
      plate_gene_dirs <- grep("*ct_median_genes.csv", recent_files, value = TRUE, invert = TRUE)
      gene_list <- lapply(plate_gene_dirs, \(f) read.csv(f, sep=","))
      gene_meta <- bind_rows(gene_list)
      #if there is a X Column, remove X column
      if ("X" %in% colnames(gene_meta)) {
        gene_meta <- gene_meta %>% select(-X)
        }
      } else {
        # cell type level
        ct_gene_dirs <- grep("*ct_median_genes.csv", recent_files, value = TRUE)
        gene_list <- lapply(ct_gene_dirs, \(f) read.csv(f, sep=","))
        gene_meta <- bind_rows(gene_list)
        #if there is a X Column, remove X column
        if ("X" %in% colnames(gene_meta)) {
          gene_meta <- gene_meta %>% select(-X)
        }
      }
  } else {
      # this is when pattern == "*median_counts.csv"
      if (plate) {
        # plate level
        plate_gene_dirs <- grep("*ct_median_counts.csv", recent_files, value = TRUE, invert = TRUE)
        gene_list <- lapply(plate_gene_dirs, function(f) {
          df <- read.csv(f, sep = ",")
          colnames(df) <- str_replace_all(colnames(df), "_counts", "_UMIs")
          df
        })
        gene_meta <- bind_rows(gene_list)
        #if there is a X Column, remove X column
        if ("X" %in% colnames(gene_meta)) {
          gene_meta <- gene_meta %>% select(-X)
        }
      } else {
        # cell type level
        ct_gene_dirs <- grep("*ct_median_counts.csv", recent_files, value = TRUE)
        gene_list <- lapply(ct_gene_dirs, function(f) {
          df <- read.csv(f, sep = ",")
          colnames(df) <- str_replace_all(colnames(df), "_counts", "_UMIs")
          df
        })
        gene_meta <- bind_rows(gene_list)
        #if there is a X Column, remove X column
        if ("X" %in% colnames(gene_meta)) {
          gene_meta <- gene_meta %>% select(-X)
        }
      }
  }

  # Reorder columns to have Sample as the first column
  if ("Sample" %in% colnames(gene_meta)) {
    # rename Org_Cell to Model_Type
    # if Model_type is NA, replace it with Org_Cell
    gene_meta$Model_type <- ifelse(is.na(gene_meta$Model_type), gene_meta$Org_Cell, gene_meta$Model_type)
    gene_meta <- gene_meta %>% select(-Org_Cell)
    gene_meta <- gene_meta %>% select(Sample, Model_type, everything())
  } else {
    stop("The 'Sample' column is missing in the data.")
  }
  return(gene_meta)
}


