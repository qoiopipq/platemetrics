#' Get a list of files matching a pattern in a directory
#' @param path The directory to search in
#' @param pattern The pattern to match files against
#' @param date The date to filter files by modification time,
#' only getting files with mtime > date
#' @import dplyr
#' @import stringr
#' @return A character vector of file paths that match the pattern and are modified after the specified date
#' @examples
#' path <- "/Volumes/bioinf/team_folders/MGC_VCFG/macseq/results"
#' pattern <- "multiqc_general_stats.txt"
#' date <- "2025-03-07 00:00:00"
#' recent_files <- get_recent_files(path = path,
#'                             pattern = pattern,
#'                              date = date)


get_recent_files <- function(path,
                             pattern = "multiqc_general_stats.txt",
                             date = Sys.time()) {
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
  countmatrix_dirs<-list.files(path,recursive = TRUE,full.names=TRUE,pattern=pattern)
  info <- file.info(countmatrix_dirs)
  cutoff <- as.POSIXct(date)
  recent_files <- rownames(info)[info$mtime > cutoff]
  return(recent_files)
}

#' Read and concatenate MultiQC general stats files
#' @param recent_files A character vector of file paths from `get_recent_files`
#' @param rename_map A named character vector for renaming samples
#' @import stringr
#' @import dplyr

# recent_files<-recent_files[grep("PMMSq040_174132776|PMMSq029_1747700619|PMMSq040_1747632181|AK_PMMSq041_1744279814|PMMSq042_1744692226|PMMSq033_1747610300|PMMSq033_1747633313|6199_6199_1749708215", recent_files, invert = TRUE)]

format_multiqc_stats <- function(recent_files,
                                 rename_map = NULL){
  # Check if recent_files is a character vector
  if (!is.character(recent_files) || length(recent_files) == 0) {
    stop("recent_files must be a non-empty character vector.")
  }
  # Check if rename_map is NULL or character vector
  if (!is.null(rename_map) && !is.character(rename_map)) {
    stop("rename_map must be a character vector or NULL.")
  }
  df_list <- lapply(recent_files, \(f) read.csv(f, sep = "\t", stringsAsFactors = FALSE))
  qc_meta <- bind_rows(df_list)
  # Rename samples by removing "-null" suffix
  qc_meta$Sample<-str_replace(qc_meta$Sample,"-null","")
  # Rename samples based on rename_map if provided
  if(!is.null(rename_map)){
    if (!all(names(rename_map) %in% qc_meta$Sample)) {
      stop("Some sample names in rename_map do not match the samples in qc_meta.")
    } else {
      qc_meta$Sample <- str_replace_all(qc_meta$Sample, rename_map)
    }
  }
  # Get R2 sequencing QC stats from each sample
  r2_stats <- qc_meta %>%
    filter(str_detect(Sample, "[0-9]_2$")) %>%
    transmute(Sample = str_replace(Sample, "_2$", ""),
              duplication_percent = FastQC_mqc.generalstats.fastqc.percent_duplicates,
              total_reads = FastQC_mqc.generalstats.fastqc.total_sequences)
  # Get rows that Sample is not _1 or _2
  # Join with r2_stats
  sample_stats <- qc_meta %>%
    filter(!str_detect(Sample, "_[12]$")) %>%
    transmute(Sample = Sample,
              uniquely_mapped_percent = STAR_mqc.generalstats.star.uniquely_mapped_percent) %>%
    left_join(r2_stats, by = "Sample")
  return(sample_stats)
}
