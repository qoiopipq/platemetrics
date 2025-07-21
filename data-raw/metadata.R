## code to prepare `metadata` dataset goes here

usethis::use_data(metadata, overwrite = TRUE)

metadata <- read.csv("/Users/liuxin/Documents/MACseq/macpie/macpieData_backup/extdata/PMMSq033/PMMSq033_metadata_drugnames.csv")
metadata$X <- NULL
# read full PMMSq033 data
project_rawdata <- paste0("/Users/liuxin/Documents/MACseq/macpie/macpieData_backup/extdata/PMMSq033/raw_matrix")
project_name <- "PMMSq033"
raw_counts <- Read10X(data.dir = project_rawdata)
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 1,
                          min.features = 1)
mac <- mac %>%
  inner_join(metadata, by = c(".cell" = "Barcode"))

# Add unique identifier
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

# Filter by read count per sample group
mac <- filter_genes_by_expression(mac,
                                  group_by = "combined_id",
                                  min_counts = 10,
                                  min_samples = 2)
metadata <- mac@meta.data

usethis::use_data(mac, overwrite = TRUE)
usethis::use_data(metadata, overwrite = TRUE)
