test_that("Test if get_plate_median_metrics works", {
  data(mac)
  plate_metrics <- get_plate_median_metrics(
    metadata = mac@meta.data,
    metrics = "nFeature_RNA",
    group_by = "Treatment_1",
    controls = c("^DMSO|^Media|^Staurosporine|^Untreated"),
    Model_type = "2D_non_adherent"
  )
  expect_true(is.data.frame(plate_metrics))
  expect_true(!"Cell_type" %in% colnames(plate_metrics))
})


test_that("Test if get_celltypes_median_metrics works with cell type",{
  data(mac)
  cell_type_metrics <- get_celltypes_medain_metrics(
    metadata = mac@meta.data,
    metrics = "nFeature_RNA",
    group_by = "Treatment_1",
    controls = c("^DMSO|^Media|^Staurosporine|^Untreated"),
    Model_type = c("3D","3D","2D_non_adherent")
  )
  expect_true(is.data.frame(cell_type_metrics))
  expect_true("Cell_type" %in% colnames(cell_type_metrics))
})
