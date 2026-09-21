test_that("mapping preserves input order", {
  x <- matrix(
    c(13, 12),
    ncol = 1,
    dimnames = list(c("13", "12"), "sample")
  )

  result <- probe_convert(x, platform = "GPL6864", quiet = TRUE)

  expect_equal(
    rownames(result),
    c("Os01g0721700", "Os01g0532600")
  )
  expect_s3_class(result, "oryza_probe_matrix")
  expect_true(is.matrix(result))
  expect_equal(as.numeric(result), c(13, 12))
})


test_that("automatic platform detection uses the complete subset", {
  x <- matrix(
    c(13, 14),
    ncol = 1,
    dimnames = list(c("13", "14"), "sample")
  )

  result <- probe_convert(x, quiet = TRUE)

  expect_equal(attr(result, "oryza_probe_platform"), "GPL6864")
  expect_equal(rownames(result), c("Os01g0721700", "Os06g0215600"))
})


test_that("multiple probes are merged after target conversion", {
  x <- matrix(
    c(2, 4, 10, 14),
    nrow = 2,
    dimnames = list(c("36144", "37310"), c("sample_1", "sample_2"))
  )

  mean_result <- probe_convert(
    x,
    platform = "GPL6864",
    merge_by = "mean",
    quiet = TRUE
  )
  max_result <- probe_convert(
    x,
    platform = "GPL6864",
    merge_by = "max",
    quiet = TRUE
  )

  expect_equal(dim(mean_result), c(1L, 2L))
  expect_equal(rownames(mean_result), "Os01g0100700")
  expect_equal(as.numeric(mean_result), c(3, 12))
  expect_equal(as.numeric(max_result), c(4, 14))
})


test_that("unmerged duplicate targets receive transparent suffixes", {
  x <- matrix(
    c(2, 4),
    ncol = 1,
    dimnames = list(c("36144", "37310"), "sample")
  )

  result <- probe_convert(x, platform = "GPL6864", quiet = TRUE)

  expect_equal(
    rownames(result),
    c("Os01g0100700", "Os01g0100700__probe1")
  )
})


test_that("unmapped and non-RAP policies are explicit", {
  unmapped_x <- matrix(1, dimnames = list("unknown_probe", "sample"))
  non_rap_x <- matrix(1, dimnames = list("35", "sample"))

  kept <- probe_convert(
    unmapped_x,
    platform = "GPL6864",
    unmapped = "keep",
    quiet = TRUE
  )
  dropped <- probe_convert(
    unmapped_x,
    platform = "GPL6864",
    unmapped = "drop",
    quiet = TRUE
  )

  expect_equal(rownames(kept), "unknown_probe")
  expect_equal(nrow(dropped), 0L)
  expect_error(
    probe_convert(
      unmapped_x,
      platform = "GPL6864",
      unmapped = "error",
      quiet = TRUE
    ),
    "no target annotation"
  )
  expect_error(
    probe_convert(
      non_rap_x,
      platform = "GPL6864",
      non_rap = "error",
      quiet = TRUE
    ),
    "non-RAP annotation"
  )
})


test_that("one-to-many mappings can be kept, selected, expanded, or dropped", {
  x <- matrix(2, dimnames = list("Os.19863.1.S1_a_at", "sample"))

  expanded <- probe_convert(
    x,
    platform = "GPL2025",
    multi_mapping = "expand",
    quiet = TRUE
  )
  first <- probe_convert(
    x,
    platform = "GPL2025",
    multi_mapping = "first",
    quiet = TRUE
  )
  dropped <- probe_convert(
    x,
    platform = "GPL2025",
    multi_mapping = "drop",
    quiet = TRUE
  )

  expect_equal(rownames(expanded), c("Os02g0558100", "Os04g0442700"))
  expect_equal(as.numeric(expanded), c(2, 2))
  expect_equal(rownames(first), "Os02g0558100")
  expect_equal(nrow(dropped), 0L)
})


test_that("GPL8852 targets are parsed and merged", {
  x <- matrix(
    c(2, 4),
    ncol = 1,
    dimnames = list(
      c("Os01g0100100|A|B|0", "Os01g0100100|C|D|0"),
      "sample"
    )
  )

  result <- probe_convert(x, platform = "GPL8852", merge_by = "mean", quiet = TRUE)

  expect_equal(rownames(result), "Os01g0100100")
  expect_equal(as.numeric(result), 3)
})


test_that("NA aggregation is defined", {
  x <- matrix(
    c(NA, NA, 2, NA),
    nrow = 2,
    dimnames = list(c("36144", "37310"), c("sample_1", "sample_2"))
  )

  kept_na <- probe_convert(
    x,
    platform = "GPL6864",
    merge_by = "mean",
    na_rm = FALSE,
    quiet = TRUE
  )
  removed_na <- probe_convert(
    x,
    platform = "GPL6864",
    merge_by = "mean",
    na_rm = TRUE,
    quiet = TRUE
  )

  expect_true(is.na(kept_na[1, 1]))
  expect_true(is.na(removed_na[1, 1]))
  expect_true(is.na(kept_na[1, 2]))
  expect_equal(removed_na[1, 2], 2)
})


test_that("conversion report and summary are auditable", {
  x <- matrix(1, dimnames = list("13", "sample"))
  result <- probe_convert(x, platform = "GPL6864", quiet = TRUE)

  report <- probe_conversion_report(result)
  summary <- probe_conversion_summary(result)

  expect_equal(report$probe_id, "13")
  expect_equal(report$target_id, "Os01g0721700")
  expect_equal(report$status, "mapped_rap")
  expect_equal(summary$platform, "GPL6864")
  expect_equal(summary$input_rows, 1L)
  expect_error(probe_conversion_report(x), "does not contain")
})


test_that("printing hides the potentially large audit attributes", {
  x <- matrix(1, dimnames = list("13", "sample"))
  result <- probe_convert(x, platform = "GPL6864", quiet = TRUE)

  printed <- capture.output(print(result))

  expect_true(any(grepl("Os01g0721700", printed, fixed = TRUE)))
  expect_false(any(grepl("oryza_probe_report", printed, fixed = TRUE)))
})


test_that("legacy interface delegates to the rewritten engine", {
  x <- matrix(
    c(2, 4),
    ncol = 1,
    dimnames = list(c("36144", "37310"), "sample")
  )

  result <- probeConvert(
    x,
    probeMerge = TRUE,
    mergeBy = "mean",
    platform = "GPL6864",
    quiet = TRUE
  )

  expect_equal(rownames(result), "Os01g0100700")
  expect_equal(as.numeric(result), 3)
})


test_that("invalid inputs fail with actionable messages", {
  x <- matrix(1, dimnames = list("13", "sample"))

  expect_error(probe_convert(unname(x)), "probe ID")
  expect_error(probe_convert(as.character(x)), "numeric matrix")
  expect_error(probe_convert(x, platform = "GPL0000"), "Unsupported platform")
  expect_error(probe_convert(x, merge_by = "median"), "merge_by")
  expect_error(probe_convert(x, na_rm = NA), "na_rm")
})


test_that("annotation data and provenance are inspectable", {
  expect_equal(supported_platforms(), c("GPL6864", "GPL2025", "GPL8852"))
  expect_equal(nrow(probe_annotation("GPL6864")), 45151L)
  expect_equal(nrow(probe_annotation("GPL2025")), 57381L)
  expect_equal(nrow(probe_annotation("GPL8852")), 0L)

  metadata <- probe_annotation_metadata()
  expect_equal(metadata$platform, supported_platforms())
  expect_match(metadata$annotation_snapshot[1], "0.1.0")
})
