library(OryzaProbe)

# Order preservation: the 0.1.0 implementation silently swapped these labels.
x <- matrix(
  c(13, 12),
  ncol = 1,
  dimnames = list(c("13", "12"), "sample")
)
ordered <- probe_convert(x, platform = "GPL6864", quiet = TRUE)
stopifnot(
  identical(rownames(ordered), c("Os01g0721700", "Os01g0532600")),
  identical(as.numeric(ordered), c(13, 12))
)

# Valid subsets no longer require a hard-coded sentinel probe.
subset_x <- matrix(
  c(13, 14),
  ncol = 1,
  dimnames = list(c("13", "14"), "sample")
)
detected <- probe_convert(subset_x, quiet = TRUE)
stopifnot(identical(attr(detected, "oryza_probe_platform"), "GPL6864"))

# Merge after conversion, including the legacy compatibility function.
merge_x <- matrix(
  c(2, 4),
  ncol = 1,
  dimnames = list(c("36144", "37310"), "sample")
)
merged <- probe_convert(
  merge_x,
  platform = "GPL6864",
  merge_by = "mean",
  quiet = TRUE
)
legacy <- probeConvert(
  merge_x,
  probeMerge = TRUE,
  platform = "GPL6864",
  quiet = TRUE
)
stopifnot(
  identical(rownames(merged), "Os01g0100700"),
  identical(as.numeric(merged), 3),
  identical(as.numeric(legacy), 3)
)

# GPL8852 composite IDs are parsed before aggregation.
gpl8852_x <- matrix(
  c(2, 4),
  ncol = 1,
  dimnames = list(
    c("Os01g0100100|A|B|0", "Os01g0100100|C|D|0"),
    "sample"
  )
)
gpl8852 <- probe_convert(
  gpl8852_x,
  platform = "GPL8852",
  merge_by = "mean",
  quiet = TRUE
)
stopifnot(
  identical(rownames(gpl8852), "Os01g0100100"),
  identical(as.numeric(gpl8852), 3)
)
