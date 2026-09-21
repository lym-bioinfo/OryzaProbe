# Rebuild OryzaProbe's internal annotation data from the audited legacy tables.
#
# Run this script from the package root. The CSV files were recovered from the
# OryzaProbe 0.1.0 `R/sysdata.rda` object. They preserve the published mapping
# exactly; they are not presented as a newly downloaded RAP-DB annotation.

read_mapping <- function(path, integer_id = FALSE) {
  column_classes <- if (integer_id) c("integer", "character") else c("character", "character")
  mapping <- utils::read.csv(
    path,
    colClasses = column_classes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(identical(names(mapping), c("ID", "ACC")))
  stopifnot(!anyNA(mapping$ID), !anyDuplicated(mapping$ID))
  mapping
}

GPL6864 <- read_mapping("data-raw/GPL6864.csv", integer_id = TRUE)
GPL2025 <- read_mapping("data-raw/GPL2025.csv")
GPL8852 <- read_mapping("data-raw/GPL8852.csv")

mapping_metadata <- data.frame(
  platform = c("GPL6864", "GPL2025", "GPL8852"),
  source_url = c(
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6864",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2025",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8852"
  ),
  annotation_snapshot = c(
    "OryzaProbe 0.1.0 legacy table",
    "OryzaProbe 0.1.0 legacy table",
    "Composite probe row name"
  ),
  snapshot_date = c("2022-10-15", "2022-10-15", "not applicable"),
  mapping_method = c(
    "Exact probe-ID lookup",
    "Exact probe-ID lookup",
    "Text before the first pipe character"
  ),
  rows = c(nrow(GPL6864), nrow(GPL2025), nrow(GPL8852)),
  stringsAsFactors = FALSE
)

save(
  GPL6864,
  GPL2025,
  GPL8852,
  mapping_metadata,
  file = "R/sysdata.rda",
  version = 2,
  compress = "xz"
)
