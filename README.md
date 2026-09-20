# OryzaProbe

OryzaProbe converts probe identifiers from the rice microarray platforms
GPL6864, GPL2025, and GPL8852 to annotated targets used by RAP-DB workflows.

Version 0.2.0 rewrites the conversion engine to preserve expression-row order,
make annotation edge cases visible, and provide an auditable conversion report.

## Installation

Install the development version from GitHub:

```r
remotes::install_github("lym-bioinfo/OryzaProbe")
```

Or install from a local package source directory:

```r
install.packages(".", repos = NULL, type = "source")
```

## Basic use

```r
library(OryzaProbe)

expression <- matrix(
  c(2, 4, 6, 8),
  nrow = 2,
  dimnames = list(c("36144", "37310"), c("sample_1", "sample_2"))
)

converted <- probe_convert(
  expression,
  platform = "GPL6864",
  merge_by = "mean"
)

converted
probe_conversion_report(converted)
probe_conversion_summary(converted)
```

Supplying `platform` explicitly is recommended for reproducible analyses.
`platform = "auto"` compares all probe IDs against the supported annotations
and raises an error when the result is weak or ambiguous.

## Annotation edge cases

The source platform annotations contain controls, miRNA identifiers, unmapped
probes, and a small number of one-to-many mappings. Their treatment is explicit:

```r
converted <- probe_convert(
  expression,
  platform = "GPL6864",
  unmapped = "keep",       # keep, drop, or error
  non_rap = "keep",        # keep, drop, or error
  multi_mapping = "keep"   # keep, first, expand, drop, or error
)
```

The result is always a numeric matrix. The attached report records the source
probe, its annotation status, how it was resolved, and its output row ID.

## Annotation provenance

The bundled GPL6864 and GPL2025 mappings are preserved from OryzaProbe 0.1.0
(the package snapshot dated 2022-10-15). They have not been silently relabeled
as a current RAP-DB release. Use `probe_annotation_metadata()` to record the
snapshot used by an analysis. Rebuild instructions and inspectable mapping
tables live in `data-raw/`.

## Compatibility

The original `probeConvert(exprMatrix, probeMerge, mergeBy)` function remains
available and delegates to the new conversion engine. New analyses should use
`probe_convert()` so all mapping policies are visible in the call.
