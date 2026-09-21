# OryzaProbe 0.2.0

- Rewrote probe conversion to preserve the input expression-row order.
- Fixed probe merging for GPL6864, GPL2025, and GPL8852.
- Replaced sentinel-probe platform detection with whole-input scoring.
- Added explicit policies for unmapped, non-RAP, and multi-target annotations.
- Added row-level conversion reports and summary metadata.
- Added public annotation and provenance inspection functions.
- Preserved `probeConvert()` as a backward-compatible interface.
- Replaced the original placeholder tests with regression and unit tests.
- Replaced the disabled vignette with a self-contained executable example.
- Made the legacy annotation tables inspectable and reproducible from
  `data-raw/`.

## Breaking changes

- Conversion failures now raise errors instead of returning `NULL`.
- Conversion functions now return a numeric matrix consistently.
- Duplicate unmerged target IDs use the explicit `__probe` suffix.
