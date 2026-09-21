# Annotation data provenance

`GPL6864.csv` and `GPL2025.csv` are lossless text exports of the annotation
tables shipped inside OryzaProbe 0.1.0. They were recovered during the 0.2.0
rewrite so the package data can be inspected, diffed, and rebuilt.

These files preserve the published 2022 mapping; they are not claimed to be a
current RAP-DB release. Before changing any biological annotation, update the
source and version information, document the matching rules, and review all
unmapped, non-RAP, and multi-target records.

Run the following command from the package root to regenerate `R/sysdata.rda`:

```sh
Rscript --vanilla data-raw/build_internal_data.R
```

The GPL8852 conversion does not use a lookup table. For that platform the
target is encoded before the first `|` character in each composite row name.
