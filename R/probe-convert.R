#' Convert rice microarray probe IDs to annotated targets
#'
#' `probe_convert()` converts probe IDs stored in the row names of an
#' expression matrix to targets from the annotation bundled with OryzaProbe.
#' It preserves input order, records an auditable row-level conversion report,
#' and can merge probes that resolve to the same target.
#'
#' @param x A numeric matrix or data frame. Row names must contain probe IDs.
#' @param platform One of `"auto"`, `"GPL6864"`, `"GPL2025"`, or
#'   `"GPL8852"`. Automatic detection uses all row names, not one sentinel
#'   probe. Supplying the platform explicitly is recommended for production
#'   analyses.
#' @param merge_by One of `"none"`, `"mean"`, `"max"`, or `"min"`.
#' @param unmapped How probes without a target should be handled: keep the
#'   original probe ID, drop the row, or raise an error.
#' @param non_rap How annotations that are not canonical RAP-DB gene IDs
#'   should be handled. Controls, miRNA annotations, and other identifiers can
#'   occur in the source platforms.
#' @param multi_mapping How targets containing multiple IDs separated by
#'   `///` should be handled: keep the combined label, use the first target,
#'   expand the expression row once per target, drop it, or raise an error.
#' @param na_rm Whether missing values should be removed during aggregation.
#'   If every value in an aggregation group is missing, the result remains
#'   `NA`.
#' @param quiet If `TRUE`, suppress the conversion summary message.
#'
#' @return A numeric matrix. Its row names contain converted target IDs and
#'   are always unique. Use [probe_conversion_report()] to inspect how every
#'   input row was handled.
#'
#' @examples
#' expression <- matrix(
#'   c(2, 4, 6, 8),
#'   nrow = 2,
#'   dimnames = list(c("36144", "37310"), c("sample_1", "sample_2"))
#' )
#'
#' converted <- probe_convert(
#'   expression,
#'   platform = "GPL6864",
#'   merge_by = "mean"
#' )
#' converted
#' probe_conversion_report(converted)
#'
#' @export
probe_convert <- function(
    x,
    platform = "auto",
    merge_by = "none",
    unmapped = "keep",
    non_rap = "keep",
    multi_mapping = "keep",
    na_rm = FALSE,
    quiet = FALSE) {
  validated <- .validate_expression(x)
  values <- validated$values
  probe_ids <- validated$probe_ids

  platform <- .resolve_platform(platform, probe_ids)
  merge_by <- .match_option(merge_by, c("none", "mean", "max", "min"), "merge_by")
  unmapped <- .match_option(unmapped, c("keep", "drop", "error"), "unmapped")
  non_rap <- .match_option(non_rap, c("keep", "drop", "error"), "non_rap")
  multi_mapping <- .match_option(
    multi_mapping,
    c("keep", "first", "expand", "drop", "error"),
    "multi_mapping"
  )
  .validate_flag(na_rm, "na_rm")
  .validate_flag(quiet, "quiet")

  initial <- .map_probe_ids(probe_ids, platform)
  report <- .resolve_mapping(
    initial,
    unmapped = unmapped,
    non_rap = non_rap,
    multi_mapping = multi_mapping
  )

  kept <- which(report$kept)
  if (length(kept) == 0L) {
    output <- matrix(
      numeric(0),
      nrow = 0L,
      ncol = ncol(values),
      dimnames = list(character(), colnames(values))
    )
  } else {
    kept_report <- report[kept, , drop = FALSE]
    kept_values <- values[kept_report$input_index, , drop = FALSE]

    if (identical(merge_by, "none")) {
      output_ids <- make.unique(kept_report$resolved_target_id, sep = "__probe")
      rownames(kept_values) <- output_ids
      output <- kept_values
      report$output_id[kept] <- output_ids
    } else {
      output <- .merge_expression(
        kept_values,
        groups = kept_report$resolved_target_id,
        method = merge_by,
        na_rm = na_rm
      )
      report$output_id[kept] <- kept_report$resolved_target_id
    }
  }

  summary <- .conversion_summary(
    initial = initial,
    report = report,
    platform = platform,
    output_rows = nrow(output),
    merge_by = merge_by
  )
  attr(output, "oryza_probe_platform") <- platform
  attr(output, "oryza_probe_report") <- report
  attr(output, "oryza_probe_summary") <- summary
  class(output) <- c("oryza_probe_matrix", "matrix", "array")

  if (!quiet) {
    message(.format_summary(summary))
  }

  output
}


#' @export
print.oryza_probe_matrix <- function(x, ...) {
  display <- x
  attr(display, "oryza_probe_platform") <- NULL
  attr(display, "oryza_probe_report") <- NULL
  attr(display, "oryza_probe_summary") <- NULL
  class(display) <- c("matrix", "array")
  print(display, ...)
  invisible(x)
}


#' Backward-compatible probe conversion interface
#'
#' `probeConvert()` preserves the original OryzaProbe function name and its
#' first three arguments. New code should prefer [probe_convert()], whose
#' options make annotation edge cases explicit.
#'
#' @param exprMatrix A numeric matrix or data frame with probe IDs as row
#'   names.
#' @param probeMerge Whether probes resolving to the same target should be
#'   merged.
#' @param mergeBy Aggregation method used when `probeMerge = TRUE`.
#' @param platform Platform name or `"auto"`.
#' @param unmapped How unmapped probes should be handled.
#' @param na.rm Whether missing values should be removed during aggregation.
#' @param quiet Whether the conversion summary should be suppressed.
#' @param ... Additional arguments passed to [probe_convert()], such as
#'   `non_rap` and `multi_mapping`.
#'
#' @return A numeric matrix; see [probe_convert()].
#'
#' @examples
#' expression <- matrix(
#'   1:4,
#'   nrow = 2,
#'   dimnames = list(c("13", "12"), c("sample_1", "sample_2"))
#' )
#' probeConvert(expression, platform = "GPL6864")
#'
#' @export
probeConvert <- function(
    exprMatrix,
    probeMerge = FALSE,
    mergeBy = "mean",
    platform = "auto",
    unmapped = "keep",
    na.rm = FALSE,
    quiet = FALSE,
    ...) {
  .validate_flag(probeMerge, "probeMerge")
  merge_by <- if (probeMerge) mergeBy else "none"

  probe_convert(
    x = exprMatrix,
    platform = platform,
    merge_by = merge_by,
    unmapped = unmapped,
    na_rm = na.rm,
    quiet = quiet,
    ...
  )
}


#' Inspect a probe conversion
#'
#' @param x A matrix returned by [probe_convert()] or [probeConvert()].
#'
#' @return A data frame with the source probe, original annotation, resolved
#'   target, status, inclusion flag, and output row ID.
#'
#' @export
probe_conversion_report <- function(x) {
  report <- attr(x, "oryza_probe_report", exact = TRUE)
  if (is.null(report)) {
    stop("`x` does not contain an OryzaProbe conversion report.", call. = FALSE)
  }
  report
}


#' Summarize a probe conversion
#'
#' @param x A matrix returned by [probe_convert()] or [probeConvert()].
#'
#' @return A named list containing platform, input and output row counts,
#'   mapping-status counts, dropped rows, and aggregation method.
#'
#' @export
probe_conversion_summary <- function(x) {
  summary <- attr(x, "oryza_probe_summary", exact = TRUE)
  if (is.null(summary)) {
    stop("`x` does not contain an OryzaProbe conversion summary.", call. = FALSE)
  }
  summary
}


.validate_expression <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(x) == 0L) {
    stop("`x` must contain at least one probe row.", call. = FALSE)
  }
  if (ncol(x) == 0L) {
    stop("`x` must contain at least one expression column.", call. = FALSE)
  }

  probe_ids <- rownames(x)
  if (is.null(probe_ids) || anyNA(probe_ids) || any(!nzchar(probe_ids))) {
    stop("Every row of `x` must have a non-empty probe ID.", call. = FALSE)
  }

  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("All expression values in `x` must be numeric.", call. = FALSE)
    }
  } else {
    numeric_columns <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      bad <- paste(names(x)[!numeric_columns], collapse = ", ")
      stop(
        sprintf("All expression columns must be numeric; check: %s.", bad),
        call. = FALSE
      )
    }
  }

  values <- as.matrix(x)
  storage.mode(values) <- "double"
  list(values = values, probe_ids = as.character(probe_ids))
}


.validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be one non-missing logical value.", name), call. = FALSE)
  }
  invisible(x)
}


.match_option <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be one character value.", name), call. = FALSE)
  }
  value <- tolower(x)
  if (!value %in% choices) {
    stop(
      sprintf("`%s` must be one of: %s.", name, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
  value
}


.resolve_platform <- function(platform, probe_ids) {
  if (!is.character(platform) || length(platform) != 1L || is.na(platform)) {
    stop("`platform` must be one character value.", call. = FALSE)
  }

  normalized <- toupper(platform)
  if (!identical(normalized, "AUTO")) {
    if (!normalized %in% supported_platforms()) {
      stop(
        sprintf(
          "Unsupported platform `%s`. Supported platforms are: %s.",
          platform,
          paste(supported_platforms(), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(normalized)
  }

  scores <- .platform_scores(probe_ids)
  ordered <- order(scores, decreasing = TRUE)
  best <- scores[ordered[1L]]
  second <- scores[ordered[2L]]

  if (best < 0.5) {
    stop(
      paste0(
        "The platform could not be detected: fewer than half of the probe IDs ",
        "matched any supported platform. Supply `platform` explicitly."
      ),
      call. = FALSE
    )
  }
  if (second >= 0.5 && (best - second) < 0.1) {
    stop(
      "Platform detection was ambiguous. Supply `platform` explicitly.",
      call. = FALSE
    )
  }

  names(scores)[ordered[1L]]
}


.platform_scores <- function(probe_ids) {
  gpl6864_ids <- as.character(GPL6864$ID)
  gpl2025_ids <- as.character(GPL2025$ID)
  gpl8852_pattern <- "^Os[0-9]{2}g[0-9]{7}\\|"

  c(
    GPL6864 = mean(probe_ids %in% gpl6864_ids),
    GPL2025 = mean(probe_ids %in% gpl2025_ids),
    GPL8852 = mean(grepl(gpl8852_pattern, probe_ids))
  )
}


.map_probe_ids <- function(probe_ids, platform) {
  if (identical(platform, "GPL8852")) {
    composite <- grepl("|", probe_ids, fixed = TRUE)
    canonical <- grepl("^Os[0-9]{2}g[0-9]{7}$", probe_ids)
    target <- rep(NA_character_, length(probe_ids))
    target[composite] <- sub("\\|.*$", "", probe_ids[composite])
    target[!composite & canonical] <- probe_ids[!composite & canonical]
  } else {
    annotation <- .annotation_table(platform)
    index <- match(probe_ids, annotation$probe_id)
    target <- annotation$target_id[index]
  }

  data.frame(
    input_index = seq_along(probe_ids),
    probe_id = probe_ids,
    target_id = target,
    status = .mapping_status(probe_ids, target),
    stringsAsFactors = FALSE
  )
}


.mapping_status <- function(probe_ids, target_ids) {
  status <- rep("non_rap", length(target_ids))
  missing <- is.na(target_ids) | !nzchar(target_ids)
  multi <- !missing & grepl("///", target_ids, fixed = TRUE)
  canonical <- !missing & grepl("^Os[0-9]{2}g[0-9]{7}$", target_ids)
  unchanged <- !missing & target_ids == probe_ids

  status[missing | (unchanged & !canonical)] <- "unmapped"
  status[canonical] <- "mapped_rap"
  status[multi] <- "multi_mapped"
  status
}


.resolve_mapping <- function(initial, unmapped, non_rap, multi_mapping) {
  rows <- vector("list", nrow(initial))

  for (i in seq_len(nrow(initial))) {
    current <- initial[i, , drop = FALSE]
    targets <- current$target_id
    resolution <- "as_annotated"
    kept <- TRUE

    if (identical(current$status, "multi_mapped")) {
      split_targets <- trimws(strsplit(targets, "///", fixed = TRUE)[[1L]])
      split_targets <- split_targets[nzchar(split_targets)]

      if (identical(multi_mapping, "error")) {
        stop(
          sprintf("Probe `%s` maps to multiple targets.", current$probe_id),
          call. = FALSE
        )
      } else if (identical(multi_mapping, "drop")) {
        kept <- FALSE
        resolution <- "dropped_multi_mapping"
      } else if (identical(multi_mapping, "first")) {
        targets <- split_targets[1L]
        resolution <- "first_target"
      } else if (identical(multi_mapping, "expand")) {
        targets <- split_targets
        resolution <- "expanded"
      } else {
        resolution <- "combined_target"
      }
    }

    if (identical(current$status, "unmapped")) {
      if (identical(unmapped, "error")) {
        stop(
          sprintf("Probe `%s` has no target annotation.", current$probe_id),
          call. = FALSE
        )
      } else if (identical(unmapped, "drop")) {
        kept <- FALSE
        resolution <- "dropped_unmapped"
      } else {
        targets <- current$probe_id
        resolution <- "kept_probe_id"
      }
    }

    if (identical(current$status, "non_rap")) {
      if (identical(non_rap, "error")) {
        stop(
          sprintf(
            "Probe `%s` maps to non-RAP annotation `%s`.",
            current$probe_id,
            current$target_id
          ),
          call. = FALSE
        )
      } else if (identical(non_rap, "drop")) {
        kept <- FALSE
        resolution <- "dropped_non_rap"
      }
    }

    if (length(targets) == 0L || anyNA(targets)) {
      targets <- current$probe_id
    }

    rows[[i]] <- data.frame(
      input_index = rep(current$input_index, length(targets)),
      probe_id = rep(current$probe_id, length(targets)),
      target_id = rep(current$target_id, length(targets)),
      resolved_target_id = targets,
      status = rep(current$status, length(targets)),
      resolution = rep(resolution, length(targets)),
      kept = rep(kept, length(targets)),
      output_id = rep(NA_character_, length(targets)),
      stringsAsFactors = FALSE
    )
  }

  report <- do.call(rbind, rows)
  rownames(report) <- NULL
  report
}


.merge_expression <- function(values, groups, method, na_rm) {
  group_order <- unique(groups)
  output <- matrix(
    NA_real_,
    nrow = length(group_order),
    ncol = ncol(values),
    dimnames = list(group_order, colnames(values))
  )

  reduce_one <- function(z) {
    if (na_rm && all(is.na(z))) {
      return(NA_real_)
    }
    switch(
      method,
      mean = mean(z, na.rm = na_rm),
      max = max(z, na.rm = na_rm),
      min = min(z, na.rm = na_rm)
    )
  }

  for (i in seq_along(group_order)) {
    block <- values[groups == group_order[i], , drop = FALSE]
    output[i, ] <- vapply(
      seq_len(ncol(block)),
      function(j) reduce_one(block[, j]),
      numeric(1)
    )
  }

  output
}


.conversion_summary <- function(initial, report, platform, output_rows, merge_by) {
  counts <- table(factor(
    initial$status,
    levels = c("mapped_rap", "multi_mapped", "unmapped", "non_rap")
  ))
  kept_input <- unique(report$input_index[report$kept])

  list(
    platform = platform,
    input_rows = nrow(initial),
    output_rows = output_rows,
    mapped_rap = unname(counts[["mapped_rap"]]),
    multi_mapped = unname(counts[["multi_mapped"]]),
    unmapped = unname(counts[["unmapped"]]),
    non_rap = unname(counts[["non_rap"]]),
    dropped_input_rows = nrow(initial) - length(kept_input),
    merge_by = merge_by
  )
}


.format_summary <- function(x) {
  sprintf(
    paste0(
      "%s: %d input row(s) -> %d output row(s); ",
      "%d RAP, %d multi-mapped, %d unmapped, %d non-RAP, %d dropped."
    ),
    x$platform,
    x$input_rows,
    x$output_rows,
    x$mapped_rap,
    x$multi_mapped,
    x$unmapped,
    x$non_rap,
    x$dropped_input_rows
  )
}
