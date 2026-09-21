.ORYZA_PLATFORMS <- c("GPL6864", "GPL2025", "GPL8852")


#' List supported rice microarray platforms
#'
#' @return A character vector of supported GEO platform accessions.
#'
#' @export
supported_platforms <- function() {
  .ORYZA_PLATFORMS
}


#' Inspect the bundled probe annotation
#'
#' @param platform One supported GEO platform accession.
#'
#' @return A data frame with `probe_id`, `target_id`, and `status`. GPL8852
#'   returns an empty table because its target is parsed from the composite row
#'   name rather than looked up in a bundled table.
#'
#' @export
probe_annotation <- function(platform) {
  platform <- .explicit_platform(platform)
  annotation <- .annotation_table(platform)

  if (nrow(annotation) == 0L) {
    annotation$status <- character()
  } else {
    annotation$status <- .mapping_status(
      annotation$probe_id,
      annotation$target_id
    )
  }

  attr(annotation, "platform") <- platform
  annotation
}


#' Inspect annotation provenance
#'
#' The bundled GPL6864 and GPL2025 tables are preserved from OryzaProbe 0.1.0.
#' This function makes that legacy provenance explicit so analyses can record
#' the exact annotation snapshot they used.
#'
#' @return A data frame describing the mapping source and method for each
#'   supported platform.
#'
#' @export
probe_annotation_metadata <- function() {
  mapping_metadata
}


.explicit_platform <- function(platform) {
  if (!is.character(platform) || length(platform) != 1L || is.na(platform)) {
    stop("`platform` must be one character value.", call. = FALSE)
  }
  platform <- toupper(platform)
  if (!platform %in% supported_platforms()) {
    stop(
      sprintf(
        "Unsupported platform `%s`. Supported platforms are: %s.",
        platform,
        paste(supported_platforms(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  platform
}


.annotation_table <- function(platform) {
  platform <- .explicit_platform(platform)
  source <- switch(
    platform,
    GPL6864 = GPL6864,
    GPL2025 = GPL2025,
    GPL8852 = GPL8852
  )

  data.frame(
    probe_id = as.character(source$ID),
    target_id = as.character(source$ACC),
    stringsAsFactors = FALSE
  )
}
