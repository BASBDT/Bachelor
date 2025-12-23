#' Edge labels for ggparty plots
#'
#' @description
#' Formats the default ggparty edge labels to a more readable mathematical
#' notation for categorical and numeric splits.
#'
#' @param ggpy A ggparty object created by ggparty(PartyTree).
#' @param digits Optional integer; number of digits to round numeric thresholds.
#' @return A list of language objects with formatted edge labels.
#' @keywords internal
get_edge_labels <- function(ggpy, digits = NULL) {
  # Default edge labels from ggparty
  breaks_labels <- ggpy$data$breaks_label

  # Differentiate between numeric and categorical splits
  numericSplit_indices <- grep("NA", breaks_labels)
  allIndices <- seq_along(breaks_labels)
  # Index 1 is an empty placeholder used by ggparty; drop safely if present
  if (length(numericSplit_indices) > 0L) {
    numericSplit_indices <- setdiff(numericSplit_indices, 1L)
  }
  categoricalSplit_indices <- setdiff(allIndices, numericSplit_indices)
  categorical_labels <- breaks_labels[categoricalSplit_indices]

  # Trim categorical label vectors to avoid overly long strings
  trimLabel <- function(x) {
    labelString <- paste(x, collapse = " ")
    changed <- FALSE
    while (nchar(labelString) > 30) {
      x <- x[-length(x)]
      labelString <- paste(x, collapse = " ")
      changed <- TRUE
    }
    if (changed) x[length(x) + 1] <- "..."
    x
  }

  categorical_labels <- lapply(categorical_labels, trimLabel)
  numeric_labels <- breaks_labels[numericSplit_indices]

  # Categorical label formatting
  if (as.numeric(R.version$major) >= 4 && as.numeric(R.version$minor) >= 2) {
    new_categorical_labels <- lapply(
      categorical_labels,
      function(x) bquote("x \u1455" ~ group("{", .(paste(x, collapse = ", ")), "}"))
    ) # \u1455 is a temporary alternative for subset sign rendering issues in some R 4.2.x
  } else {
    new_categorical_labels <- lapply(
      categorical_labels,
      function(x) bquote('x' %subset% group('{', .(paste(x, collapse = ', ')), '}'))
    )
  }
  breaks_labels[categoricalSplit_indices] <- new_categorical_labels

  # Numeric label formatting (override partykit defaults)
  numericBreaksFormatting <- function(x, digits) {
    operator <- gsub("[^<>=]+", "", x)
    value <- gsub("[^0-9.]", "", x)
    if (!is.null(digits)) value <- round(as.numeric(value), digits)
    switch(
      operator,
      "<=" = bquote("x \U2264" ~ .(value)),  # ≤
      "<"  = bquote("x < " ~ .(value)),
      ">=" = bquote("x \U2265" ~ .(value)),  # ≥
      ">"  = bquote("x >" ~ .(value)),
      stop("Invalid operator")
    )
  }
  new_numeric_labels <- lapply(numeric_labels, function(x) numericBreaksFormatting(x, digits))
  breaks_labels[numericSplit_indices] <- new_numeric_labels

  breaks_labels
}

#' Deprecated alias for get_edge_labels
#' @keywords internal
getEdgeLabels <- function(ggpy, digits = NULL) {
  .Deprecated("get_edge_labels")
  get_edge_labels(ggpy = ggpy, digits = digits)
}
  