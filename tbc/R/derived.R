#' Materialize derived (bivariate) columns used by the model/tree
#'
#' Creates columns only for bivariate splits actually used by the final tree:
#' for each used (i, j, a) it adds the column
#' \deqn{D(i,j,a) = X[, i] + a * X[, j].}
#'
#' @param X data.frame or matrix with the original features. New columns will be
#'   appended; existing columns are never renamed.
#' @param model list or object that carries either:
#'   \itemize{
#'     \item \code{$DerivedSpec} and \code{$is_derived} (from TreeList), or
#'     \item \code{$TreeList} with those fields.
#'   }
#' @param add_pretty logical; add human-readable names as column names using the
#'   \strong{nice} convention \code{"Xi ± c*Xj"} (e.g., \code{"SepalLen + 0.50*PetalLen"}).
#'   This replaces the older \code{"A0.50"} token.
#' @param add_safe logical; also add a machine-safe alias (e.g., \code{"C1_P050_C2"}).
#'   Leave \code{FALSE} unless you need alias columns for programmatic access.
#' @param overwrite logical; if \code{FALSE}, existing columns with the same name
#'   are left unchanged.
#' @return \code{X} with added columns for each unique used (i,j,a) combination.
#' @keywords internal
tbc_materialize_derived <- function(
  X, model,
  add_pretty = TRUE,
  add_safe   = FALSE,
  overwrite  = FALSE
) {
  if (is.null(X)) return(X)

  # Extract TreeList fields from model or pass-through list
  TL <- NULL
  if (is.list(model) && !is.null(model$TreeList)) {
    TL <- model$TreeList
  } else if (is.list(model)) {
    TL <- model
  } else {
    stop("tbc_materialize_derived: 'model' must be a tbc_model or a TreeList.")
  }

  if (is.null(TL$DerivedSpec) || is.null(TL$is_derived)) return(X)

  specs <- TL$DerivedSpec
  flags <- TL$is_derived
  if (length(specs) == 0L || !any(isTRUE(flags))) return(X)

  Xdf   <- if (!is.data.frame(X)) as.data.frame(X, check.names = FALSE) else X
  Names <- TL$VarNames %||% colnames(Xdf) %||% character()

  # unique keying (i|j|a) to avoid recomputation
  uniq <- list(); keyset <- character()
  add_spec <- function(s) {
    if (is.null(s)) return()
    i <- as.integer(s$i); j <- as.integer(s$j); a <- as.numeric(s$a)
    if (!is.finite(a) || is.na(i) || is.na(j)) return()
    k <- paste(i, j, format(a, trim = TRUE), sep = "|")
    if (k %in% keyset) return()
    keyset <<- c(keyset, k)
    uniq[[length(uniq) + 1L]] <<- list(i = i, j = j, a = a)
  }
  for (n in seq_along(specs)) if (isTRUE(flags[[n]])) add_spec(specs[[n]])
  if (length(uniq) == 0L) return(Xdf)

  # Safe alias tag: P for positive, M for negative, 2-digit decimals without dot
  fmt_safe <- function(a) {
    s <- if (a >= 0) "P" else "M"
    v <- sprintf("%0.2f", abs(a))
    v <- gsub("\\.", "", v)  # 0.50 -> 050
    paste0(s, v)
  }

  # Pretty label using the internal helper (ensures the "+"/"-" convention)
  pretty_label <- function(i, j, a) {
    # prefer internal pretty helper if available
    if (exists(".tbc_bivar_pretty", mode = "function")) {
      return(.tbc_bivar_pretty(Names, i, j, a))
    }
    # fallback: local format
    ni <- Names[i] %||% paste0("C", i)
    nj <- Names[j] %||% paste0("C", j)
    sign <- if (a < 0) "-" else "+"
    coef <- format(abs(a), trim = TRUE)
    coef <- sub("0+$", "", sub("(\\.\\d*?)0+$", "\\1", coef))
    paste0(ni, " ", sign, " ", coef, "*", nj)
  }

  for (sp in uniq) {
    i <- sp$i; j <- sp$j; a <- sp$a
    ni <- Names[i] %||% paste0("C", i)
    nj <- Names[j] %||% paste0("C", j)

    Xi <- Xdf[[ni]] %||% Xdf[[i]]
    Xj <- Xdf[[nj]] %||% Xdf[[j]]
    if (is.null(Xi) || is.null(Xj)) next

    val    <- Xi + a * Xj
    pretty <- pretty_label(i, j, a)
    safe   <- paste0("C", i, "_", fmt_safe(a), "_C", j)

    if (add_pretty && (overwrite || is.null(Xdf[[pretty]]))) Xdf[[pretty]] <- val
    if (add_safe   && (overwrite || is.null(Xdf[[safe]])))   Xdf[[safe]]   <- val
  }

  Xdf
}
