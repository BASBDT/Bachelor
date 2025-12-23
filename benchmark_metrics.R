# Benchmark metrics and graphics generator
#
# Reads all CSV files from tools/benchmark, computes metrics
# (MCC, Accuracy, F1-macro/micro, Balanced Accuracy) from flattened
# confusion matrix columns (cm_i_j), and writes augmented CSVs to
# tools/benchmark_metrics. Produces per-dataset graphics in
# tools/benchmark_graphics/<dataset>/.
#
# Main (primary) graphics per dataset:
#   - Scatter: MCC vs Maximum depth (raw points)
#   - Scatter: MCC vs Leaf count (raw points)
#   - Scatter: MCC vs Maximum depth (folds averaged within each repeat)
#   - Scatter: MCC vs Leaf count (folds averaged within each repeat)
#   - Scatter: MCC vs Maximum depth (mean MCC per x value within algo/criterion)
#   - Scatter: MCC vs Leaf count (mean MCC per x value within algo/criterion)
#   - MD plot: MCC (raw, split by algorithm/criterion)
#   - MD plot: MCC (folds averaged within each repeat, split by algorithm/criterion)
#
# All other plots are written into a secondary folder:
#   tools/benchmark_graphics/<dataset>/secondary_plots/
#
# Notes
# - Only rows with mode == "predict" are used for metrics and plots.
# - Confusion matrices of varying size are supported; cm_i_j columns are detected
#   automatically to reconstruct the K x K matrix.
# - Optional MD plots require CRAN package DataVisualizations (Thrun). If not
#   available, MD plots are skipped (no install attempted).

# --- Package handling ---------------------------------------------------------

.ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE, quiet = TRUE)
  }
}

.ensure_packages <- function(pkgs) for (p in pkgs) .ensure_package(p)

# Required packages
.ensure_packages(c("ggplot2", "readr", "dplyr", "stringr"))
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

quietly_load <- function(pkg) try(suppressPackageStartupMessages(library(pkg, character.only = TRUE)), silent = TRUE)
quietly_load("DataVisualizations")

# --- Lightweight logging helpers ---------------------------------------------

.log_ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_info <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message("[INFO  ", .log_ts(), "] ", paste0(..., collapse = ""))
    try(flush.console(), silent = TRUE)
  }
}

log_warn <- function(...) {
  warning(paste0("[WARN  ", .log_ts(), "] ", paste0(..., collapse = "")), call. = FALSE, immediate. = TRUE)
  try(flush.console(), silent = TRUE)
}

log_error <- function(...) {
  message("[ERROR ", .log_ts(), "] ", paste0(..., collapse = ""))
  try(flush.console(), silent = TRUE)
}

# --- Algorithm / criterion labeling + ordering --------------------------------

# Build combined label. Criterion "(none)" when empty.
add_algo_crit <- function(df) {
  if (all(c("algorithm", "criterion") %in% names(df))) {
    df %>%
      mutate(
        criterion = ifelse(is.na(criterion) | criterion == "", "(none)", as.character(criterion)),
        algo_crit  = paste0(as.character(algorithm), " / ", as.character(criterion))
      )
  } else if ("algorithm" %in% names(df) && !"algo_crit" %in% names(df)) {
    df$algo_crit <- as.character(df$algorithm)
    df
  } else {
    df
  }
}

# Alphabetical by algorithm then criterion (stable across datasets).
apply_algo_crit_levels <- function(df) {
  if (!"algo_crit" %in% names(df)) return(df)

  # If algorithm+criterion present, sort pairs. Otherwise sort labels.
  if (all(c("algorithm", "criterion") %in% names(df))) {
    tmp <- df %>%
      transmute(
        algo_crit = as.character(algo_crit),
        algorithm = as.character(algorithm),
        criterion = as.character(criterion)
      ) %>%
      distinct() %>%
      arrange(tolower(algorithm), tolower(criterion), tolower(algo_crit))

    levels <- tmp$algo_crit
  } else {
    levels <- sort(unique(as.character(df$algo_crit)), method = "radix")
  }

  df$algo_crit <- factor(as.character(df$algo_crit), levels = levels)
  df
}

# --- Confusion matrix helpers -------------------------------------------------

detect_cm_spec <- function(df) {
  cm_cols <- grep("^cm_\\d+_\\d+$", names(df), value = TRUE)
  if (length(cm_cols) == 0) return(list(cols = character(), K = 0))
  matches_list <- stringr::str_match_all(cm_cols, "cm_(\\d+)_(\\d+)")
  idx <- do.call(rbind, matches_list)
  if (is.null(idx) || nrow(idx) == 0) return(list(cols = cm_cols, K = 0))
  i_vals <- suppressWarnings(as.integer(idx[, 2]))
  j_vals <- suppressWarnings(as.integer(idx[, 3]))
  ij <- c(i_vals, j_vals)
  ij <- ij[!is.na(ij) & is.finite(ij)]
  if (length(ij) == 0) return(list(cols = cm_cols, K = 0))
  list(cols = cm_cols, K = max(ij))
}

row_to_cm <- function(row, K, cm_cols) {
  if (K <= 0) return(NULL)
  m <- matrix(0, nrow = K, ncol = K)
  for (nm in cm_cols) {
    m_idx <- str_match(nm, "cm_(\\d+)_(\\d+)")
    i <- as.integer(m_idx[2])
    j <- as.integer(m_idx[3])
    if (is.na(i) || is.na(j) || i < 1 || j < 1 || i > K || j > K) next
    val <- suppressWarnings(as.numeric(row[[nm]]))
    if (is.na(val)) val <- 0
    m[i, j] <- val
  }
  storage.mode(m) <- "double"
  m
}

# Multi-class MCC (Gorodkin 2004)
mcc_from_cm <- function(C, den0_as_zero = FALSE) {
  if (is.null(C)) return(NA_real_)
  if (!is.matrix(C)) {
    C2 <- tryCatch(as.matrix(C), error = function(...) NULL)
    if (is.null(C2) || !is.matrix(C2)) return(NA_real_)
    C <- C2
  }
  if (all(dim(C) == c(1, 1))) return(NA_real_)
  if (!is.numeric(C)) storage.mode(C) <- "double"
  if (!is.finite(sum(C))) return(NA_real_)
  s <- sum(C)
  if (s == 0) return(NA_real_)
  p <- colSums(C)
  t <- rowSums(C)
  csum <- sum(diag(C))
  sum_pt <- sum(p * t)
  num <- (csum * s) - sum_pt
  den <- sqrt((s^2 - sum(p^2)) * (s^2 - sum(t^2)))
  if (den == 0) {
    if (nrow(C) == 2 && ncol(C) == 2) {
      tp <- C[1, 1]; fn <- C[1, 2]; fp <- C[2, 1]; tn <- C[2, 2]
      num2 <- (tp * tn) - (fp * fn)
      den2 <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      if (den2 > 0) return(as.numeric(num2 / den2))
    }
    return(if (isTRUE(den0_as_zero)) 0 else NA_real_)
  }
  as.numeric(num / den)
}

accuracy_from_cm <- function(C) {
  if (is.null(C) || length(C) == 0) return(NA_real_)
  if (!is.matrix(C)) {
    C2 <- tryCatch(as.matrix(C), error = function(...) NULL)
    if (is.null(C2) || !is.matrix(C2)) return(NA_real_)
    C <- C2
  }
  s <- sum(C)
  if (s == 0) return(NA_real_)
  sum(diag(C)) / s
}

f1_macro_from_cm <- function(C) {
  if (is.null(C) || length(C) == 0) return(NA_real_)
  if (!is.matrix(C)) {
    C2 <- tryCatch(as.matrix(C), error = function(...) NULL)
    if (is.null(C2) || !is.matrix(C2)) return(NA_real_)
    C <- C2
  }
  K <- suppressWarnings(as.integer(nrow(C)))
  if (is.na(K) || !is.finite(K) || K <= 0) return(NA_real_)
  f1s <- rep(NA_real_, K)
  p_k <- colSums(C)
  t_k <- rowSums(C)
  for (k in seq_len(K)) {
    tp <- C[k, k]
    prec <- if (p_k[k] > 0) tp / p_k[k] else NA_real_
    rec  <- if (t_k[k] > 0) tp / t_k[k] else NA_real_
    f1s[k] <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec)
  }
  mean(f1s, na.rm = TRUE)
}

f1_micro_from_cm <- function(C) {
  if (is.null(C) || length(C) == 0) return(NA_real_)
  if (!is.matrix(C)) {
    C2 <- tryCatch(as.matrix(C), error = function(...) NULL)
    if (is.null(C2) || !is.matrix(C2)) return(NA_real_)
    C <- C2
  }
  tp <- sum(diag(C))
  fp <- sum(colSums(C)) - tp
  fn <- sum(rowSums(C)) - tp
  denom <- (2 * tp + fp + fn)
  if (denom == 0) return(NA_real_)
  2 * tp / denom
}

bal_acc_from_cm <- function(C) {
  if (is.null(C) || length(C) == 0) return(NA_real_)
  if (!is.matrix(C)) {
    C2 <- tryCatch(as.matrix(C), error = function(...) NULL)
    if (is.null(C2) || !is.matrix(C2)) return(NA_real_)
    C <- C2
  }
  t_k <- rowSums(C)
  K <- nrow(C)
  recalls <- vapply(seq_len(K), function(k) if (t_k[k] > 0) C[k, k] / t_k[k] else NA_real_, numeric(1))
  mean(recalls, na.rm = TRUE)
}

compute_metrics <- function(df, mcc_strategy = c("gorodkin_zero", "gorodkin_strict", "ovr_macro")) {
  mcc_strategy <- match.arg(mcc_strategy)
  spec <- detect_cm_spec(df)
  K <- spec$K
  cm_cols <- spec$cols
  if (K == 0) {
    return(df %>% mutate(MCC = NA_real_, Accuracy = NA_real_, F1_macro = NA_real_, F1_micro = NA_real_, BalancedAccuracy = NA_real_))
  }

  suppressWarnings({
    for (nm in cm_cols) {
      if (nm %in% names(df)) {
        df[[nm]] <- as.numeric(df[[nm]])
        df[[nm]][is.na(df[[nm]])] <- 0
      }
    }
  })

  n <- nrow(df)
  cm_list <- vector("list", n)
  for (i in seq_len(n)) {
    C <- row_to_cm(df[i, , drop = FALSE], K, cm_cols)
    if (is.matrix(C)) {
      dimnames(C) <- NULL
      storage.mode(C) <- "double"
    }
    cm_list[[i]] <- C
  }

  .bin_mcc <- function(tp, fp, fn, tn, zero_on_den0 = TRUE) {
    num2 <- (tp * tn) - (fp * fn)
    den2 <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if (den2 == 0) return(if (zero_on_den0) 0 else NA_real_)
    as.numeric(num2 / den2)
  }

  .mcc_ovr_macro <- function(C, zero_on_den0 = TRUE) {
    if (is.null(C) || !is.matrix(C)) return(NA_real_)
    K <- nrow(C)
    if (K < 2) return(NA_real_)
    s <- sum(C)
    if (s == 0) return(NA_real_)
    vals <- numeric(K)
    for (k in seq_len(K)) {
      tp <- C[k, k]
      fp <- sum(C[-k, k, drop = FALSE])
      fn <- sum(C[k, -k, drop = FALSE])
      tn <- s - tp - fp - fn
      vals[k] <- .bin_mcc(tp, fp, fn, tn, zero_on_den0 = zero_on_den0)
    }
    mean(vals, na.rm = TRUE)
  }

  safe_apply <- function(f) {
    vapply(cm_list, function(C) {
      val <- tryCatch(f(C), error = function(e) NA_real_)
      if (is.na(val) || !is.finite(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
  }

  mcc_vec <- switch(
    mcc_strategy,
    gorodkin_zero = vapply(cm_list, function(C) tryCatch(mcc_from_cm(C, den0_as_zero = TRUE), error = function(e) NA_real_), numeric(1)),
    gorodkin_strict = vapply(cm_list, function(C) tryCatch(mcc_from_cm(C, den0_as_zero = FALSE), error = function(e) NA_real_), numeric(1)),
    ovr_macro = vapply(cm_list, function(C) tryCatch(.mcc_ovr_macro(C, zero_on_den0 = TRUE), error = function(e) NA_real_), numeric(1))
  )

  if (mcc_strategy == "gorodkin_zero") {
    mcc_strict <- vapply(cm_list, function(C) tryCatch(mcc_from_cm(C, den0_as_zero = FALSE), error = function(e) NA_real_), numeric(1))
    df$.mcc_converted_count <- sum(is.na(mcc_strict) & !is.na(mcc_vec) & mcc_vec == 0)
  }

  df$MCC <- mcc_vec
  df$Accuracy <- safe_apply(accuracy_from_cm)
  df$F1_macro <- safe_apply(f1_macro_from_cm)
  df$F1_micro <- safe_apply(f1_micro_from_cm)
  df$BalancedAccuracy <- safe_apply(bal_acc_from_cm)
  df
}

# --- Plot helpers -------------------------------------------------------------

.pretty_labels <- list(
  MCC = "Matthews correlation coefficient (MCC)",
  diag_max_depth = "Maximum depth",
  diag_leaf_nodes = "Leaf count",
  diag_min_depth = "Minimum depth",
  diag_avg_depth = "Average depth",
  diag_avg_weighted_depth = "Average weighted depth",
  diag_decision_nodes = "Decision-node count",
  diag_unique_vars = "Unique-variable count"
)

pretty_label <- function(var) {
  if (!is.null(.pretty_labels[[var]])) .pretty_labels[[var]] else var
}

# Detect a runtime column (optional). Prefers an exact 'runtime' match.
detect_runtime_col <- function(df) {
  nms <- names(df)
  if ("runtime" %in% nms) return("runtime")
  # Common alternatives: elapsed, time_sec, fit_time, predict_time, runtime_sec, seconds, duration
  cand <- grep("(runtime|elapsed|duration|time(_sec|_s)?|seconds|sec)$", nms, ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, grep("^cm_\\d+_\\d+$", nms, value = TRUE))
  if (length(cand) == 0) return(NULL)
  cand[[1]]
}


scatter_mcc <- function(df, x_var, dataset_name, subtitle, out_fp, alpha = 0.45, size = 1.6) {
  df <- df %>% filter(!is.na(.data[[x_var]]), !is.na(.data$MCC), !is.na(.data$algo_crit))
  if (nrow(df) == 0) return(invisible(FALSE))

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data$MCC, color = .data$algo_crit)) +
    geom_point(alpha = alpha, size = size) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(x_var),
      y = pretty_label("MCC"),
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}

scatter_mcc_repeat_avg <- function(df_repeat, x_mean_col, dataset_name, subtitle, out_fp) {
  df_repeat <- df_repeat %>% filter(!is.na(.data[[x_mean_col]]), !is.na(.data$MCC_mean), !is.na(.data$algo_crit))
  if (nrow(df_repeat) == 0) return(invisible(FALSE))

  p <- ggplot(df_repeat, aes(x = .data[[x_mean_col]], y = .data$MCC_mean, color = .data$algo_crit)) +
    geom_point(alpha = 0.65, size = 2.0) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(str_replace(x_mean_col, "_mean$", "")),
      y = "Mean MCC (folds averaged within repeat)",
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}

scatter_mcc_groupmeans_by_x <- function(df, x_var, dataset_name, subtitle, out_fp) {
  df <- df %>% filter(!is.na(.data[[x_var]]), !is.na(.data$MCC), !is.na(.data$algo_crit))
  if (nrow(df) == 0) return(invisible(FALSE))

  means <- df %>%
    group_by(algo_crit, .data[[x_var]]) %>%
    summarize(mean_mcc = mean(MCC, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

  if (nrow(means) == 0) return(invisible(FALSE))

  p <- ggplot(means, aes(x = .data[[x_var]], y = .data$mean_mcc, color = .data$algo_crit)) +
    geom_point(size = 2.3) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(x_var),
      y = "Mean MCC (per x value)",
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}


# Variants of scatter plots with trend lines (geom_smooth), keeping originals unchanged.
scatter_mcc_smooth <- function(df, x_var, dataset_name, subtitle, out_fp,
                               alpha = 0.45, size = 1.6, method = "lm") {
  df <- df %>% filter(!is.na(.data[[x_var]]), !is.na(.data$MCC), !is.na(.data$algo_crit))
  if (nrow(df) == 0) return(invisible(FALSE))

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data$MCC, color = .data$algo_crit)) +
    geom_point(alpha = alpha, size = size) +
    geom_smooth(se = FALSE, method = method) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(x_var),
      y = pretty_label("MCC"),
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}

scatter_mcc_repeat_avg_smooth <- function(df_repeat, x_mean_col, dataset_name, subtitle, out_fp, method = "lm") {
  df_repeat <- df_repeat %>% filter(!is.na(.data[[x_mean_col]]), !is.na(.data$MCC_mean), !is.na(.data$algo_crit))
  if (nrow(df_repeat) == 0) return(invisible(FALSE))

  p <- ggplot(df_repeat, aes(x = .data[[x_mean_col]], y = .data$MCC_mean, color = .data$algo_crit)) +
    geom_point(alpha = 0.65, size = 2.0) +
    geom_smooth(se = FALSE, method = method) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(str_replace(x_mean_col, "_mean$", "")),
      y = "Mean MCC (folds averaged within repeat)",
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}

scatter_mcc_groupmeans_by_x_smooth <- function(df, x_var, dataset_name, subtitle, out_fp, method = "lm") {
  df <- df %>% filter(!is.na(.data[[x_var]]), !is.na(.data$MCC), !is.na(.data$algo_crit))
  if (nrow(df) == 0) return(invisible(FALSE))

  means <- df %>%
    group_by(algo_crit, .data[[x_var]]) %>%
    summarize(mean_mcc = mean(MCC, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

  if (nrow(means) == 0) return(invisible(FALSE))

  p <- ggplot(means, aes(x = .data[[x_var]], y = .data$mean_mcc, color = .data$algo_crit)) +
    geom_point(size = 2.3) +
    geom_smooth(se = FALSE, method = method) +
    labs(
      title = paste0(dataset_name, " — ", subtitle),
      x = pretty_label(x_var),
      y = "Mean MCC (per x value)",
      color = "Algorithm / split criterion"
    ) +
    theme_minimal()

  ggsave(filename = out_fp, plot = p, width = 8, height = 5.5, dpi = 150)
  TRUE
}


build_algo_matrix <- function(df, target_var) {
  if (!target_var %in% names(df)) return(NULL)
  if (!"algo_crit" %in% names(df)) return(NULL)

  df2 <- df %>% filter(!is.na(.data[[target_var]]), !is.na(algo_crit))
  if (nrow(df2) == 0) return(NULL)

  # Respect factor order if already set; otherwise sort labels.
  if (!is.factor(df2$algo_crit)) {
    df2$algo_crit <- factor(as.character(df2$algo_crit), levels = sort(unique(as.character(df2$algo_crit)), method = "radix"))
  }

  split_list <- split(df2[[target_var]], df2$algo_crit)
  m_len <- max(vapply(split_list, length, integer(1)))
  m <- do.call(cbind, lapply(split_list, function(v) { length(v) <- m_len; v }))
  colnames(m) <- names(split_list)
  m
}

try_md_plot_matrix <- function(mat, title, out_fp, verbose = TRUE) {
  if (!requireNamespace("DataVisualizations", quietly = TRUE)) {
    log_warn("Package 'DataVisualizations' not available. Skipping MD plot: ", title)
    return(invisible(FALSE))
  }
  M <- tryCatch(as.matrix(mat), error = function(e) NULL)
  if (is.null(M)) return(invisible(FALSE))

  tmp_png <- tempfile(fileext = ".png")
  grDevices::png(filename = tmp_png, width = 1200, height = 800, res = 96)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  # Ordering = 'Columnwise' to respect the matrix column order
  res <- try(DataVisualizations::MDplot(M, Ordering = 'Columnwise', OnlyPlotOutput = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) {
    log_warn("DataVisualizations::MDplot failed for ", title, ": ", as.character(res))
    return(invisible(FALSE))
  }

  p <- if (is.list(res) && !is.null(res$ggplotObj)) res$ggplotObj else res
  if (!inherits(p, "ggplot")) return(invisible(FALSE))

  p <- p + ggplot2::ggtitle(title)
  ggplot2::ggsave(filename = out_fp, plot = p, width = 12, height = 8, dpi = 150)
  log_info("Saved MD plot to ", out_fp, verbose = verbose)
  invisible(TRUE)
}

# --- Main workflow ------------------------------------------------------------

run_benchmark_metrics <- function(
  benchmark_dir = file.path(getwd(), "tools", "benchmark"),
  metrics_out_dir = file.path(getwd(), "tools", "benchmark_metrics"),
  graphics_out_dir = file.path(getwd(), "tools", "benchmark_graphics"),
  verbose = TRUE,
  continue_on_error = TRUE
) {
  log_info("Starting benchmark metrics run", verbose = verbose)
  log_info("  Working dir: ", getwd(), verbose = verbose)
  log_info("  R version: ", R.version.string, verbose = verbose)

  if (!dir.exists(benchmark_dir)) stop(paste("Benchmark dir not found:", benchmark_dir))
  if (!dir.exists(metrics_out_dir)) dir.create(metrics_out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(graphics_out_dir)) dir.create(graphics_out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(benchmark_dir, pattern = "_cv_benchmark\\.csv$", full.names = TRUE)
  log_info("Found ", length(files), " benchmark file(s).", verbose = verbose)
  if (length(files) == 0) stop("No benchmark CSV files found.")

  t_start <- Sys.time()
  n_ok <- 0L; n_skip <- 0L; n_err <- 0L

  # For optional global plots later (kept as before, but now secondary-level only)
  global_mcc_list <- list()
  global_acc_list <- list()

  # Dataset-balanced per-algorithm/split aggregates (each dataset contributes equally)
  dataset_algo_means_list <- list()


  for (fp in files) {
    ds_name_base <- tools::file_path_sans_ext(basename(fp))
    one_ok <- FALSE

    tryCatch({
      df <- suppressMessages(readr::read_csv(fp, show_col_types = FALSE))

# Exclude LDATree from all processing (drop any rows where algorithm contains 'ldatree')
if ("algorithm" %in% names(df)) {
  .n_before <- nrow(df)
  df <- df %>%
    dplyr::mutate(algorithm = as.character(algorithm)) %>%
    dplyr::filter(!stringr::str_detect(tolower(.data$algorithm), "ldatree"))
  .n_after <- nrow(df)
  if (isTRUE(verbose) && (.n_before != .n_after)) {
    log_info(ds_name_base, ": excluded ", (.n_before - .n_after),
             " row(s) with algorithm matching 'ldatree'.", verbose = TRUE)
  }
}

      df_pred <- df %>% filter(tolower(mode) == "predict")
      df_fit <- df %>% filter(tolower(mode) == "fit")
      if (nrow(df_pred) == 0) {
        log_warn("No predict rows found in ", basename(fp), "; skipping.")
        n_skip <- n_skip + 1L
        next
      }

      ds_name <- if ("dataset" %in% names(df_pred)) as.character(df_pred$dataset[[1]]) else ds_name_base
      if (is.na(ds_name) || ds_name == "") ds_name <- ds_name_base

      # Compute metrics
      df_pred_m <- compute_metrics(df_pred)

      if (".mcc_converted_count" %in% names(df_pred_m)) {
        log_info(ds_name, ": converted ", df_pred_m$.mcc_converted_count[1], " undefined MCC values to 0.", verbose = verbose)
        df_pred_m$.mcc_converted_count <- NULL
      }

      # Add algo_crit and enforce alphabetical ordering across plots
      df_pred_m <- df_pred_m %>% add_algo_crit() %>% apply_algo_crit_levels()

      # Collect dataset-balanced means per (algorithm, criterion) for overall summaries later
      # Collect dataset-balanced means per (algorithm, criterion) for overall summaries later
      # Runtime is taken from FIT rows (mode == "fit"), not from predict rows.
      runtime_col_fit <- if (nrow(df_fit) > 0) detect_runtime_col(df_fit) else NULL

      # IMPORTANT: normalize criterion in FIT rows the same way as in PREDICT rows
      # so joins on (algorithm, criterion) work even when criterion is NA/empty.
      df_fit_norm <- if (nrow(df_fit) > 0) {
        df_fit %>%
          mutate(
            algorithm = as.character(algorithm),
            criterion = ifelse(is.na(criterion) | criterion == "", "(none)", as.character(criterion))
          )
      } else {
        df_fit
      }

      runtime_fit_tbl <- if (!is.null(runtime_col_fit) && runtime_col_fit %in% names(df_fit_norm) && nrow(df_fit_norm) > 0) {
        df_fit_norm %>%
          group_by(algorithm, criterion) %>%
          summarize(runtime = mean(.data[[runtime_col_fit]], na.rm = TRUE), .groups = "drop")
      } else {
        if (is.null(runtime_col_fit)) {
          log_warn(ds_name_base, ": no runtime column detected in FIT rows (mode == 'fit'). Runtime summaries/plots may be incomplete.")
        } else {
          log_warn(ds_name_base, ": runtime column '", runtime_col_fit, "' not found/usable in FIT rows. Runtime summaries/plots may be incomplete.")
        }
        dplyr::tibble(algorithm = character(0), criterion = character(0), runtime = numeric(0))
      }

      ds_means <- df_pred_m %>%
        group_by(algorithm, criterion) %>%
        summarize(
          dataset = ds_name,
          MCC = if ("MCC" %in% names(.)) mean(MCC, na.rm = TRUE) else NA_real_,
          diag_max_depth = if ("diag_max_depth" %in% names(.)) mean(diag_max_depth, na.rm = TRUE) else NA_real_,
          diag_leaf_nodes = if ("diag_leaf_nodes" %in% names(.)) mean(diag_leaf_nodes, na.rm = TRUE) else NA_real_,
          .groups = "drop"
        ) %>%
        left_join(runtime_fit_tbl, by = c("algorithm", "criterion")) %>%
        add_algo_crit() %>%
        apply_algo_crit_levels()


      # Diagnostics: if some algo/criterion pairs have runtime in FIT rows but end up NA after join,
      # print the affected combinations and a quick check of what keys exist.
      if (nrow(runtime_fit_tbl) > 0) {
        missing_rt <- ds_means %>%
          filter(is.na(runtime)) %>%
          distinct(algorithm, criterion, algo_crit)

        if (nrow(missing_rt) > 0) {
          log_warn(ds_name, ": runtime missing after join for ", nrow(missing_rt),
                   " algorithm/split combination(s). This usually means key mismatch between FIT and PREDICT rows.")
          if (isTRUE(verbose)) {
            print(missing_rt)
            log_info(ds_name, ": FIT runtime keys (first 30):", verbose = TRUE)
            print(head(runtime_fit_tbl %>% mutate(algo_crit = paste0(algorithm, " / ", criterion)) %>% arrange(algo_crit), 30))
          }
        }
      }
dataset_algo_means_list[[ds_name]] <- ds_means


      # Save augmented CSV
      out_csv <- file.path(metrics_out_dir, paste0(ds_name, "_cv_benchmark_metrics.csv"))
      suppressMessages(readr::write_csv(df_pred_m, out_csv))

      # Aggregations
      metrics_cols <- c("MCC", "Accuracy", "F1_macro", "F1_micro", "BalancedAccuracy")
      diag_cols <- grep("^diag_", names(df_pred_m), value = TRUE)
      all_agg_cols <- intersect(c(metrics_cols, diag_cols), names(df_pred_m))

      # Mean per (algorithm, criterion, fold) across repeats
      table_fold <- df_pred_m %>%
        group_by(algorithm, criterion, fold) %>%
        summarize(across(all_of(all_agg_cols), list(mean = ~mean(.x, na.rm = TRUE))), .groups = "drop") %>%
        add_algo_crit() %>%
        apply_algo_crit_levels()
      readr::write_csv(table_fold, file.path(metrics_out_dir, paste0(ds_name, "_metrics_by_fold.csv")))

      # Mean per (algorithm, criterion, repeat) across folds  <-- required for "folds averaged within repeat"
      if ("repeat" %in% names(df_pred_m)) {
        table_repeat <- df_pred_m %>%
          group_by(algorithm, criterion, `repeat`) %>%
          summarize(across(all_of(all_agg_cols), list(mean = ~mean(.x, na.rm = TRUE))), .groups = "drop") %>%
          add_algo_crit() %>%
          apply_algo_crit_levels()
        readr::write_csv(table_repeat, file.path(metrics_out_dir, paste0(ds_name, "_metrics_by_repeat.csv")))
      } else {
        table_repeat <- NULL
        log_warn(ds_name, ": column 'repeat' not found. Repeat-averaged plots will be skipped.")
      }

      # Scientific summary per algorithm/split pair
      summary_vars <- intersect(c("MCC", "diag_max_depth", "diag_leaf_nodes"), names(df_pred_m))
      if (length(summary_vars) > 0) {
        scientific_summary <- df_pred_m %>%
          group_by(algorithm, criterion) %>%
          summarize(across(all_of(summary_vars), list(
            mean   = ~mean(.x, na.rm = TRUE),
            sd     = ~sd(.x, na.rm = TRUE),
            min    = ~min(.x, na.rm = TRUE),
            q25    = ~quantile(.x, 0.25, na.rm = TRUE),
            median = ~median(.x, na.rm = TRUE),
            q75    = ~quantile(.x, 0.75, na.rm = TRUE),
            max    = ~max(.x, na.rm = TRUE),
            var    = ~var(.x, na.rm = TRUE)
          )), .groups = "drop") %>%
          add_algo_crit() %>%
          apply_algo_crit_levels()
        readr::write_csv(scientific_summary, file.path(metrics_out_dir, paste0(ds_name, "_scientific_summary.csv")))
      }

      # Graphics dirs
      ds_graphics_dir <- file.path(graphics_out_dir, ds_name)
      secondary_dir <- file.path(ds_graphics_dir, "secondary_plots")
      if (!dir.exists(ds_graphics_dir)) dir.create(ds_graphics_dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(secondary_dir)) dir.create(secondary_dir, recursive = TRUE, showWarnings = FALSE)

      # --- PRIMARY PLOTS (requested) -----------------------------------------

      # Raw scatter (all datapoints)
      if (all(c("MCC", "diag_max_depth", "algo_crit") %in% names(df_pred_m))) {
        scatter_mcc(
          df_pred_m, "diag_max_depth", ds_name,
          "MCC vs maximum depth (all folds, all repeats)",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_raw.png"))
        )

        # Smooth trend variant (no CI)
        scatter_mcc_smooth(
          df_pred_m, "diag_max_depth", ds_name,
          "MCC vs maximum depth (all folds, all repeats) — trend",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_raw_smooth.png"))
        )
      }
      if (all(c("MCC", "diag_leaf_nodes", "algo_crit") %in% names(df_pred_m))) {
        scatter_mcc(
          df_pred_m, "diag_leaf_nodes", ds_name,
          "MCC vs leaf count (all folds, all repeats)",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_raw.png"))
        )

        # Smooth trend variant (no CI)
        scatter_mcc_smooth(
          df_pred_m, "diag_leaf_nodes", ds_name,
          "MCC vs leaf count (all folds, all repeats) — trend",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_raw_smooth.png"))
        )
      }

      # Scatter with folds averaged within each repeat (one point per algo/crit/repeat)
      if (!is.null(table_repeat)) {
        if (all(c("MCC_mean", "diag_max_depth_mean", "algo_crit") %in% names(table_repeat))) {
          scatter_mcc_repeat_avg(
            table_repeat, "diag_max_depth_mean", ds_name,
            "MCC vs maximum depth (folds averaged within repeat)",
            file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_repeat_avg.png"))
          )

          # Smooth trend variant (no CI)
          scatter_mcc_repeat_avg_smooth(
            table_repeat, "diag_max_depth_mean", ds_name,
            "MCC vs maximum depth (folds averaged within repeat) — trend",
            file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_repeat_avg_smooth.png"))
          )
        }
        if (all(c("MCC_mean", "diag_leaf_nodes_mean", "algo_crit") %in% names(table_repeat))) {
          scatter_mcc_repeat_avg(
            table_repeat, "diag_leaf_nodes_mean", ds_name,
            "MCC vs leaf count (folds averaged within repeat)",
            file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_repeat_avg.png"))
          )

          # Smooth trend variant (no CI)
          scatter_mcc_repeat_avg_smooth(
            table_repeat, "diag_leaf_nodes_mean", ds_name,
            "MCC vs leaf count (folds averaged within repeat) — trend",
            file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_repeat_avg_smooth.png"))
          )
        }
      }

      # Scatter with mean MCC per x value within algo/crit
      if (all(c("MCC", "diag_max_depth", "algo_crit") %in% names(df_pred_m))) {
        scatter_mcc_groupmeans_by_x(
          df_pred_m, "diag_max_depth", ds_name,
          "MCC vs maximum depth (mean MCC per x value)",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_x_means.png"))
        )

        # Smooth trend variant (no CI)
        scatter_mcc_groupmeans_by_x_smooth(
          df_pred_m, "diag_max_depth", ds_name,
          "MCC vs maximum depth (mean MCC per x value) — trend",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_x_means_smooth.png"))
        )
      }
      if (all(c("MCC", "diag_leaf_nodes", "algo_crit") %in% names(df_pred_m))) {
        scatter_mcc_groupmeans_by_x(
          df_pred_m, "diag_leaf_nodes", ds_name,
          "MCC vs leaf count (mean MCC per x value)",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_x_means.png"))
        )

        # Smooth trend variant (no CI)
        scatter_mcc_groupmeans_by_x_smooth(
          df_pred_m, "diag_leaf_nodes", ds_name,
          "MCC vs leaf count (mean MCC per x value) — trend",
          file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_x_means_smooth.png"))
        )
      }

      # NEW: mean-only scatter with runtime-sized points (per dataset)
      # Uses ds_means: one row per (algorithm, criterion) with dataset-balanced means + FIT runtime.
      if (exists("ds_means") && nrow(ds_means) > 0 &&
          all(c("MCC", "diag_leaf_nodes", "diag_max_depth", "runtime", "algo_crit") %in% names(ds_means))) {

        df_sc <- ds_means %>%
          dplyr::filter(
            is.finite(.data$MCC),
            is.finite(.data$diag_leaf_nodes),
            is.finite(.data$diag_max_depth),
            is.finite(.data$runtime),
            .data$runtime > 0
          )

        if (nrow(df_sc) > 0) {

          p_leaf <- ggplot(df_sc, aes(x = .data$diag_leaf_nodes, y = .data$MCC, color = .data$algo_crit, size = .data$runtime)) +
            geom_point(alpha = 0.85) +
            scale_size_continuous(trans = "log10") +
            labs(
              title = paste0(ds_name, " — MCC vs leaf count (means; size = runtime)"),
              x = "Leaf count (mean)",
              y = "Matthews correlation coefficient (MCC, mean)",
              color = "Algorithm / split criterion",
              size = "Runtime (seconds, log10 scale)"
            ) +
            theme_minimal()

          ggsave(
            filename = file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_size_runtime_means.png")),
            plot = p_leaf, width = 9, height = 6, dpi = 150
          )

# Smooth trend variant (overall trend; no CI)
p_leaf_smooth <- p_leaf +
  geom_smooth(aes(x = diag_leaf_nodes, y = MCC), inherit.aes = FALSE, se = FALSE, method = "lm")

ggsave(
  filename = file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_leaf_count_size_runtime_means_smooth.png")),
  plot = p_leaf_smooth, width = 9, height = 6, dpi = 150
)


          p_depth <- ggplot(df_sc, aes(x = .data$diag_max_depth, y = .data$MCC, color = .data$algo_crit, size = .data$runtime)) +
            geom_point(alpha = 0.85) +
            scale_size_continuous(trans = "log10") +
            labs(
              title = paste0(ds_name, " — MCC vs maximum depth (means; size = runtime)"),
              x = "Maximum depth (mean)",
              y = "Matthews correlation coefficient (MCC, mean)",
              color = "Algorithm / split criterion",
              size = "Runtime (seconds, log10 scale)"
            ) +
            theme_minimal()

          ggsave(
            filename = file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_size_runtime_means.png")),
            plot = p_depth, width = 9, height = 6, dpi = 150
          )

# Smooth trend variant (overall trend; no CI)
p_depth_smooth <- p_depth +
  geom_smooth(aes(x = diag_max_depth, y = MCC), inherit.aes = FALSE, se = FALSE, method = "lm")

ggsave(
  filename = file.path(ds_graphics_dir, paste0(ds_name, "_scatter_mcc_vs_max_depth_size_runtime_means_smooth.png")),
  plot = p_depth_smooth, width = 9, height = 6, dpi = 150
)

        } else if (isTRUE(verbose)) {
          log_info(ds_name, ": runtime-sized mean scatter skipped (no finite runtime>0 / MCC data).", verbose = TRUE)
        }
      }

      # MD plots (MCC) — raw and repeat-averaged
      mcc_algo_mat <- build_algo_matrix(df_pred_m, "MCC")
      if (!is.null(mcc_algo_mat)) {
        try_md_plot_matrix(
          mcc_algo_mat,
          paste0(ds_name, " — MD plot: MCC (all folds, all repeats)"),
          file.path(ds_graphics_dir, paste0(ds_name, "_MD_MCC_raw.png")),
          verbose = FALSE
        )
      }
      if (!is.null(table_repeat)) {
        mcc_repeat_mat <- build_algo_matrix(table_repeat, "MCC_mean")
        if (!is.null(mcc_repeat_mat)) {
          try_md_plot_matrix(
            mcc_repeat_mat,
            paste0(ds_name, " — MD plot: MCC (folds averaged within repeat)"),
            file.path(ds_graphics_dir, paste0(ds_name, "_MD_MCC_repeat_avg.png")),
            verbose = FALSE
          )
        }
      }

      # --- SECONDARY PLOTS (everything else) ---------------------------------

      # Secondary scatter plots for other diagnostics vs MCC (raw points)
      scatter_vars_secondary <- c(
        diag_min_depth = "Minimum depth",
        diag_avg_depth = "Average depth",
        diag_avg_weighted_depth = "Average weighted depth",
        diag_decision_nodes = "Decision-node count",
        diag_unique_vars = "Unique-variable count"
      )

      for (v_name in names(scatter_vars_secondary)) {
        if (!all(c(v_name, "MCC", "algo_crit") %in% names(df_pred_m))) next
        subtitle <- paste0("MCC vs ", tolower(scatter_vars_secondary[[v_name]]), " (all folds, all repeats)")
        out_fp <- file.path(secondary_dir, paste0(ds_name, "_scatter_mcc_vs_", str_replace(v_name, "^diag_", ""), "_raw.png"))
        scatter_mcc(df_pred_m, v_name, ds_name, subtitle, out_fp, alpha = 0.40, size = 1.5)
      }

      # Secondary MD plots (Accuracy + basic tree stats), if present
      acc_algo_mat <- build_algo_matrix(df_pred_m, "Accuracy")
      if (!is.null(acc_algo_mat)) {
        try_md_plot_matrix(
          acc_algo_mat,
          paste0(ds_name, " — MD plot: Accuracy (all folds, all repeats)"),
          file.path(secondary_dir, paste0(ds_name, "_MD_Accuracy_raw.png")),
          verbose = FALSE
        )
      }

      depth_algo_mat <- build_algo_matrix(df_pred_m, "diag_max_depth")
      if (!is.null(depth_algo_mat)) {
        try_md_plot_matrix(
          depth_algo_mat,
          paste0(ds_name, " — MD plot: maximum depth (all folds, all repeats)"),
          file.path(secondary_dir, paste0(ds_name, "_MD_max_depth_raw.png")),
          verbose = FALSE
        )
      }

      leaf_algo_mat <- build_algo_matrix(df_pred_m, "diag_leaf_nodes")
      if (!is.null(leaf_algo_mat)) {
        try_md_plot_matrix(
          leaf_algo_mat,
          paste0(ds_name, " — MD plot: leaf count (all folds, all repeats)"),
          file.path(secondary_dir, paste0(ds_name, "_MD_leaf_count_raw.png")),
          verbose = FALSE
        )
      }

      # Keep globals for optional global plots
      global_mcc_list[[ds_name]] <- df_pred_m$MCC
      if ("Accuracy" %in% names(df_pred_m)) global_acc_list[[ds_name]] <- df_pred_m$Accuracy

      one_ok <- TRUE
    }, error = function(e) {
      n_err <<- n_err + 1L
      log_error("Error while processing ", basename(fp), ": ", conditionMessage(e))
      if (!continue_on_error) stop(e)
    })

    if (one_ok) n_ok <- n_ok + 1L
  }


  # --- OVERALL (across datasets) summaries + boxplots ------------------------

  if (length(dataset_algo_means_list) > 0) {
    overall_df <- dplyr::bind_rows(dataset_algo_means_list) %>%
      add_algo_crit() %>%
      apply_algo_crit_levels()

    overall_metrics_dir <- file.path(metrics_out_dir)
    overall_graphics_dir <- file.path(graphics_out_dir, "_overall_summary")
    if (!dir.exists(overall_graphics_dir)) dir.create(overall_graphics_dir, recursive = TRUE, showWarnings = FALSE)

    .summarize_one_metric <- function(df, metric_col) {
      if (!metric_col %in% names(df)) return(NULL)
      df %>%
        group_by(algorithm, criterion) %>%
        summarize(
          n_datasets = sum(!is.na(.data[[metric_col]])),
          mean   = mean(.data[[metric_col]], na.rm = TRUE),
          median = median(.data[[metric_col]], na.rm = TRUE),
          sd     = sd(.data[[metric_col]], na.rm = TRUE),
          min    = min(.data[[metric_col]], na.rm = TRUE),
          q25    = as.numeric(quantile(.data[[metric_col]], 0.25, na.rm = TRUE)),
          q75    = as.numeric(quantile(.data[[metric_col]], 0.75, na.rm = TRUE)),
          max    = max(.data[[metric_col]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        add_algo_crit() %>%
        apply_algo_crit_levels()
    }

    # Combined table: key metrics side-by-side (dataset-balanced, by algo/criterion)
    overall_combined <- overall_df %>%
      group_by(algorithm, criterion) %>%
      summarize(
        n_datasets_MCC = sum(!is.na(MCC)),
        MCC_mean = mean(MCC, na.rm = TRUE),
        MCC_median = median(MCC, na.rm = TRUE),
        MCC_sd = sd(MCC, na.rm = TRUE),
        MCC_q25 = as.numeric(quantile(MCC, 0.25, na.rm = TRUE)),
        MCC_q75 = as.numeric(quantile(MCC, 0.75, na.rm = TRUE)),

        n_datasets_runtime = sum(!is.na(runtime)),
        runtime_mean = mean(runtime, na.rm = TRUE),
        runtime_median = median(runtime, na.rm = TRUE),
        runtime_sd = sd(runtime, na.rm = TRUE),
        runtime_q25 = as.numeric(quantile(runtime, 0.25, na.rm = TRUE)),
        runtime_q75 = as.numeric(quantile(runtime, 0.75, na.rm = TRUE)),

        n_datasets_max_depth = sum(!is.na(diag_max_depth)),
        max_depth_mean = mean(diag_max_depth, na.rm = TRUE),
        max_depth_median = median(diag_max_depth, na.rm = TRUE),
        max_depth_sd = sd(diag_max_depth, na.rm = TRUE),
        max_depth_q25 = as.numeric(quantile(diag_max_depth, 0.25, na.rm = TRUE)),
        max_depth_q75 = as.numeric(quantile(diag_max_depth, 0.75, na.rm = TRUE)),

        n_datasets_leaf_count = sum(!is.na(diag_leaf_nodes)),
        leaf_count_mean = mean(diag_leaf_nodes, na.rm = TRUE),
        leaf_count_median = median(diag_leaf_nodes, na.rm = TRUE),
        leaf_count_sd = sd(diag_leaf_nodes, na.rm = TRUE),
        leaf_count_q25 = as.numeric(quantile(diag_leaf_nodes, 0.25, na.rm = TRUE)),
        leaf_count_q75 = as.numeric(quantile(diag_leaf_nodes, 0.75, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      add_algo_crit() %>%
      apply_algo_crit_levels()

    readr::write_csv(overall_combined, file.path(overall_metrics_dir, "overall_dataset_balanced_summary_by_algo_crit.csv"))

    # Per-metric tables
    metrics_to_summarize <- c("MCC", "runtime", "diag_max_depth", "diag_leaf_nodes")
    for (mc in metrics_to_summarize) {
      tbl <- .summarize_one_metric(overall_df, mc)
      if (is.null(tbl)) next
      readr::write_csv(tbl, file.path(overall_metrics_dir, paste0("overall_summary_", mc, ".csv")))
    }

    # Boxplots across datasets (distribution = dataset-level means per algo/criterion)
    .boxplot_metric <- function(df, metric_col, out_fp, title_prefix = "Overall") {
      if (!metric_col %in% names(df)) return(invisible(FALSE))
      d <- df %>% filter(!is.na(.data[[metric_col]]), !is.na(algo_crit))
      if (nrow(d) == 0) return(invisible(FALSE))

      y_lab <- pretty_label(metric_col)
      if (metric_col == "runtime") y_lab <- "Runtime (log10 scale)"
      if (metric_col == "diag_leaf_nodes") y_lab <- "Leaf count (log10 scale)"
      if (metric_col == "diag_max_depth") y_lab <- "Maximum depth"


      log_scale <- FALSE
      if (metric_col == "runtime") {
        # Runtime often spans orders of magnitude; use log10 scale for readability.
        d <- d %>% filter(.data[[metric_col]] > 0)
        log_scale <- TRUE
      }
      if (metric_col == "diag_leaf_nodes") {
        # Leaf count can span orders of magnitude; use log10 scale for readability.
        d <- d %>% filter(.data[[metric_col]] > 0)
        log_scale <- TRUE
      }
      if (nrow(d) == 0) return(invisible(FALSE))

      # Reorder algo_crit by the mean of the metric
      d <- d %>%
        mutate(algo_crit = reorder(algo_crit, .data[[metric_col]], FUN = function(x) mean(x, na.rm = TRUE)))

      p <- ggplot(d, aes(x = .data$algo_crit, y = .data[[metric_col]])) +
        geom_boxplot(outlier.alpha = 0.35) +
        labs(
          title = paste0(title_prefix, " — ", y_lab, " across datasets (dataset-level means)"),
          x = "Algorithm / split criterion",
          y = y_lab
        ) +
        (if (log_scale) scale_y_log10() else NULL) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(filename = out_fp, plot = p, width = 12, height = 6.5, dpi = 150)
      invisible(TRUE)
    }

    .boxplot_metric(overall_df, "MCC", file.path(overall_graphics_dir, "Boxplot_MCC_dataset_means.png"))
    .boxplot_metric(overall_df, "runtime", file.path(overall_graphics_dir, "Boxplot_runtime_dataset_means.png"))
    .boxplot_metric(overall_df, "diag_max_depth", file.path(overall_graphics_dir, "Boxplot_max_depth_dataset_means.png"))
    .boxplot_metric(overall_df, "diag_leaf_nodes", file.path(overall_graphics_dir, "Boxplot_leaf_count_dataset_means.png"))
  }

# --- OVERALL scatter (means only): MCC vs complexity with runtime sizing ----
# Uses dataset-balanced means per (algorithm, criterion): one point per algo/split pair.
if (exists("overall_df") && exists("overall_graphics_dir") && nrow(overall_df) > 0) {

  overall_pair_means <- overall_df %>%
    group_by(algorithm, criterion) %>%
    summarize(
      MCC        = mean(MCC, na.rm = TRUE),
      runtime    = mean(runtime, na.rm = TRUE),
      max_depth  = mean(diag_max_depth, na.rm = TRUE),
      leaf_count = mean(diag_leaf_nodes, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    add_algo_crit() %>%
    apply_algo_crit_levels()

  # runtime used for point sizing; require positive finite values for log sizing
  overall_pair_means <- overall_pair_means %>%
    filter(is.finite(runtime), runtime > 0, is.finite(MCC))

  .scatter_mcc_runtime <- function(df, x_col, x_lab, out_fp, title_suffix) {
    if (!x_col %in% names(df)) return(invisible(FALSE))
    d <- df %>% filter(is.finite(.data[[x_col]]))
    if (nrow(d) == 0) return(invisible(FALSE))

    p <- ggplot(
      d,
      aes(
        x = .data[[x_col]],
        y = .data$MCC,
        color = .data$algo_crit,
        size = .data$runtime
      )
    ) +
      geom_point(alpha = 0.85) +
      geom_text(aes(label = .data$algo_crit), check_overlap = TRUE, show.legend = FALSE, vjust = -0.6) +
      scale_size_continuous(trans = "log10") +
      labs(
        title = paste0("Overall — MCC vs ", title_suffix, " (means only; size = runtime)"),
        x = x_lab,
        y = "Matthews correlation coefficient (MCC)",
        color = "Algorithm / split criterion",
        size = "Runtime (s, log10 size)"
      ) +
      theme_minimal()

    ggsave(filename = out_fp, plot = p, width = 12, height = 7, dpi = 150)
    invisible(TRUE)
  }

  .scatter_mcc_runtime(
    overall_pair_means,
    x_col = "leaf_count",
    x_lab = "Leaf count (mean across datasets)",
    out_fp = file.path(overall_graphics_dir, "Scatter_MCC_vs_leaf_count_size_runtime_means.png"),
    title_suffix = "leaf count"
  )

  .scatter_mcc_runtime(
    overall_pair_means,
    x_col = "max_depth",
    x_lab = "Maximum depth (mean across datasets)",
    out_fp = file.path(overall_graphics_dir, "Scatter_MCC_vs_max_depth_size_runtime_means.png"),
    title_suffix = "maximum depth"
  )
}

}


# --- End-of-file convenience runner -------------------------------------------
# When sourced interactively, run with defaults so it "does something" from console.
if (interactive()) {
  log_info("benchmark_metrics.R sourced interactively — running run_benchmark_metrics() with defaults...", verbose = TRUE)
  tryCatch(
    run_benchmark_metrics(),
    error = function(e) log_error("run_benchmark_metrics() failed: ", conditionMessage(e))
  )
} else {
  message("Loaded benchmark_metrics.R. To run: run_benchmark_metrics()")
}
