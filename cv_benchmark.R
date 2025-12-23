## -----------------------------------------------------------------------------
## LRN/CLS CV Benchmark Harness
## - Scans TBC/data for paired .LRN (numeric matrix) and .CLS (class vector)
## - For each dataset: 10-fold stratified CV, repeated 100 times
## - Runs all engines/criteria like penguins_sdt_harness plus TBC (bivariate = FALSE/TRUE)
## - Measures time and (best-effort) RAM for fit and predict
## - Computes unified tree diagnostics via get_tree_rule_stats()
## - Writes one CSV per dataset with rows for both fit and predict of each fold
## -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

hr <- function(x) {
  cat("\n", paste(rep("-", nchar(x) + 8), collapse = ""), "\n", sep = "")
  cat("--- ", x, " ---\n", sep = "")
}

## Try to make the TBC functions available
load_TBC <- function() {
  if (suppressWarnings(requireNamespace("TBC", quietly = TRUE))) {
    suppressPackageStartupMessages(library(TBC))
    return(invisible(TRUE))
  }
  if (suppressWarnings(requireNamespace("devtools", quietly = TRUE))) {
    ok <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  # Fallback: source R/*.R directly
  rdir <- file.path(getwd(), "R")
  rfiles <- list.files(rdir, pattern = "\\.R$", full.names = TRUE)
  for (f in rfiles) {
    try(suppressWarnings(sys.source(f, envir = .GlobalEnv)), silent = TRUE)
  }
  invisible(TRUE)
}

## Explicitly load packages so missing ones are visible ------------------------
load_required_packages <- function() {
  core_pkgs <- c(
    "TBC",          # in-package API
    "cdbt.DataIO", # ReadLRN / ReadCLS provider
    "partykit",
    "rpart",
    "tree",
    "C50",
    "RWeka",
    "evtree",
    "LDATree",
    "CORElearn"
  )
  for (p in core_pkgs) {
    tryCatch({
      library(p, character.only = TRUE)
      message(sprintf("Loaded package: %s", p))
    }, error = function(e) {
      message(sprintf("Package '%s' not available: %s", p, conditionMessage(e)))
    })
  }

  # Parallel backends (optional)
  for (p in c("future", "future.apply")) {
    tryCatch({
      library(p, character.only = TRUE)
      message(sprintf("Loaded package: %s", p))
    }, error = function(e) {
      message(sprintf("Package '%s' not available (parallel fallback to sequential): %s", p, conditionMessage(e)))
    })
  }
}

## Build placeholder rows on error --------------------------------------------
.error_rows_for_case <- function(ds_name, rep_id, fold_id, case, levels_ref, n_train = NA_integer_, n_test = NA_integer_, err = "") {
  cm_template <- as.list(setNames(rep(0, length(levels_ref)^2), paste0("cm_", rep(levels_ref, each = length(levels_ref)), "_", rep(levels_ref, times = length(levels_ref)))))
  diag_na <- list(
    diag_min_depth = NA, diag_avg_depth = NA, diag_max_depth = NA,
    diag_avg_weighted_depth = NA, diag_decision_nodes = NA, diag_leaf_nodes = NA, diag_unique_vars = NA
  )
  fit_row <- c(
    list(mode = "fit", dataset = ds_name, algorithm = case$algorithm, criterion = case$criterion,
         `repeat` = rep_id, fold = fold_id, n_train = n_train, n_test = n_test,
         time_sec = NA_real_, mem_mb = NA_real_, error = as.character(err)),
    diag_na,
    cm_template
  )
  pred_row <- c(
    list(mode = "predict", dataset = ds_name, algorithm = case$algorithm, criterion = case$criterion,
         `repeat` = rep_id, fold = fold_id, n_train = n_train, n_test = n_test,
         time_sec = NA_real_, mem_mb = NA_real_, error = as.character(err)),
    diag_na,
    cm_template
  )
  list(fit_row = fit_row, pred_row = pred_row)
}

## Ensure tree statistics function is available
ensure_stats_fn <- function() {
  if (exists("get_tree_rule_stats", mode = "function")) return(invisible(TRUE))
  td_path <- file.path(getwd(), "tools", "tree_diagnostics.R")
  if (file.exists(td_path)) {
    try(suppressWarnings(sys.source(td_path, envir = .GlobalEnv)), silent = TRUE)
  }
  invisible(TRUE)
}

## Dataset discovery -----------------------------------------------------------
list_lrn_cls_datasets <- function(data_dir = file.path(getwd(), "data")) {
  lrn <- list.files(data_dir, pattern = "\\.LRN$", full.names = TRUE, ignore.case = TRUE)
  cls <- list.files(data_dir, pattern = "\\.CLS$", full.names = TRUE, ignore.case = TRUE)
  base <- function(x) tools::file_path_sans_ext(basename(x))
  m <- intersect(tolower(base(lrn)), tolower(base(cls)))
  out <- lapply(m, function(b) {
    lrn_fp <- lrn[tolower(base(lrn)) == b][1L]
    cls_fp <- cls[tolower(base(cls)) == b][1L]
    list(name = b, lrn = lrn_fp, cls = cls_fp)
  })
  out
}

read_lrn_cls <- function(lrn_fp, cls_fp, data_dir = dirname(lrn_fp)) {
  if (!exists("ReadLRN", mode = "function") || !exists("ReadCLS", mode = "function")) {
    stop("ReadLRN/ReadCLS not found. Please ensure the appropriate reader functions are available.")
  }
  # Suppress noisy version messages from readers
  lrn_obj <- suppressWarnings(suppressMessages(
    ReadLRN(FileName = basename(lrn_fp), InDirectory = data_dir)
  ))
  cls_obj <- suppressWarnings(suppressMessages(
    ReadCLS(FileName = basename(cls_fp), InDirectory = data_dir)
  ))
  X <- as.data.frame(lrn_obj$Data)
  y <- cls_obj$Cls
  # coerce response to factor
  if (!is.factor(y)) y <- factor(y)
  # drop rows with NA in X or y
  keep <- stats::complete.cases(X, y)
  X <- X[keep, , drop = FALSE]
  y <- droplevels(y[keep])
  # ensure column names
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  list(X = X, y = y)
}

## Stratified CV ---------------------------------------------------------------
make_stratified_folds <- function(y, k = 10L, repeats = 100L, seed = 1L) {
  set.seed(seed)
  y <- droplevels(factor(y))
  n <- length(y)
  cls <- levels(y)
  out <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    # indices per class
    idx_by_class <- lapply(cls, function(cl) which(y == cl))
    # shuffle within class
    idx_by_class <- lapply(idx_by_class, function(v) if (length(v) > 0L) sample(v) else integer(0))
    # assign folds round-robin per class, ensuring exactly k lists
    folds <- vector("list", k)
    for (i in seq_len(k)) folds[[i]] <- integer(0)
    for (c_i in seq_along(idx_by_class)) {
      v <- idx_by_class[[c_i]]
      if (length(v) == 0L) next
      parts <- vector("list", k)
      for (i in seq_along(v)) {
        j <- ((i - 1L) %% k) + 1L
        parts[[j]] <- c(parts[[j]], v[i])
      }
      for (i in seq_len(k)) folds[[i]] <- c(folds[[i]], parts[[i]])
    }
    # shuffle each fold to mix classes
    folds <- lapply(folds, function(v) if (length(v) > 0L) sample(v) else v)
    out[[r]] <- folds
  }
  out
}

## Engine registry -------------------------------------------------------------
engine_grid <- function() {
  # Returns a list of cases with: engine, algorithm, criterion, engine_args, control_overrides
  cases <- list(
    list(engine = "rpart",   algorithm = "rpart",   criterion = "gini",         engine_args = list(parms = list(split = "gini")),          control_overrides = list()),
    list(engine = "rpart",   algorithm = "rpart",   criterion = "information",  engine_args = list(parms = list(split = "information")),   control_overrides = list()),
    list(engine = "ctree",   algorithm = "ctree",   criterion = "",             engine_args = list(),                                        control_overrides = list()),
    list(engine = "tree",    algorithm = "tree",    criterion = "deviance",     engine_args = list(split = "deviance"),                    control_overrides = list()),
    list(engine = "J48",     algorithm = "J48",     criterion = "",             engine_args = list(),                                        control_overrides = list()),
    list(engine = "C50",     algorithm = "C5.0",    criterion = "",             engine_args = list(),                                        control_overrides = list()),
    list(engine = "evtree",  algorithm = "evtree",  criterion = "",             engine_args = list(),                                        control_overrides = list()),
    list(engine = "ldatree", algorithm = "ldatree", criterion = "",             engine_args = list(),                                        control_overrides = list()),
    # CoreLearn variants
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "Gini",               engine_args = list(model = "tree", selectionEstimator = "Gini"),               control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "InfGain",            engine_args = list(model = "tree", selectionEstimator = "InfGain"),            control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "GainRatio",         engine_args = list(model = "tree", selectionEstimator = "GainRatio"),         control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "MDL",               engine_args = list(model = "tree", selectionEstimator = "MDL"),               control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "Accuracy",          engine_args = list(model = "tree", selectionEstimator = "Accuracy"),          control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "ImpurityHellinger", engine_args = list(model = "tree", selectionEstimator = "ImpurityHellinger"), control_overrides = list()),
    list(engine = "CoreLearn", algorithm = "CoreLearn", criterion = "ImpurityEuclid",    engine_args = list(model = "tree", selectionEstimator = "ImpurityEuclid"),    control_overrides = list()),
    # TBC (in-package) variants
    list(engine = "tbc", algorithm = "tbc", criterion = "axis",    engine_args = list(tbc_control = list(bivariate = FALSE)), control_overrides = list()),
    list(engine = "tbc", algorithm = "tbc", criterion = "bivar",   engine_args = list(tbc_control = list(bivariate = TRUE)),  control_overrides = list())
  )
  cases
}

## Measurement helpers --------------------------------------------------------
now_time <- function() proc.time()[["elapsed"]]

mem_used_mb <- function() {
  # Best-effort platform-dependent memory snapshot; may return NA
  mb <- NA_real_
  if ("pryr" %in% rownames(installed.packages())) {
    # Not mem_change (needs expr); we approximate by mem_used()
    m <- try(pryr::mem_used(), silent = TRUE)
    if (!inherits(m, "try-error")) mb <- as.numeric(m) / (1024^2)
  }
  mb
}

measure <- function(expr) {
  gc()
  mem_before <- mem_used_mb()
  t0 <- now_time()
  val <- try(eval.parent(substitute(expr)), silent = TRUE)
  t1 <- now_time()
  mem_after <- mem_used_mb()
  list(value = val,
       time_sec = as.numeric(t1 - t0),
       mem_mb = if (is.na(mem_before) || is.na(mem_after)) NA_real_ else max(0, mem_after - mem_before))
}

## Confusion matrix flattening -------------------------------------------------
flatten_cm <- function(y_true, y_pred, levels_ref = NULL) {
  y_true <- factor(y_true)
  if (is.null(levels_ref)) levels_ref <- levels(y_true)
  y_pred <- factor(y_pred, levels = levels_ref)
  cm <- table(truth = y_true, pred = y_pred)
  # Ensure full grid
  grid <- expand.grid(truth = levels_ref, pred = levels_ref, stringsAsFactors = FALSE)
  out <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    tr <- grid$truth[i]; pr <- grid$pred[i]
    out[i] <- if (!is.null(cm[tr, pr])) cm[tr, pr] else 0
  }
  names(out) <- paste0("cm_", grid$truth, "_", grid$pred)
  out
}

## Fit/predict/diagnostics runner --------------------------------------------
run_one_case <- function(X_train, y_train, X_test, y_test, case) {
  # Build control via sdt for external engines, or tbc via sdt_tbc engine
  # Fit
  ctrl <- do.call(sdt_control, c(list(engine = case$engine, engine_args = case$engine_args), case$control_overrides))
  fit_meas <- measure({ sdt(X_train, y_train, control = ctrl) })

  # Assemble fit row data
  diag_vals <- list()
  diag_err <- NULL
  ensure_stats_fn()
  if (!inherits(fit_meas$value, "try-error")) {
    train_df <- data.frame(.y = y_train, X_train, check.names = FALSE)
    dtry <- try(get_tree_rule_stats(fit_meas$value$model, engine = case$engine, train_data = train_df), silent = TRUE)
    if (!inherits(dtry, "try-error")) diag_vals <- dtry else diag_err <- conditionMessage(attr(dtry, "condition") %||% dtry)
  } else {
    diag_err <- conditionMessage(attr(fit_meas$value, "condition") %||% fit_meas$value)
  }

  # Predict
  pred_meas <- list(value = NULL, time_sec = NA_real_, mem_mb = NA_real_)
  pred_vec <- NULL
  pred_err <- NULL
  if (!inherits(fit_meas$value, "try-error")) {
    pred_meas <- measure({ predict(fit_meas$value, X_test, type = "class") })
    if (!inherits(pred_meas$value, "try-error")) {
      pred_vec <- as.factor(pred_meas$value)
    } else {
      pred_err <- conditionMessage(attr(pred_meas$value, "condition") %||% pred_meas$value)
    }
  }

  list(fit = fit_meas, diag_vals = diag_vals, diag_err = diag_err,
       pred = pred_meas, pred_vec = pred_vec, pred_err = pred_err)
}

## Row assembler --------------------------------------------------------------
assemble_rows <- function(dataset, case, repeat_id, fold_id, res, levels_ref, y_test = NULL) {
  # Common fields
  base <- list(
    dataset   = dataset,
    algorithm = case$algorithm,
    criterion = case$criterion,
    `repeat`  = repeat_id,
    fold      = fold_id
  )

  # Diagnostics columns (prefix diag_)
  diag_cols <- c(
    diag_min_depth = NA_real_,
    diag_avg_depth = NA_real_,
    diag_max_depth = NA_real_,
    diag_avg_weighted_depth = NA_real_,
    diag_decision_nodes = NA_real_,
    diag_leaf_nodes = NA_real_,
    diag_unique_vars = NA_real_
  )
  if (is.list(res$diag_vals) && length(res$diag_vals)) {
    dv <- res$diag_vals
    diag_cols["diag_min_depth"]          <- dv$min_depth %||% NA
    diag_cols["diag_avg_depth"]          <- dv$avg_depth %||% NA
    diag_cols["diag_max_depth"]          <- dv$max_depth %||% NA
    diag_cols["diag_avg_weighted_depth"] <- dv$avg_weighted_depth %||% NA
    diag_cols["diag_decision_nodes"]     <- dv$decision_nodes %||% NA
    diag_cols["diag_leaf_nodes"]         <- dv$leaf_nodes %||% NA
    diag_cols["diag_unique_vars"]        <- dv$unique_vars %||% NA
  }

  # Fit row
  fit_row <- c(
    list(mode = "fit"),
    base,
    list(time_sec = res$fit$time_sec, mem_mb = res$fit$mem_mb, error = if (!is.null(res$diag_err)) as.character(res$diag_err) else ""),
    as.list(diag_cols)
  )

  # Predict row (confusion matrix computed if y_test provided and predictions available)
  pred_row <- c(
    list(mode = "predict"),
    base,
    list(time_sec = res$pred$time_sec, mem_mb = res$pred$mem_mb, error = if (!is.null(res$pred_err)) as.character(res$pred_err) else "")
  )
  list(fit_row = fit_row, pred_row = pred_row)
}

## CSV writer ---------------------------------------------------------------
write_rows <- function(fp, rows_df, append = FALSE) {
  utils::write.table(rows_df, file = fp, sep = ",", row.names = FALSE, col.names = !append, append = append, qmethod = "double")
}

## Main run -------------------------------------------------------------------
run_all <- function(data_dir = file.path(getwd(), "data"), out_dir = file.path(getwd(), "tools", "benchmark"),
                    k = 10L, repeats = 100L, seed = 20251210L, verbose = TRUE, n_workers = -1L) {
  load_TBC(); ensure_stats_fn(); load_required_packages()
  datasets <- list_lrn_cls_datasets(data_dir)
  if (length(datasets) == 0L) stop("No .LRN/.CLS pairs found in ", data_dir)
  cases <- engine_grid()
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Decide on workers and set a future plan once per run, if available
  # change / remove if different hardware used
  use_future <- requireNamespace("future.apply", quietly = TRUE) && requireNamespace("future", quietly = TRUE)
  workers_used <- 1L
  if (use_future) {
    # Force a specific worker count by default so no console setup is needed.
    # Priority: explicit n_workers argument -> env TBC_CV_WORKERS -> hard default 28
    workers <- n_workers
    if (is.null(workers) || is.na(workers) || workers == -1L) {
      env_workers <- suppressWarnings(as.integer(Sys.getenv("TBC_CV_WORKERS", unset = "")))
      if (!is.na(env_workers) && length(env_workers) == 1L) {
        workers <- env_workers
      } else {
        workers <- 28L
      }
    }
    workers_used <- as.integer(workers)
    # Always set/override plan explicitly to ensure desired concurrency
    future::plan(future::multisession, workers = workers_used)
    message(sprintf("Parallel backend: future multisession with %d workers (forced)", workers_used))
  } else {
    message("Parallel backend: sequential (1 worker); install packages 'future' and 'future.apply' to enable parallelism")
  }

  # Sort datasets by size of LRN data (fastest first)
  get_dataset_size <- function(ds) {
    # Prefer number of rows in LRN; fall back to file size if reader not available
    n <- NA_real_
    # Try quick read of LRN header/data
    n <- try({
      lobj <- suppressWarnings(suppressMessages(
        ReadLRN(FileName = basename(ds$lrn), InDirectory = dirname(ds$lrn))
      ))
      if (!is.null(lobj$Data)) nrow(lobj$Data) else NA_real_
    }, silent = TRUE)
    if (inherits(n, "try-error") || is.na(n)) {
      fi <- try(file.info(ds$lrn)$size, silent = TRUE)
      if (!inherits(fi, "try-error")) return(as.numeric(fi))
      return(Inf)
    }
    as.numeric(n)
  }

  sizes <- vapply(datasets, get_dataset_size, numeric(1))
  ord <- order(sizes, na.last = TRUE)
  datasets <- datasets[ord]
  # Informative listing
  message("Dataset run order (fastest first):")
  for (i in seq_along(datasets)) {
    nm <- datasets[[i]]$name
    sz <- sizes[ord][i]
    msg <- if (is.finite(sz)) sprintf("  %2d) %s  (rows or size: %g)", i, nm, sz) else sprintf("  %2d) %s", i, nm)
    message(msg)
  }

  # Resume/start from a specific dataset name (case-insensitive) — default 'coimbra'
  # Remove/change if run on different datasets, helps when something crashes to pick up later
  start_from_name <- Sys.getenv("TBC_CV_START_AT", unset = "coimbra")
  if (!is.null(start_from_name) && nzchar(start_from_name)) {
    names_lower <- vapply(datasets, function(ds) tolower(ds$name), character(1))
    pos <- match(tolower(start_from_name), names_lower)
    if (!is.na(pos)) {
      if (pos > 1L) {
        message(sprintf("Starting from dataset '%s' at position %d of %d", datasets[[pos]]$name, pos, length(datasets)))
        datasets <- datasets[pos:length(datasets)]
      }
    }
  }

  for (ds in datasets) {
    ds_name <- ds$name
    hr(paste0("dataset: ", ds_name))
    dat <- try(read_lrn_cls(ds$lrn, ds$cls, data_dir = dirname(ds$lrn)), silent = TRUE)
    if (inherits(dat, "try-error")) {
      message("ERROR reading ", ds_name, ": ", conditionMessage(attr(dat, "condition") %||% dat))
      next
    }
    X <- dat$X; y <- dat$y
    levels_ref <- levels(y)
    folds_rep <- make_stratified_folds(y, k = k, repeats = repeats, seed = seed)
    out_fp <- file.path(out_dir, paste0(ds_name, "_cv_benchmark.csv"))
    wrote_header <- FALSE

    # Precompute a task list across all repeats, folds, and cases
    tasks <- list()
    for (rep_id in seq_len(repeats)) {
      folds <- folds_rep[[rep_id]]
      for (fold_id in seq_len(k)) {
        tasks[[length(tasks) + 1L]] <- list(rep_id = rep_id, fold_id = fold_id)
      }
    }

    # Helper to process one (rep, fold, case)
    process_case <- function(rep_id, fold_id, case) {
      folds <- folds_rep[[rep_id]]
      test_idx <- folds[[fold_id]]
      train_idx <- setdiff(seq_len(nrow(X)), test_idx)
      X_train <- X[train_idx, , drop = FALSE]
      y_train <- droplevels(y[train_idx])
      X_test  <- X[test_idx, , drop = FALSE]
      y_test  <- droplevels(y[test_idx])

      cm_template <- as.list(setNames(rep(0, length(levels_ref)^2), paste0("cm_", rep(levels_ref, each = length(levels_ref)), "_", rep(levels_ref, times = length(levels_ref)))))

      res <- try(run_one_case(X_train, y_train, X_test, y_test, case), silent = TRUE)
      if (inherits(res, "try-error")) {
        err_msg <- conditionMessage(attr(res, "condition") %||% res)
        return(.error_rows_for_case(ds_name, rep_id, fold_id, case, levels_ref, n_train = length(y_train), n_test = length(y_test), err = err_msg))
      }
      fit_row <- c(
        list(mode = "fit", dataset = ds_name, algorithm = case$algorithm, criterion = case$criterion, `repeat` = rep_id, fold = fold_id,
             n_train = length(y_train), n_test = length(y_test), time_sec = res$fit$time_sec, mem_mb = res$fit$mem_mb,
             error = if (!is.null(res$diag_err)) as.character(res$diag_err) else ""),
        if (is.list(res$diag_vals) && length(res$diag_vals)) {
          list(
            diag_min_depth = res$diag_vals$min_depth %||% NA,
            diag_avg_depth = res$diag_vals$avg_depth %||% NA,
            diag_max_depth = res$diag_vals$max_depth %||% NA,
            diag_avg_weighted_depth = res$diag_vals$avg_weighted_depth %||% NA,
            diag_decision_nodes = res$diag_vals$decision_nodes %||% NA,
            diag_leaf_nodes = res$diag_vals$leaf_nodes %||% NA,
            diag_unique_vars = res$diag_vals$unique_vars %||% NA
          )
        } else list(
          diag_min_depth = NA, diag_avg_depth = NA, diag_max_depth = NA,
          diag_avg_weighted_depth = NA, diag_decision_nodes = NA, diag_leaf_nodes = NA, diag_unique_vars = NA
        ),
        cm_template
      )
      cm_vals <- cm_template
      pred_error <- if (!is.null(res$pred_err)) as.character(res$pred_err) else ""
      if (!is.null(res$pred_vec)) {
        cmf <- flatten_cm(y_test, res$pred_vec, levels_ref)
        for (nm in names(cmf)) cm_vals[[nm]] <- as.integer(cmf[[nm]])
      }
      pred_row <- c(
        list(mode = "predict", dataset = ds_name, algorithm = case$algorithm, criterion = case$criterion, `repeat` = rep_id, fold = fold_id,
             n_train = length(y_train), n_test = length(y_test), time_sec = res$pred$time_sec, mem_mb = res$pred$mem_mb,
             error = pred_error),
        if (is.list(res$diag_vals) && length(res$diag_vals)) {
          list(
            diag_min_depth = res$diag_vals$min_depth %||% NA,
            diag_avg_depth = res$diag_vals$avg_depth %||% NA,
            diag_max_depth = res$diag_vals$max_depth %||% NA,
            diag_avg_weighted_depth = res$diag_vals$avg_weighted_depth %||% NA,
            diag_decision_nodes = res$diag_vals$decision_nodes %||% NA,
            diag_leaf_nodes = res$diag_vals$leaf_nodes %||% NA,
            diag_unique_vars = res$diag_vals$unique_vars %||% NA
          )
        } else list(
          diag_min_depth = NA, diag_avg_depth = NA, diag_max_depth = NA,
          diag_avg_weighted_depth = NA, diag_decision_nodes = NA, diag_leaf_nodes = NA, diag_unique_vars = NA
        ),
        cm_vals
      )
      list(fit_row = fit_row, pred_row = pred_row)
    }

    # Build the full list of (rep, fold, case) tasks
    rf_pairs <- tasks
    # function that takes one rf pair and returns list of results for all cases
    process_rf_pair <- function(rf) {
      rep_id <- rf$rep_id; fold_id <- rf$fold_id
      # Ensure worker session has what it needs (multisession launches clean R sessions)
      # This is lightweight and only runs once per task on each worker.
      # Also ensure working directory points to project root so relative paths work.
      try(load_TBC(), silent = TRUE)
      try(ensure_stats_fn(), silent = TRUE)
      try(load_required_packages(), silent = TRUE)
      try(setwd("", warn = FALSE), silent = TRUE)  # no-op if not allowed
      tryCatch({
        lapply(cases, function(case) process_case(rep_id, fold_id, case))
      }, error = function(e) {
        # Build error rows for all cases to keep CSV structure complete
        msg <- conditionMessage(e)
        message(sprintf("Error in dataset %s rep=%d fold=%d: %s", ds_name, rep_id, fold_id, msg))
        folds <- folds_rep[[rep_id]]
        test_idx <- folds[[fold_id]]
        train_idx <- setdiff(seq_len(nrow(X)), test_idx)
        y_train <- droplevels(y[train_idx]); y_test <- droplevels(y[test_idx])
        lapply(cases, function(case) .error_rows_for_case(ds_name, rep_id, fold_id, case, levels_ref, n_train = length(y_train), n_test = length(y_test), err = msg))
      })
    }

    # Execute in parallel across rf pairs; within each, iterate cases
    rf_results <- if (use_future) {
      future.apply::future_lapply(rf_pairs, process_rf_pair, future.seed = TRUE)
    } else {
      lapply(rf_pairs, process_rf_pair)
    }

    # Flatten and write results
    cols_order <- c(
      "mode","dataset","algorithm","criterion","repeat","fold","n_train","n_test","time_sec","mem_mb","error",
      "diag_min_depth","diag_avg_depth","diag_max_depth","diag_avg_weighted_depth","diag_decision_nodes","diag_leaf_nodes","diag_unique_vars",
      paste0("cm_", rep(levels_ref, each = length(levels_ref)), "_", rep(levels_ref, times = length(levels_ref)))
    )

    for (rf in rf_results) {
      for (rr in rf) {
        fit_df  <- as.data.frame(as.list(rr$fit_row), optional = TRUE, stringsAsFactors = FALSE)
        pred_df <- as.data.frame(as.list(rr$pred_row), optional = TRUE, stringsAsFactors = FALSE)
        for (nm in cols_order) { if (!nm %in% names(fit_df))  fit_df[[nm]]  <- NA }
        for (nm in cols_order) { if (!nm %in% names(pred_df)) pred_df[[nm]] <- NA }
        fit_df  <- fit_df[ , cols_order]
        pred_df <- pred_df[, cols_order]
        append_mode <- wrote_header
        write_rows(out_fp, fit_df,  append = append_mode)
        write_rows(out_fp, pred_df, append = TRUE)
        wrote_header <- TRUE
        if (verbose) cat(sprintf("%s rep %d fold %d: %s/%s done\n",
                                 ds_name, fit_df[["repeat"]][1], fit_df$fold[1], fit_df$algorithm[1], fit_df$criterion[1]))
      }
    }
    cat("Wrote:", out_fp, "\n")

    # --- Cleanup after each dataset to reduce memory footprint -------------
    # Drop large objects and trigger garbage collection.
    # Note: worker processes (multisession) manage their own memory; this clears master.
    rm(list = c(
      "X", "y", "levels_ref", "folds_rep", "tasks", "rf_pairs",
      "rf_results", "process_case", "process_rf_pair"
    ), envir = environment())
    invisible(gc())
  }
}

## Run if sourced interactively ----------------------------------------------
if (identical(environment(), .GlobalEnv)) {
  # Defaults: allow override via env vars
  k <- as.integer(Sys.getenv("TBC_CV_K", unset = 10L))
  repeats <- as.integer(Sys.getenv("TBC_CV_REPEATS", unset = 100L))
  seed <- as.integer(Sys.getenv("TBC_CV_SEED", unset = 20251210L))
  n_workers <- as.integer(Sys.getenv("TBC_CV_WORKERS", unset = -1L))
  run_all(k = k, repeats = repeats, seed = seed, n_workers = n_workers)
}
