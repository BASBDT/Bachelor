#' @keywords internal
.tbc_conf_mat <- function(y_true, y_pred) {
  # Ensure both factors share the same levels union
  lv <- union(levels(factor(y_true)), levels(factor(y_pred)))
  y_true <- factor(y_true, levels = lv)
  y_pred <- factor(y_pred, levels = lv)

  tab <- table(y_true, y_pred)
  # Attach levels for downstream logic if needed
  attr(tab, "levels") <- lv
  tab
}

#' Accuracy
#' @keywords internal
.tbc_acc <- function(cm) {
  total <- sum(cm)
  if (total == 0L) return(NA_real_)
  sum(diag(cm)) / total
}

#' Balanced Accuracy (Macro-Averaged Recall)
#' @keywords internal
.tbc_bacc <- function(cm) {
  n_classes <- nrow(cm)
  recalls   <- numeric(n_classes)

  for (i in seq_len(n_classes)) {
    tp    <- cm[i, i]
    total <- sum(cm[i, , drop = TRUE]) # Row sum (True positives + False negatives)
    recalls[i] <- if (total > 0) tp / total else NA_real_
  }

  mean(recalls, na.rm = TRUE)
}

#' Macro-F1 Score
#' @keywords internal
.tbc_macro_f1 <- function(cm) {
  n_classes <- nrow(cm)
  f1_scores <- numeric(n_classes)

  for (i in seq_len(n_classes)) {
    tp <- cm[i, i]
    fp <- sum(cm[, i, drop = TRUE]) - tp # Col sum - TP
    fn <- sum(cm[i, , drop = TRUE]) - tp # Row sum - TP

    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

    if (is.na(prec) || is.na(rec) || (prec + rec) == 0) {
      f1_scores[i] <- NA_real_
    } else {
      f1_scores[i] <- 2 * prec * rec / (prec + rec)
    }
  }

  mean(f1_scores, na.rm = TRUE)
}

#' Matthews Correlation Coefficient (Multiclass)
#' @keywords internal
.tbc_mcc <- function(cm) {
  # Based on Gorodkin (2004) / R_K statistic for multiclass confusion matrix
  c_sum <- sum(diag(cm))
  s_sum <- sum(cm)

  # Row sums (true counts per class) and Col sums (predicted counts per class)
  t_k <- rowSums(cm)
  p_k <- colSums(cm)

  numerator <- (c_sum * s_sum) - sum(t_k * p_k)

  denom_sq_1 <- (s_sum^2) - sum(p_k^2)
  denom_sq_2 <- (s_sum^2) - sum(t_k^2)

  if (denom_sq_1 <= 0 || denom_sq_2 <= 0) return(NA_real_)

  numerator / sqrt(denom_sq_1 * denom_sq_2)
}

#' Calculate Metric Vector
#'
#' @description
#' Convenience wrapper to compute a set of metrics. Accepts either a confusion
#' matrix directly, or truth/pred vectors.
#'
#' @keywords internal
.tbc_metric_vec <- function(cm, pred = NULL) {
  if (!is.null(pred)) {
    # Called as .tbc_metric_vec(truth, pred)
    cm <- .tbc_conf_mat(y_true = cm, y_pred = pred)
  }

  c(
    acc      = .tbc_acc(cm),
    bacc     = .tbc_bacc(cm),
    macro_f1 = .tbc_macro_f1(cm),
    mcc      = .tbc_mcc(cm)
  )
}