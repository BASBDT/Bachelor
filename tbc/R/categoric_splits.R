#' Categorical split search (R) with exact mode for small L
#' Returns exactly the fields LearnPlanarDecision() expects.
#' CatBestCrit: numeric scalar (-Inf if no split)
#' CatBestVar:  integer scalar (0L if none)
#' CatBestCutA/B: list of character levels (or "0"/"0" when none)
#' @keywords internal
findCurrentBestCatDataCut <- function(
  CDIM,
  CatDataInNode,            # matrix (character) of categorical data for rows in node
  NodeXClsMembership,       # N x K logical/binary matrix
  priorRatio,               # length-K numeric
  sum_weighted_class_prob,  # unused, kept for signature stability
  ShuffleFeatures,          # logical
  Debug,                    # logical
  CatExactMaxLevels = 8L    # NEW: exact enumeration threshold
) {
  # default "no split" result
  no_split <- list(
    CatBestCrit = -Inf,
    CatBestVar  = 0L,
    CatBestCutA = "0",
    CatBestCutB = "0"
  )

  if (CDIM <= 0L || is.null(CatDataInNode) || ncol(as.matrix(CatDataInNode)) == 0L) {
    return(no_split)
  }

  Xc <- as.matrix(CatDataInNode)
  storage.mode(Xc) <- "character"

  # Permute categorical features if requested
  pcat  <- ncol(Xc)
  order <- if (isTRUE(ShuffleFeatures)) sample.int(pcat) else seq_len(pcat)

  # Parent risk (aligned to numeric/bivariate)
  K            <- length(priorRatio)
  y            <- max.col(NodeXClsMembership, ties.method = "first")
  parentCounts <- tabulate(y, nbins = K)
  parentMass   <- sum(priorRatio * parentCounts)
  if (!(is.finite(parentMass) && parentMass > 0)) return(no_split)
  parentMaxP   <- max((priorRatio * parentCounts) / parentMass)
  parentRisk   <- parentMass * (1 - parentMaxP)

  bestCrit <- -Inf
  bestVar  <- 0L
  bestA    <- "0"
  bestB    <- "0"

  for (jj in order) {
    xj   <- Xc[, jj]
    lev  <- sort(unique(xj))
    nlev <- length(lev)
    if (nlev < 2L) next

    # route-by-majority among KNOWN (tie -> left); training-time decision only
    # (actual prediction routing is per-node via TL$NA_to_left)
    route_majority <- function(leftMask, rightMask, unknown) {
      nLk <- sum(leftMask)
      nRk <- sum(rightMask)
      if (nLk == 0L || nRk == 0L) return(NULL) # empty child -> invalid candidate
      (nLk >= nRk)
    }

    eval_split <- function(A, B) {
      leftMask  <- xj %in% A
      rightMask <- xj %in% B
      unknown   <- !(leftMask | rightMask)

      route_left <- route_majority(leftMask, rightMask, unknown)
      if (is.null(route_left)) return(NULL)

      leftAll  <- tabulate(y[leftMask],  nbins = K)
      rightAll <- tabulate(y[rightMask], nbins = K)
      if (any(unknown)) {
        unkCounts <- tabulate(y[unknown], nbins = K)
        if (route_left) leftAll <- leftAll + unkCounts else rightAll <- rightAll + unkCounts
      }

      massL <- sum(priorRatio * leftAll)
      massR <- sum(priorRatio * rightAll)
      riskL <- if (massL > 0) massL * (1 - max((priorRatio * leftAll)  / massL)) else 0
      riskR <- if (massR > 0) massR * (1 - max((priorRatio * rightAll) / massR)) else 0
      parentRisk - (riskL + riskR)
    }

    if (nlev <= as.integer(CatExactMaxLevels)) {
      # ---- exact two-set enumeration (unique partitions) --------------------
      # masks from 1 .. 2^(nlev-1)-1 to avoid duplicates via complements
      Lhalf <- nlev - 1L
      for (m in seq_len(bitwShiftL(1L, Lhalf) - 1L)) {
        A <- lev[which(as.logical(intToBits(m))[seq_len(nlev)])]
        # ensure non-empty and not full (B must be non-empty)
        if (length(A) == 0L || length(A) == nlev) next
        B <- setdiff(lev, A)
        crit <- eval_split(A, B)
        if (!is.null(crit) && crit > bestCrit) {
          bestCrit <- crit
          bestVar  <- jj
          bestA    <- A
          bestB    <- B
        }
      }
    } else {
      # ---- heuristic: ordered sweep by class index mean ---------------------
      level_score <- tapply(y, factor(xj, levels = lev), mean, na.rm = TRUE)
      ord         <- order(level_score, na.last = NA)
      lev_ord     <- lev[ord]

      for (t in 1:(nlev - 1L)) {
        A <- lev_ord[1:t]
        B <- lev_ord[(t + 1L):nlev]
        crit <- eval_split(A, B)
        if (!is.null(crit) && crit > bestCrit) {
          bestCrit <- crit
          bestVar  <- jj
          bestA    <- A
          bestB    <- B
        }
      }
    }
  }

  if (!is.finite(bestCrit)) bestCrit <- -Inf
  if (bestVar == 0L) {
    return(no_split)
  } else {
    return(list(
      CatBestCrit = bestCrit,
      CatBestVar  = as.integer(bestVar),
      CatBestCutA = bestA,
      CatBestCutB = bestB
    ))
  }
}

#' Internal shim for legacy name used inside tbc_fit
#' @keywords internal
.tbc_find_best_cat_cut <- function(CDIM, CatDataInNode, NodeXClsMembership,
                                  priorRatio, sum_weighted_class_prob = NULL,
                                  ShuffleFeatures = TRUE, Debug = FALSE,
                                  CatExactMaxLevels = 8L) {
  findCurrentBestCatDataCut(
    CDIM = CDIM,
    CatDataInNode = CatDataInNode,
    NodeXClsMembership = NodeXClsMembership,
    priorRatio = priorRatio,
    sum_weighted_class_prob = sum_weighted_class_prob,
    ShuffleFeatures = ShuffleFeatures,
    Debug = Debug,
    CatExactMaxLevels = CatExactMaxLevels
  )
}
