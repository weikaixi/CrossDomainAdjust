# Feature encoders are fitted on training samples and then frozen.  These
# helpers intentionally have a module-specific prefix.
fe_numeric_matrix <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("Every feature column must be numeric.", call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L) {
    stop("x must be a nonempty numeric matrix, data frame, or named numeric vector.",
         call. = FALSE)
  }
  feature_names <- colnames(x)
  if (is.null(feature_names) || anyNA(feature_names) ||
      any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
    stop("x must have unique, nonempty feature names.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x must contain only finite values; missing values are not supported.",
         call. = FALSE)
  }
  x
}

fe_positive_integer <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop(name, " must be a positive integer", if (allow_null) " or NULL" else "",
         ".", call. = FALSE)
  }
  x
}

fe_pair_indices <- function(pairs, feature_names) {
  if (is.data.frame(pairs)) pairs <- as.matrix(pairs)
  if (!is.matrix(pairs) || !is.character(pairs) || ncol(pairs) != 2L ||
      nrow(pairs) < 1L || anyNA(pairs) || any(!nzchar(pairs))) {
    stop("pairs must be a nonempty two-column character matrix or data frame of gene names.",
         call. = FALSE)
  }
  first <- match(pairs[, 1L], feature_names)
  second <- match(pairs[, 2L], feature_names)
  if (anyNA(first) || anyNA(second)) {
    stop("All genes named in pairs must occur in the training feature names.",
         call. = FALSE)
  }
  if (any(first == second)) stop("Self pairs are not allowed.", call. = FALSE)
  pair_first <- pmin(first, second)
  pair_second <- pmax(first, second)
  if (anyDuplicated(paste(pair_first, pair_second, sep = ":"))) {
    stop("Duplicate or reversed redundant pairs are not allowed.", call. = FALSE)
  }
  list(first = pair_first, second = pair_second)
}

#' Fit a frozen feature representation
#'
#' @export
fit_feature_encoder <- function(x, representation = c("raw", "rank", "pair"),
                                frequency = c(0.2, 0.8), pairs = NULL,
                                max_pairs = NULL, max_candidates = 1000000,
                                chunk_size = 1000) {
  x <- fe_numeric_matrix(x)
  representation <- match.arg(representation)
  if (!is.numeric(frequency) || length(frequency) != 2L ||
      any(!is.finite(frequency)) || frequency[1L] < 0 || frequency[2L] > 1 ||
      frequency[1L] > frequency[2L]) {
    stop("frequency must contain ordered lower and upper bounds within [0, 1].",
         call. = FALSE)
  }
  max_pairs <- fe_positive_integer(max_pairs, "max_pairs", allow_null = TRUE)
  max_candidates <- fe_positive_integer(max_candidates, "max_candidates")
  chunk_size <- fe_positive_integer(chunk_size, "chunk_size")
  if (chunk_size > .Machine$integer.max) {
    stop("chunk_size exceeds the supported integer range.", call. = FALSE)
  }
  feature_names <- colnames(x)
  p <- ncol(x)
  if (representation != "raw" && p < 2L) {
    stop("rank and pair representations require at least two genes.", call. = FALSE)
  }
  if (representation != "pair" && !is.null(pairs)) {
    stop("Explicit pairs are only used with representation = 'pair'.", call. = FALSE)
  }
  encoder <- list(
    representation = representation, feature_names = feature_names,
    output_names = feature_names, pairs = NULL, frequencies = NULL,
    candidate_count = p, retained_count = p, filtered_count = 0,
    frequency = unname(frequency), frequency_pass_count = p,
    frequency_filtered_count = 0, capped_count = 0,
    max_pairs = max_pairs, training_samples = nrow(x),
    tie_rule = if (representation == "rank") "average" else if
      (representation == "pair") "zero" else NA_character_
  )
  if (representation != "pair") {
    return(structure(encoder, class = "feature_encoder"))
  }

  explicit <- !is.null(pairs)
  if (explicit) {
    pair_index <- fe_pair_indices(pairs, feature_names)
    candidate_count <- length(pair_index$first)
  } else {
    # Compute and check the count before constructing any quadratic object.
    candidate_count <- as.double(p) * (p - 1) / 2
  }
  if (candidate_count > max_candidates) {
    stop("Candidate pair count (", format(candidate_count, scientific = FALSE),
         ") exceeds max_candidates (", format(max_candidates, scientific = FALSE),
         "). Preselect genes, provide explicit pairs, or deliberately raise the limit.",
         call. = FALSE)
  }

  retained_blocks <- list()
  best <- list(first = integer(), second = integer(), frequency = numeric(),
               candidate = numeric())
  frequency_pass_count <- 0
  candidate_offset <- 0
  consume_block <- function(first, second) {
    f <- colMeans(x[, first, drop = FALSE] > x[, second, drop = FALSE])
    candidate <- candidate_offset + seq_along(first)
    candidate_offset <<- candidate_offset + length(first)
    keep <- f >= frequency[1L] & f <= frequency[2L]
    frequency_pass_count <<- frequency_pass_count + sum(keep)
    if (!any(keep)) return(invisible(NULL))
    block <- list(first = first[keep], second = second[keep],
                  frequency = f[keep], candidate = candidate[keep])
    if (is.null(max_pairs)) {
      retained_blocks[[length(retained_blocks) + 1L]] <<- block
    } else {
      joined <- Map(c, best, block)
      # Prefer balanced pairs; ties are resolved by their original enumeration.
      selected <- order(abs(joined$frequency - 0.5), joined$candidate)
      selected <- utils::head(selected, max_pairs)
      best <<- lapply(joined, `[`, selected)
    }
    invisible(NULL)
  }

  if (explicit) {
    for (start in seq.int(1, candidate_count, by = chunk_size)) {
      at <- seq.int(start, min(candidate_count, start + chunk_size - 1))
      consume_block(pair_index$first[at], pair_index$second[at])
    }
  } else {
    # Only one comparison block is materialized, never the n x all-pairs matrix.
    for (first in seq_len(p - 1L)) {
      for (start in seq.int(first + 1L, p, by = chunk_size)) {
        second <- seq.int(start, min(p, start + chunk_size - 1))
        consume_block(rep.int(first, length(second)), second)
      }
    }
  }
  if (frequency_pass_count == 0) {
    stop("No pair features remain within the requested training frequency bounds.",
         call. = FALSE)
  }
  if (is.null(max_pairs)) {
    best <- lapply(names(best), function(nm) {
      unlist(lapply(retained_blocks, `[[`, nm), use.names = FALSE)
    })
    names(best) <- c("first", "second", "frequency", "candidate")
  }
  # Keep output order independent of block size and of the cap-selection loop.
  at <- order(best$candidate)
  best <- lapply(best, `[`, at)
  output_names <- sprintf("pair_%06.0f", best$candidate)
  retained_pairs <- cbind(first = feature_names[best$first],
                          second = feature_names[best$second])
  rownames(retained_pairs) <- output_names
  encoder$output_names <- output_names
  encoder$pairs <- retained_pairs
  encoder$frequencies <- stats::setNames(best$frequency, output_names)
  encoder$candidate_indices <- best$candidate
  encoder$candidate_count <- candidate_count
  encoder$retained_count <- length(output_names)
  encoder$filtered_count <- candidate_count - length(output_names)
  encoder$frequency_pass_count <- frequency_pass_count
  encoder$frequency_filtered_count <- candidate_count - frequency_pass_count
  encoder$capped_count <- frequency_pass_count - length(output_names)
  encoder$chunk_size <- chunk_size
  structure(encoder, class = "feature_encoder")
}

#' Apply a fitted feature encoder to new samples
#'
#' @export
encode_features <- function(encoder, x) {
  if (!inherits(encoder, "feature_encoder")) {
    stop("encoder must be produced by fit_feature_encoder().", call. = FALSE)
  }
  x <- fe_numeric_matrix(x)
  absent <- setdiff(encoder$feature_names, colnames(x))
  if (length(absent)) {
    stop("Missing fitted input genes: ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  x <- x[, encoder$feature_names, drop = FALSE]
  if (encoder$representation == "raw") return(x)
  if (encoder$representation == "rank") {
    p <- ncol(x)
    out <- t(vapply(seq_len(nrow(x)), function(i) {
      (rank(x[i, ], ties.method = "average") - 1) / (p - 1)
    }, numeric(p)))
    dimnames(out) <- list(rownames(x), encoder$output_names)
    return(out)
  }
  if (encoder$representation != "pair") {
    stop("Unrecognized representation in encoder.", call. = FALSE)
  }
  first <- match(encoder$pairs[, "first"], colnames(x))
  second <- match(encoder$pairs[, "second"], colnames(x))
  out <- matrix(0L, nrow(x), length(first),
                dimnames = list(rownames(x), encoder$output_names))
  # The retained output must exist, but comparison temporaries stay bounded.
  for (start in seq.int(1L, length(first), by = encoder$chunk_size)) {
    at <- seq.int(start, min(length(first), start + encoder$chunk_size - 1))
    out[, at] <- 1L * (x[, first[at], drop = FALSE] >
                        x[, second[at], drop = FALSE])
  }
  out
}
