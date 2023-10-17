asym_dtw_alignment <- function(query, reference, index_query, index_ref) {

  dtw_fit <- tryCatch({
    dtw::dtw(query, reference, step = dtw::asymmetric, open.begin = TRUE, open.end = TRUE, keep=TRUE)
  },
  error = function(e) {
    print(e)
    return(NULL)
  })

  if (is.null(dtw_fit)) {
    return(list(ndist = Inf, alignment = NA))
  }

  alignment_matrix <- cbind(dtw_fit$index1, dtw_fit$index2)
  alignment_matrix[, 1] <- index_query[alignment_matrix[, 1]]
  alignment_matrix[, 2] <- index_ref[alignment_matrix[, 2]]

  return(list(ndist = dtw_fit$normalizedDistance,
              alignment = alignment_matrix))
}

dtw_neighbors <- function(x, k, center=TRUE){
  centers <- lapply(x, function(xi) colMeans(xi[,-1,drop=FALSE]))

  x <- if (center){
    lapply(x, function(xi){
      xi[,-1] <- xi[,-1] - colMeans(xi[,-1,drop=FALSE])
      return(xi)
    })
  }
  dist_matrix <- matrix(nrow=length(x),ncol=length(x))
  alignments <- list()
  for (i in seq_along(x)){
    alignments <- c(alignments, list(list()))
    for (j in seq_along(x)){
      if (nrow(x[[j]]) >= nrow(x[[i]]) & i!=j){
        dtw_alignment <-
          asym_dtw_alignment(x[[i]][, -1], x[[j]][, -1], x[[i]][, 1], x[[j]][, 1])
        dist_matrix[i,j] <- dtw_alignment$ndist
        alignments[[i]] <- c(alignments[[i]], list(dtw_alignment$alignment))
      }
      else{
        dist_matrix[i,j] <- Inf
        alignments[[i]] <- c(alignments[[i]],list(NA))
      }
    }
  }

  neighbors <- sapply(seq_along(x), function(i){
    idx <- sort(dist_matrix[i,], index.return = TRUE)$ix
    n_idx <- idx[1:min(length(idx),k[i])]

    return(n_idx)
  }, simplify = FALSE)
  dists_neighbors <- sapply(seq_along(x), function(i) dist_matrix[i, neighbors[[i]]], simplify=FALSE)
  alignments_neighbors <- sapply(seq_along(x), function(i) alignments[[i]][neighbors[[i]]], simplify=FALSE)

  return(list(data=x,
              neighbors=neighbors,
              dists=dists_neighbors,
              alignments=alignments_neighbors,
              centers=centers))
}

adba <- function(data, init_cent, niter=10, return_steps=FALSE){
  result <- list(list(avg=init_cent))

  for (i in 1:niter) {
    prev_avg <- result[[length(result)]]$avg

    total_al_matrix <- data.frame()
    dist_to_avg <- numeric()
    for (d in data){
      ali_d <-
        asym_dtw_alignment(
          d[,-1],
          prev_avg[,-1],
          index_query = d[,1],
          index_ref = prev_avg[,1]
        )
      dist_to_avg <- c(dist_to_avg, ali_d$ndist)
      total_al_matrix <- if (!is.infinite(ali_d$ndist)){
        rbind(total_al_matrix, cbind(ali_d$alignment,d[d[, 1] %in% ali_d$alignment[, 1], -1]))
      }
    }

    avg_indices <- unique(total_al_matrix[,2])
    new_avg <- prev_avg
    for (j in avg_indices){
      new_avg[new_avg[,1] == j, -1] <- t(colMeans(total_al_matrix[total_al_matrix[,2] == j,-c(1,2),drop=FALSE]))
    }

    result[[i]] <- c(result[[i]], dist_to_avg=list(dist_to_avg))
    result <- c(result, list(list(avg=new_avg)))
  }

  # do one last alignment to get distances
  dist_to_avg <- numeric()
  for (d in data){
    ali_d <-
      asym_dtw_alignment(
        d[,-1],
        result[[length(result)]]$avg[,-1],
        index_query = d[,1],
        index_ref = result[[length(result)]]$avg[,1]
      )
    dist_to_avg <- c(dist_to_avg, ali_d$ndist)
  }
  result[[niter+1]] <- c(result[[niter+1]],dist_to_avg=list(dist_to_avg))

  if (return_steps == TRUE) {
    return(result)
  } else {
    return(result[[niter + 1]])
  }
}
