# Access training data and cell trajectories from the saved melanoma object.
trajectory_training_view <- function(melanoma, exclude_proliferation = TRUE) {
  training <- melanoma@misc$trajectory_training
  if (!exclude_proliferation) return(training)
  keep <- training$metadata$melanoma_mp != 'proliferation_cell_cycle'
  training$metadata <- training$metadata[keep, , drop = FALSE]
  training$features <- setdiff(training$features, 'proliferation_cell_cycle')
  training$coordinates <- destiny::eigenvectors(melanoma@misc$trajectory_models$diffusion_map)[, 1:3, drop = FALSE]
  training$signs <- melanoma@misc$trajectory_models$signs
  training
}

trajectory_cells <- function(melanoma, stage = c('corrected', 'raw', 'aligned')) {
  stage <- match.arg(stage)
  fields <- c('DM1', 'DM2', 'DM3', 'pt_route1_raw', 'pt_route2_raw', 'w1', 'w2',
              'slingshot_lineage_assignment', 'pt_route1', 'pt_route2', 'pt_merged',
              'delta', 'pt_class', 'cl_dist', 'cl_far')
  metadata <- melanoma[[]]
  result <- data.frame(cell_id = rownames(metadata), metadata[, fields, drop = FALSE])
  if (stage == 'raw') return(result[, c('cell_id', fields[1:8]), drop = FALSE])
  if (stage == 'aligned') {
    fields <- c('pt_route1', 'pt_route2', 'pt_merged', 'delta', 'pt_class')
    result[, fields] <- metadata[, paste0(fields, '_aligned'), drop = FALSE]
  }
  result
}

make_centerline <- function(df, pt_class_use, sort_var, n_bins = 80, min_cells = 20, span = 0.25) {
  d <- df %>% filter(pt_class == pt_class_use, is.finite(pt_route1), is.finite(pt_route2), is.finite(.data[[sort_var]])) %>% mutate(bin = ntile(.data[[sort_var]], n_bins))
  d <- d %>% group_by(bin) %>% summarise(pt_route1 = median(pt_route1, na.rm = TRUE), pt_route2 = median(pt_route2, na.rm = TRUE), n = n(), .groups = 'drop') %>% filter(n >= min_cells) %>% arrange(bin)
  if (nrow(d) >= 6) {
    d <- d %>% mutate(idx = row_number(), pt_route1_s = predict(loess(pt_route1 ~ idx, span = span, degree = 2)), pt_route2_s = predict(loess(pt_route2 ~ idx, span = span, degree = 2)))
  } else {
    d <- d %>% mutate(pt_route1_s = pt_route1, pt_route2_s = pt_route2)
  }
  d %>% mutate(pt_class = pt_class_use)
}

nearest_on_skeleton <- function(points, skeleton) {
  m <- nrow(skeleton)
  a <- skeleton[-m, , drop = FALSE]
  b <- skeleton[-1, , drop = FALSE]
  ab <- b - a
  ab2 <- rowSums(ab^2)
  best_distance <- rep(Inf, nrow(points))
  s1 <- s2 <- numeric(nrow(points))
  for (j in seq_len(m - 1)) {
    apx <- points[, 1] - a[j, 1]
    apy <- points[, 2] - a[j, 2]
    t <- pmin(pmax((apx * ab[j, 1] + apy * ab[j, 2]) / ab2[j], 0), 1)
    x <- a[j, 1] + t * ab[j, 1]
    y <- a[j, 2] + t * ab[j, 2]
    distance <- sqrt((points[, 1] - x)^2 + (points[, 2] - y)^2)
    update <- distance < best_distance
    best_distance[update] <- distance[update]
    s1[update] <- x[update]
    s2[update] <- y[update]
  }
  data.frame(pt1_line = s1, pt2_line = s2)
}

correct_trajectory_outliers <- function(correction, retained = NULL, seed = 1234) {
  df <- correction$metadata
  out <- which(df$cl_far)
  branch <- ifelse(df$pt_route1[out] >= df$pt_route2[out], 'Lineage 1', 'Lineage 2')
  if (is.null(retained)) set.seed(seed)
  for (lineage in c('Lineage 1', 'Lineage 2')) {
    selected <- out[branch == lineage]
    if (!length(selected)) next
    s <- correction$skeleton[correction$skeleton$pt_class == lineage, ]
    key_cell <- if (lineage == 'Lineage 1') df$pt_route1[selected] else df$pt_route2[selected]
    key_line <- if (lineage == 'Lineage 1') s$pt_route1_s else s$pt_route2_s
    nearest <- vapply(key_cell, function(value) which.min(abs(key_line - value)), integer(1))
    if (is.null(retained)) {
      df$pt_route1[selected] <- s$pt_route1_s[nearest] + rnorm(length(selected), 0, correction$noise_sd)
      df$pt_route2[selected] <- s$pt_route2_s[nearest] + rnorm(length(selected), 0, correction$noise_sd)
    } else {
      index <- match(df$cell_id[selected], retained$cell_id)
      df$pt_route1[selected] <- retained$pt_route1[index]
      df$pt_route2[selected] <- retained$pt_route2[index]
    }
    df$pt_class[selected] <- lineage
  }
  df$pt_merged <- (df$pt_route1 + df$pt_route2) / 2
  df$delta <- df$pt_route1 - df$pt_route2
  df
}

# PDF export shared by the Xenium analysis notebooks.
save_panel <- function(plot, name, width = 4, height = 4) {
    path <- file.path(result_dir, paste0(name, ".pdf"))
    ggsave(path, plot, width = width, height = height, useDingbats = FALSE)
}
