# Functions shared by the Visium analysis notebooks.

make_centerline <- function(df, sort_var, pt_class_use = NULL, n_bins = 70, min_cells = 20, span = 0.25) {
    d <- df
    if (!is.null(pt_class_use)) 
        d <- d %>% filter(pt_class == pt_class_use)
    d <- d %>% filter(is.finite(pt_route1), is.finite(pt_route2), is.finite(.data[[sort_var]])) %>% mutate(bin = ntile(.data[[sort_var]], 
        n_bins))
    d <- d %>% group_by(bin) %>% summarise(pt_route1 = median(pt_route1), pt_route2 = median(pt_route2), n = n(), 
        .groups = "drop") %>% filter(n >= min_cells) %>% arrange(bin)
    d$pt_route1_s <- d$pt_route1
    d$pt_route2_s <- d$pt_route2
    if (nrow(d) >= 6) {
        d$idx <- seq_len(nrow(d))
        d$pt_route1_s <- predict(loess(pt_route1 ~ idx, data = d, span = span, degree = 2))
        d$pt_route2_s <- predict(loess(pt_route2 ~ idx, data = d, span = span, degree = 2))
    }
    if (!is.null(pt_class_use)) 
        d$pt_class <- pt_class_use
    d
}

save_panel <- function(plot, filename, width, height) {
    path <- file.path(figure_dir, filename)
    ggsave(path, plot = plot, width = width, height = height, units = "in", limitsize = FALSE)
}

assemble_section_expression <- function(sections, assay = 'SCT', layer = 'data', genes = NULL, spots = NULL) {
    if (is.null(genes)) genes <- Reduce(intersect, lapply(sections, function(seu) rownames(seu[[assay]])))
    blocks <- lapply(sections, function(seu) {
        values <- LayerData(seu[[assay]], layer = layer)
        selected <- if (is.null(spots)) colnames(values) else intersect(colnames(values), spots)
        values[genes, selected, drop = FALSE]
    })
    expression <- do.call(cbind, blocks)
    if (!is.null(spots)) expression <- expression[, spots, drop = FALSE]
    expression
}

assign_spot_meta <- function(seu, data, columns) {
    index <- match(colnames(seu), data$spot)
    for (column in columns) seu[[column]] <- setNames(data[[column]][index], colnames(seu))
    seu
}
