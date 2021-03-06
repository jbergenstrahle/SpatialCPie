#' @importFrom data.table as.data.table
#' @importFrom digest sha1
#' @importFrom dplyr
#' arrange filter first group_by inner_join mutate n rename select summarize
#' ungroup
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggplot2
#' aes_ aes_string coord_fixed element_blank geom_segment ggplot ggtitle guides
#' guide_legend
#' labs
#' theme theme_bw theme_minimal
#' scale_color_manual scale_fill_manual scale_size
#' scale_x_continuous scale_y_continuous
#' @importFrom grid unit
#' @importFrom methods is
#' @importFrom purrr
#' %>% %||% accumulate array_branch lift invoke keep map map_dbl map_int partial
#' reduce transpose
#' @importFrom readr read_file write_file
#' @importFrom rlang !! := .data sym
#' @importFrom shiny debounce observeEvent reactive
#' @importFrom shinyjs hideElement
#' @importFrom shinyWidgets radioGroupButtons materialSwitch
#' @importFrom stats dist kmeans setNames sd
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr gather separate spread unite
#' @importFrom tidyselect everything quo
#' @importFrom tools toTitleCase
#' @importFrom utils head tail
#' @importFrom utils str
#' @importFrom zeallot %<-%
"_PACKAGE"


## Pre-declare all NSE variables as global in order to appease R CMD check
## (ref: https://stackoverflow.com/a/12429344)
globalVariables(c(
    ".",
    "cluster",
    "count",
    "name",
    "otherMargin",
    "resolution",
    "spot",
    "xcoord",
    "ycoord",
    NULL
))


#' Logsumexp
#'
#' Adapted from https://stat.ethz.ch/pipermail/r-help/2011-February/269205.html
#' @param xs input vector
#' @return log of summed exponentials
#' @keywords internal
.logsumexp <- function(xs) {
    idx <- which.max(xs)
    log1p(sum(exp(xs[-idx] - xs[idx]))) + xs[idx]
}


#' Likeness score
#'
#' @param d distance vector.
#' @param c log multiplier.
#' @return vector of scores.
#' @keywords internal
.likeness <- function(
    d,
    c = 1.0
) {
    exp(-c * d - .logsumexp(-c * d))
}

#' Z-score
#'
#' @param xs vector of observations
#' @return `xs`, z-normalized. if all elements of `xs` are equal, a vector of
#'     zeros will be returned instead.
#' @keywords internal
.zscore <- function(xs) {
    std <- sd(xs)
    std <- if (std == 0.0) 1 else std
    (xs - mean(xs)) / std
}


#' Maximize overlap
#'
#' @param xss list of lists of labels.
#' @return `xss`, relabeled so as to maximize the overlap between labels in
#' consecutive label lists.
#' @keywords internal
.maximizeOverlap <- function(
    xss
) {
    maximumOverlap <- function(xs, ys) {
        setNames(nm = sort(unique(xs))) %>%
            map(function(x)
                setNames(nm = sort(unique(ys))) %>%
                    map_dbl(function(y) sum(`*`(xs == x, ys == y)))
            ) %>%
            invoke(rbind, .) %>%
            (function(overlaps) {
                all <- union(rownames(overlaps), colnames(overlaps))
                n <- length(all)

                ## Zero-pad overlap matrix so that all labels are represented in
                ## both the to and from dimensions
                paddedOverlaps <-
                    overlaps %>%
                    rbind(do.call(
                        rbind,
                        rep(list(rep(0, n)), n - nrow(overlaps))
                    )) %>%
                    cbind(do.call(
                        cbind,
                        rep(list(rep(0, n)), n - ncol(overlaps))
                    ))
                rownames(paddedOverlaps)[rownames(paddedOverlaps) == ""] <-
                    setdiff(all, rownames(paddedOverlaps))
                colnames(paddedOverlaps)[colnames(paddedOverlaps) == ""] <-
                    setdiff(all, colnames(paddedOverlaps))

                ## Solve the assignment problem to maximize the overlap
                lpSolve::lp.assign(-paddedOverlaps)$solution %>%
                    array_branch(2) %>%
                    map(~colnames(paddedOverlaps)[which.max(.)]) %>%
                    setNames(nm = rownames(paddedOverlaps))
            })
    }

    ## Convert cluster labels to natural numbers
    xss <- map(
        xss,
        function(x) setNames(as.character(as.integer(as.factor(x))), names(x))
    )

    ## Compute reassignment map between each label pair
    reassignments <-
        list(unname(head(xss, -1)), unname(tail(xss, -1))) %>%
        transpose %>%
        map(lift(maximumOverlap))

    ## Sync reassignments by propagating them forward
    reassignments <- accumulate(
        reassignments,
        function(prev, cur) {
            list(lapply(cur, function(x) {
                if (x %in% names(prev[[1]])) prev[[1]][[x]]
                else x
            }))
        },
        .init = list(setNames(nm = unique(xss[[1]])))
    )

    ## Apply reassignments
    list(xss, reassignments) %>%
        transpose() %>%
        map(lift(function(xs, reassignment) {
            vapply(xs, function(x) reassignment[[x]], character(1))
        })) %>%
        setNames(names(xss))
}


#' Tidy assignments
#'
#' @param assignments list of assignment vectors.
#' @return a \code{\link[base]{data.frame}} containing the `assignments`, with
#' the data relabeled so that the overlap between consecutive assignment
#' vectors is maximized. Additionally, a "root" resolution is added.
#' @keywords internal
.tidyAssignments <- function(
    assignments
) {
    if (length(assignments) == 0) {
        stop("Need at least one resolution")
    }

    ## Add "root" resolution
    units <- names(assignments[[1]])
    assignments <- c(
        list("root" = setNames(rep(1, length(assignments[[1]])), nm = units)),
        assignments
    )

    ## Relabel the data to maximize overlap between labels in consecutive
    ## resolutions
    message("Maximizing label overlap in consecutive resolutions")
    assignments <- .maximizeOverlap(assignments)

    ## Concatenate assignments to `data.frame`
    assignments <-
        list(names(assignments), assignments) %>%
        transpose() %>%
        map(lift(function(res, xs)
            data.frame(
                name = sprintf("resolution %s, cluster %s", res, xs),
                resolution = res,
                cluster = xs,
                stringsAsFactors = TRUE
            ) %>%
            tibble::rownames_to_column("unit")
        )) %>%
        invoke(rbind, .)

    assignments
}


#' Compute cluster colors
#'
#' Computes colors so that dissimilar clusters are far away in color space.
#' @param clusterMeans matrix of size `(n, K)` representing the `n` feature
#' means for each of the `K` clusters.
#' @return vector of cluster colors.
#' @keywords internal
.computeClusterColors <- function(
    clusterMeans
) {
    clusterLoadings <- stats::prcomp(
        t(clusterMeans),
        rank = 2,
        center = TRUE
    )$x
    minLoading <- apply(clusterLoadings, 2, min)
    maxLoading <- apply(clusterLoadings, 2, max)

    clusterColors <- cbind(
        50,
        200 * t(
            (t(clusterLoadings) - minLoading)
            / (maxLoading - minLoading + 1e-10))
        - 100
    )

    colorspace::LAB(clusterColors) %>%
        colorspace::hex(fixup = TRUE)
}


#' Preprocess data
#'
#' Preprocesses input data for \code{\link{.makeServer}}.
#' @param counts count matrix. `rownames` should correspond to genes and
#' `colnames` should correspond to spot coordinates.
#' @param margin which margin of the count matrix to cluster. Valid values are
#' `c("spot", "sample", "gene", "feature")`.
#' @param resolutions vector of resolutions to cluster.
#' @param assignmentFunction function to compute cluster assignments. The
#' function should have the following signature: integer (number of clusters) ->
#' (m, n) feature matrix -> m-length vector (cluster assignment of each data
#' point).
#' @param coordinates optional \code{\link[base]{data.frame}} with pixel
#' coordinates for each spot. `rownames` should correspond to the `colnames` of
#' `counts` and the columns `x` and `y` should specify the pixel coordinates of
#' the spots.
#' @return list with the following elements:
#' - `$assignments`: tidy assignments
#' - `$means`: cluster means
#' - `$scores`: cluster scores for each spot in each resolution
#' - `$colors`: cluster colors
#' - `$coordinates`: spot coordinates, either from `coordinates` or parsed from
#' `assignments`
#' - `$featureName`: name of the clustered feature (the "opposite" of `margin`)
#' @keywords internal
.preprocessData <- function(
    counts,
    margin,
    resolutions,
    assignmentFunction,
    coordinates = NULL
) {
    spotNames <- c("spot", "sample")
    geneNames <- c("gene", "feature")
    c(margin, otherMargin) %<-% {
        if (margin %in% spotNames) list("spot", "gene")
        else if (margin %in% geneNames) list("gene", "spot")
        else stop(sprintf(
            "invalid margin '%s' (must be one of: %s)",
            margin,
            paste(c(spotNames, geneNames), collapse = ", ")
        ))
    }

    spots <- colnames(counts)
    if (!is.null(coordinates)) {
        spots <- intersect(spots, rownames(coordinates))
        counts <- counts[, spots]
    } else {
        c(xcoord, ycoord) %<-% {
            strsplit(spots, "x") %>%
                transpose %>%
                map(as.numeric)
        }
        coordinates <- as.data.frame(cbind(x = xcoord, y = ycoord))
        rownames(coordinates) <- spots
    }

    assignments <-
        resolutions %>%
        map(function(r) {
            message(sprintf("Clustering resolution %s", deparse(r)))
            assignmentFunction(
                r,
                if (margin == "spot") t(counts)
                else {
                    log(as.matrix(counts) + 1) %>%
                        prop.table(margin = 2) %>%
                        apply(1, .zscore) %>%
                        t()
                }
            )
        }) %>%
        setNames(resolutions) %>%
        .tidyAssignments() %>%
        rename(!! sym(margin) := .data$unit)

    longCounts <-
        counts %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        gather("spot", "count", -.data$gene)

    clusterMeans <-
        assignments %>%
        inner_join(longCounts, by = margin) %>%
        group_by(
            .data$name,
            .data$resolution,
            .data$cluster,
            !! sym(otherMargin)
        ) %>%
        summarize(mean = mean(.data$count)) %>%
        ungroup()

    colors <-
        clusterMeans %>%
        select(
            .data$name,
            .data$mean,
            !! sym(otherMargin)
        ) %>%
        spread(.data$name, .data$mean) %>%
        as.data.frame() %>%
        column_to_rownames(otherMargin) %>%
        .computeClusterColors()

    message("Scoring spot-cluster affinity")
    scores <-
        if (margin == "spot") {
            countsAndMeans <-
                longCounts %>%
                inner_join(clusterMeans, by = "gene")
            distances <- as.data.table(countsAndMeans)[,
                .(distance = sqrt(mean((count - mean) ^ 2))),
                by = .(resolution, cluster, spot, name)
            ]
            distances %>%
                group_by(.data$resolution, .data$spot) %>%
                mutate(
                    score = .likeness(.data$distance / sum(.data$distance),
                    c = 40.
                )) %>%
                ungroup() %>%
                select(-.data$distance)
        } else {
            normalizedCounts <-
                longCounts %>%
                mutate(count = log(.data$count + 1)) %>%
                group_by(.data$spot) %>%
                mutate(count = .data$count / sum(.data$count)) %>%
                group_by(.data$gene) %>%
                mutate(count = .data$count / sum(.data$count)) %>%
                ungroup()
            assignments %>%
                inner_join(normalizedCounts, by = "gene") %>%
                group_by(
                    .data$resolution,
                    .data$spot,
                    .data$cluster,
                    .data$name
                ) %>%
                summarize(score = mean(.data$count)) %>%
                ungroup()
        }

    normalizedScores <-
        scores %>%
        group_by(.data$resolution, .data$spot) %>%
        mutate(score = .data$score / max(.data$score)) %>%
        ungroup()

    list(
        assignments = assignments %>% rename(unit = !! sym(margin)),
        counts = longCounts,
        means = clusterMeans,
        scores = normalizedScores,
        colors = colors,
        coordinates = coordinates,
        featureName = otherMargin
    )
}


#' SVG barplot
#'
#' @param xs named vector with observations
#' @return \code{\link{character}} SVG barplot
#' @keywords internal
.SVGBarplot <- function(xs) {
    invoke(
        paste,
        sprintf(
            paste0(
                "<svg width=\"20em\" height=\"1.5em\">",
                paste0(
                    "<rect width=\"%f%%\" height=\"1.5em\" ",
                    "style=\"fill:rgb(125,125,125)\"></rect>"
                ),
                paste0(
                    "<text x=\"2%%\" y=\"50%%\" fill=\"black\"",
                    "dominant-baseline=\"central\">%s</text>"
                ),
                paste0(
                    "<text x=\"%f%%\" y=\"50%%\" fill=\"white\"",
                    "dominant-baseline=\"central\" >%.2f</text>"
                ),
                "</svg>"
            ),
            70 * xs / max(xs),
            names(xs),
            70 * xs / max(xs) + 2,
            xs
        ),
        sep="\n"
    )
}


#' Array pie plot
#'
#' @param scores \code{\link[base]{data.frame}} with cluster scores for each
#' spot containing the columns `"spot"`, `"name"`, and `"score"`.
#' @param coordinates \code{\link[base]{data.frame}} with `rownames` matching
#' those in `scores` and columns `"x"` and `"y"` specifying the plotting
#' position of each observation.
#' @param image a \code{\link[grid]{grid.grob}} to use as background to the
#' plots.
#' @param scoreMultiplier log multiplication factor applied to the score vector.
#' @param spotScale pie chart size.
#' @param spotOpacity pie chart opacity.
#' @return \code{\link[ggplot2]{ggplot}} object of the pie plot.
#' @keywords internal
.arrayPlot <- function(
    scores,
    coordinates,
    counts = NULL,
    image = NULL,
    scoreMultiplier = 1.0,
    spotScale = 1,
    spotOpacity = 1,
    numTopGenes = 5
) {
    spots <- intersect(scores$spot, rownames(coordinates))

    r <- spotScale * min(dist(coordinates[spots, ])) / 2

    c(ymin, ymax) %<-% range(coordinates$y)
    c(xmin, xmax) %<-% range(coordinates$x)
    c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~. - 3 * r) }
    c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~. + 3 * r) }

    if (!is.null(image)) {
        ymin <- max(ymin, 1)
        ymax <- min(ymax, nrow(image$raster))
        xmin <- max(xmin, 1)
        xmax <- min(xmax, ncol(image$raster))

        image$raster <- image$raster[ymin:ymax, xmin:xmax]
        annotation <- ggplot2::annotation_custom(image, -Inf, Inf, -Inf, Inf)
    } else {
        annotation <- NULL
    }

    coordinates$y <- ymax - coordinates$y + ymin

    df <-
        coordinates %>%
        rownames_to_column("spot") %>%
        inner_join(scores, by="spot") %>%
        mutate(score = .data$score ^ scoreMultiplier) %>%
        mutate(tooltip = .data$spot)

    if (!is.null(counts)) {
        topGenes <-
            counts %>%
            group_by(.data$spot) %>%
            mutate(rank = rank(-.data$count, ties.method = "first")) %>%
            filter(.data$rank <= numTopGenes) %>%
            arrange(-.data$count) %>%
            summarize(topGenes = paste(
                .SVGBarplot(setNames(.data$count, .data$gene))
            ))
        df <-
            df %>%
            inner_join(topGenes, by = "spot") %>%
            mutate(tooltip = paste(sep = "<br />",
                .data$tooltip,
                .data$topGenes
            )) %>%
            select(-.data$topGenes)
    }

    ggplot() +
        annotation +
        geom_scatterpie_interactive(
            mapping = ggplot2::aes_string(
                x0 = "x", y0 = "y", r = "r", amount = "score", fill = "name",
                tooltip = "tooltip"
            ),
            data = df,
            alpha = spotOpacity,
            n = 64
        ) +
        coord_fixed() +
        scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
        scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
        theme_minimal() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()
        )
}


#' Cluster graph
#'
#' @param assignments \code{\link[base]{data.frame}} with columns `"name"`,
#' `"resolution"`, and `"cluster"`.
#' @param clusterMeans \code{\link[base]{data.frame}} with columns `"name"`,
#' `"resolution"`, `"cluster"`, `featureName`, and `"mean"`.
#' @param featureName \code{\link[base]{character}} with the name of the
#' clustered feature.
#' @param transitionProportions how to compute the transition proportions.
#' Possible values are:
#' - `"From"`: based on the total number of assignments in the lower-resolution
#' cluster
#' - `"To"`: based on the total number of assignments in the higher-resolution
#' cluster
#' @param transitionLabels \code{\link[base]{logical}} specifying whether to
#' show edge labels.
#' @param transitionThreshold hide edges with transition proportions below this
#' threshold.
#' @param numTopFeatures \code{\link[base]{integer}} specifying the number of
#' features to show in the hover tooltips.
#' @return \code{\link[ggplot2]{ggplot}} object of the cluster graph.
#' @keywords internal
.clusterGraph <- function(
    assignments,
    clusterMeans,
    featureName,
    transitionProportions = "To",
    transitionLabels = FALSE,
    transitionThreshold = 0.0,
    numTopFeatures = 10
) {
    transitionSym <-
        if (transitionProportions == "To") "toNode"
        else if (transitionProportions == "From") "node"
        else stop(sprintf(
            "Invalid value `transitionProportions`: %s",
            str(transitionProportions)
        ))

    data <-
        assignments %>%
        mutate(resolution = as.numeric(.data$resolution)) %>%
        rename(node = .data$name)

    graph <- igraph::graph_from_data_frame(
        d = data %>%
            mutate(toResolution = .data$resolution + 1) %>%
            (function(x) inner_join(
                x,
                x %>%
                    select(
                        everything(),
                        toCluster = .data$cluster,
                        toNode = .data$node
                    ),
                by = c("unit", "toResolution" = "resolution")
            )) %>%
            group_by(
                .data$node,
                .data$toNode,
                .data$cluster,
                .data$toCluster
            ) %>%
            summarize(transCount = n()) %>%
            group_by(!! sym(transitionSym)) %>%
            mutate(transProp = .data$transCount / sum(.data$transCount)) %>%
            ungroup() %>%

            group_by(.data$toNode) %>%
            filter(
                .data$transProp == max(.data$transProp)
                | .data$transProp > transitionThreshold
            ) %>%
            ungroup() %>%
            # ^ filter edges with transition proportions (weights) below
            #   threshold but always keep the incident edge with the highest
            #   weight (since the graph would become disconnected if that edge
            #   also were removed)

            select(.data$node, .data$toNode, everything()),
        vertices = data %>%
            group_by(.data$node, .data$resolution, .data$cluster) %>%
            summarize(size = n())
    )

    vertices <- cbind(
        igraph::layout_as_tree(graph, flip.y = FALSE) %>%
            `colnames<-`(c("y", "x")),
        igraph::get.vertex.attribute(graph) %>%
            as.data.frame(stringsAsFactors = FALSE)
    )

    edges <- c(
        igraph::get.edgelist(graph) %>%
            array_branch(2) %>%
            `names<-`(c("from", "to")),
        igraph::get.edge.attribute(graph)
    ) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        inner_join(
            vertices %>%
                select(.data$name, .data$x, .data$y),
            by = c("from" = "name")
        ) %>%
        inner_join(
            vertices %>%
                select(.data$name, xend = .data$x, yend = .data$y),
            by = c("to" = "name")
        )

    resolutionLabels <-
        vertices %>%
        select(.data$resolution, .data$x, .data$y) %>%
        filter(.data$resolution != 1) %>%
        mutate(ymin = min(.data$y), ymax = max(.data$y)) %>%
        group_by(.data$resolution) %>%
        summarize(
            x = mean(.data$x),
            y = first(.data$ymax) +
                0.1 * (first(.data$ymax) - first(.data$ymin))
        ) %>%
        mutate(
            label = as.character(
                levels(assignments$resolution)[.data$resolution])
        )

    tooltips <-
        clusterMeans %>%
        mutate(name = as.character(.data$name)) %>%
        group_by(.data$name) %>%
        mutate(rank = rank(-.data$mean, ties.method = "first")) %>%
        filter(.data$rank <= numTopFeatures) %>%
        arrange(-.data$mean) %>%
        summarize(tooltip = .SVGBarplot(setNames(
            mean,
            nm = !! sym(featureName)
        )))
    vertices <-
        vertices %>%
        inner_join(tooltips, by = "name") %>%
        mutate(tooltip = paste(sep = "<br />",
            toTitleCase(.data$name),
            sprintf("Size: %d", .data$size),
            .data$tooltip
        ))

    ggplot() +
        geom_segment(
            aes_string(
                "x", "y",
                xend = "xend", yend = "yend",
                alpha = "transProp"
            ),
            col = "black",
            data = edges
        ) +
        geom_point_interactive(
            aes_(
                ~x, ~y,
                color = ~name,
                size = ~size,
                tooltip = ~tooltip
            ),
            data = vertices %>% filter(.data$resolution != 1)
        ) +
        {
            if (isTRUE(transitionLabels))
                ggrepel::geom_label_repel(
                    aes_(
                        x = ~(x + xend) / 2,
                        y = ~(y + yend) / 2,
                        color =
                            if (transitionProportions == "To") ~as.factor(to)
                            else ~as.factor(from),
                        label = ~round(transProp, 2)
                    ),
                    data = edges,
                    show.legend = FALSE
                )
            else NULL
        } +
        ggplot2::geom_text(
            aes_string("x", "y", label = "label"),
            data = resolutionLabels
        ) +
        labs(alpha = "Proportion", color = "Cluster") +
        scale_size(guide = "none", range = c(2, 7)) +
        scale_x_continuous(expand = c(0.1, 0.1)) +
        guides(alpha = FALSE, color = FALSE) +
        theme_bw() +
        theme(
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
        )
}


#' SpatialCPie server
#'
#' @param assignments \code{\link[base]{data.frame}} with cluster assignments
#' containing the columns `"unit"` (name of the observational unit; either a
#' gene name or a spot name), `"resolution"`, `"cluster"`, and `"name"` (a
#' unique identifier of the (resolution, cluster) pair).
#' @param clusterMeans \code{\link[base]{data.frame}} with columns `"name"`,
#' `"resolution"`, `"cluster"`, `featureName`, and `"mean"`.
#' @param scores \code{\link[base]{data.frame}} with cluster scores for each
#' spot in each resolution containing the columns `"spot"`, `"resolution"`,
#' `"cluster"`, `"name"`, and `"score"`.
#' @param colors vector of colors for each cluster. Names should match the
#' `"name"` columns of the `assignments` and `scores`.
#' @param image background image for the array plots, passed to
#' \code{\link[grid]{grid.raster}}.
#' @param coordinates \code{\link[base]{data.frame}} with `rownames` matching
#' the \code{\link[base]{names}} in `scores` and columns `"x"` and `"y"`
#' specifying the plotting position of each observation.
#' @param featureName \code{\link[base]{character}} with the name of the
#' clustered feature.
#' @return server function, to be passed to \code{\link[shiny]{shinyApp}}.
#' @keywords internal
.makeServer <- function(
    assignments,
    clusterMeans,
    counts,
    scores,
    colors,
    image,
    coordinates,
    featureName
) {
    resolutions <-
        levels(assignments$resolution) %>%
        keep(~. != "root")

    function(input, output, session) {
        if (is.null(image)) {
            hideElement("showImage")
            hideElement("spotOpacity")
        }

        ###
        ## INPUTS
        edgeProportions <- reactive({ input$edgeProportions })
        edgeThreshold   <- reactive({ input$edgeThreshold   }) %>% debounce(1000)
        edgeLabels      <- reactive({ input$edgeLabels      })
        scoreMultiplier <- reactive({ input$scoreMultiplier }) %>% debounce(1000)
        showImage       <- reactive({ input$showImage       })
        spotOpacity     <- reactive({ input$spotOpacity     }) %>% debounce(1000)
        spotSize        <- reactive({ input$spotSize        }) %>% debounce(1000)

        ###
        ## CLUSTER GRAPH
        clusterGraph <- reactive({
            p <- .clusterGraph(
                assignments,
                clusterMeans,
                transitionProportions = edgeProportions(),
                transitionLabels = edgeLabels(),
                transitionThreshold = edgeThreshold(),
                featureName = featureName
            ) +
                scale_color_manual(values = colors)
        })

        output$clusterGraph <- ggiraph::renderGirafe({
            plot <- ggiraph::girafe_options(
                x = ggiraph::girafe(ggobj = clusterGraph()),
                ggiraph::opts_toolbar(saveaspng = FALSE)
            )
            plot
        })

        ###
        ## ARRAY PLOT
        arrayName <- function(r) sprintf("array%s", sha1(r))

        for (r in resolutions) {
            shiny::insertUI("#array", "beforeEnd",
                shiny::div(class = "array", "data-resolution" = r,
                    ggiraph::girafeOutput(arrayName(r)) %>%
                    shinycssloaders::withSpinner()
                ),
                immediate = TRUE
            )
            ## We evaluate the below block in a new frame (with anonymous
            ## function call) in order to protect the value of `r`, which
            ## will have changed when the reactive expressions are
            ## evaluated.
            (function() {
                r_ <- r
                scores_ <-
                    scores %>%
                    filter(.data$resolution == r_)
                assign(envir = parent.frame(), arrayName(r_), reactive(
                    .arrayPlot(
                        scores = scores_ %>%
                            select(.data$spot, .data$name, .data$score),
                        coordinates = coordinates,
                        counts = counts,
                        image =
                            if (!is.null(image) && !is.null(coordinates) &&
                                    showImage())
                                grid::rasterGrob(
                                    image,
                                    width = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    interpolate = TRUE
                                )
                            else NULL,
                        scoreMultiplier = 2 ** scoreMultiplier(),
                        spotScale = spotSize() / 5,
                        spotOpacity = spotOpacity()
                    ) +
                        guides(fill = guide_legend(title = "Cluster")) +
                        scale_fill_manual(
                            values = colors,
                            labels = unique(scores_$cluster)
                        )
                ))
                output[[arrayName(r_)]] <- ggiraph::renderGirafe(
                    {
                        ggiraph::girafe_options(
                            x = ggiraph::girafe(
                                ggobj = eval(call(arrayName(r_))),
                                xml_reader_options = list(options = "HUGE")),
                            ggiraph::opts_toolbar(saveaspng = FALSE),
                            ggiraph::opts_zoom(max = 5)
                        )
                    }
                )
            })()
        }

        ###
        ## EXPORT
        outputs <- reactive({
            list(
                clusters = assignments %>% select(-.data$name),
                clusterGraph = clusterGraph(),
                arrayPlot = lapply(
                    setNames(nm = resolutions),
                    function(x) eval(call(arrayName(x)))
                )
            )
        })
        observeEvent(input$done, shiny::stopApp(returnValue = outputs()))
    }
}


#' SpatialCPie UI
#'
#' @return SpatialCPie UI, to be passed to \code{\link[shiny]{shinyApp}}.
#' @keywords internal
.makeUI <- function() {
    shiny::htmlTemplate(system.file(
        "www", "default", "index.html", package = "SpatialCPie"))
}


#' SpatialCPie App
#'
#' @param image background image.
#' @param ... arguments passed to \code{\link{.preprocessData}}.
#' @return SpatialCPie \code{\link[shiny]{shinyApp}} object.
#' @keywords internal
.makeApp <- function(image, ...) {
    data <- .preprocessData(...)
    shiny::shinyApp(
        ui = .makeUI(),
        server = .makeServer(
            assignments = data$assignments,
            clusterMeans = data$means,
            counts = data$counts,
            scores = data$scores,
            colors = data$colors,
            image = image,
            coordinates = data$coordinates,
            featureName = data$featureName
        )
    )
}


#' Run SpatialCPie
#'
#' Runs the SpatialCPie gadget.
#' @param counts gene count matrix or a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object
#' containing count values.
#' @param image image to be used as background to the plot.
#' @param spotCoordinates \code{\link[base]{data.frame}} with pixel coordinates.
#' The rows should correspond to the columns (spatial areas) in the count file.
#' @param margin which margin to cluster.
#' @param resolutions \code{\link[base]{numeric}} vector specifying the
#' clustering resolutions.
#' @param assignmentFunction function to compute cluster assignments.
#' @param view \code{\link[shiny]{viewer}} object.
#' @return a list with the following items:
#' - `"clusters"`: Cluster assignments (may differ from `assignments`)
#' - `"clusterGraph"`: The cluster tree ggplot object
#' - `"arrayPlot"`: The pie plot ggplot objects
#' @export
#' @examples
#' if (interactive()) {
#'     options(device.ask.default = FALSE)
#'
#'     ## Set up coordinate system
#'     coordinates <- as.matrix(expand.grid(1:10, 1:10))
#'
#'     ## Generate data set with three distinct genes generated by three
#'     ## distinct cell types
#'     profiles <- diag(rep(1, 3)) + runif(9)
#'     centers <- cbind(c(5, 2), c(2, 8), c(8, 2))
#'     mixes <- apply(coordinates, 1, function(x) {
#'         x <- exp(-colSums((centers - x) ^ 2) / 50)
#'         x / sum(x)
#'     })
#'     means <- 100 * profiles %*% mixes
#'     counts <- matrix(rpois(prod(dim(means)), means), nrow = nrow(profiles))
#'     colnames(counts) <- apply(
#'         coordinates,
#'         1,
#'         function(x) do.call(paste, c(as.list(x), list(sep = "x")))
#'     )
#'     rownames(counts) <- paste("gene", 1:nrow(counts))
#'
#'     ## Run SpatialCPie
#'     runCPie(counts)
#' }
runCPie <- function(
    counts,
    image = NULL,
    spotCoordinates = NULL,
    margin = "spot",
    resolutions = 2:4,
    assignmentFunction = function(k, x) kmeans(x, centers = k)$cluster,
    view = NULL
) {
    if (is(counts, "SummarizedExperiment")) {
        counts <- as.data.frame(SummarizedExperiment::assay(counts))
    }
    shiny::runGadget(
        app = .makeApp(
            image = image,
            counts = counts,
            coordinates = spotCoordinates,
            margin = margin,
            resolutions = resolutions,
            assignmentFunction = assignmentFunction
        ),
        viewer = view %||% shiny::paneViewer()
    )
}
