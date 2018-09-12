#' @importFrom dplyr
#' filter group_by inner_join mutate n select summarize ungroup
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggplot2
#' aes coord_fixed element_blank geom_segment ggplot ggtitle guides guide_legend
#' labs
#' theme theme_bw
#' scale_color_manual scale_fill_manual scale_size
#' scale_x_continuous scale_y_continuous
#' @importFrom grid unit
#' @importFrom purrr
#' %>% %||% array_branch lift invoke keep map map_dbl map_int partial reduce
#' transpose
#' @importFrom readr read_file write_file
#' @importFrom shiny debounce observeEvent reactive
#' @importFrom stats dist setNames
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather spread
#' @importFrom tidyselect everything quo
#' @importFrom utils head tail
#' @importFrom utils str
#' @importFrom zeallot %<-%
"_PACKAGE"


## Pre-declare all NSE variables as global in order to appease R CMD check
## (ref: https://stackoverflow.com/a/12429344)
globalVariables(c(
    ".",
    "accAssignments",
    "accReassignments",
    "assignment",
    "cluster",
    "gene",
    "label",
    "name",
    "node",
    "resolution",
    "size",
    "spot",
    "toCluster",
    "toNode",
    "transCount",
    "transProp",
    "x",
    "xcoord",
    "xend",
    "y",
    "yend",
    "ycoord",
    NULL
))


#' Likeness score
#'
#' @param d (n, K) distance matrix
#' @param c log multiplier
#' @return (n, K) scoring matrix
#' @keywords internal
.likeness <- function(
    d,
    c = 1.0
) {
    score <- exp(-c * d)
    score / rowSums(score)
}


#' Pair-wise distances
#'
#' @param a (n, d) matrix of points
#' @param b (m, d) matrix of points
#' @param d distance function (with signature \[Num\] -> \[\[Num\]\] -> \[Num\])
#' @return (m, n) matrix of distances
#' @keywords internal
.pairwiseDistance <- function(
    a,
    b,
    d = function(x, ys) sqrt(colSums((t(ys) - x) ^ 2))
) {
    apply(a, 1, partial(d, ys = b))
}


#' Maximize overlap
#'
#' @param xs list of lists of labels
#' @return `xs`, relabeled so as to maximize the overlap between labels in
#' consecutive label lists
#' @keywords internal
.maximizeOverlap <- function(
    xs
) {
    reassignments <- list(
        unname(head(xs, -1)),
        unname(tail(xs, -1))
    ) %>%
        transpose %>%
        map(lift(function(prev, cur) {
            setNames(nm = sort(unique(prev))) %>%
                map(function(x)
                    setNames(nm = sort(unique(cur))) %>%
                        map_dbl(function(y) sum(`*`(prev == x, cur == y)))
                ) %>%
                invoke(rbind, .) %>%
                (function(x) {
                    nfrom <- nrow(x)
                    nto <- ncol(x)
                    n <- max(nfrom, nto)
                    x_ <- x
                    x_ <- do.call(
                        rbind,
                        c(list(x_), rep(list(rep(0, n)), n - nfrom))
                    )
                    x_ <- do.call(
                        cbind,
                        c(list(x_), rep(list(rep(0, n)), n - nto))
                    )
                    lpSolve::lp.assign(-x_)$solution[
                        seq_len(nfrom), seq_len(nto), drop = FALSE
                    ] %>%
                        `colnames<-`(colnames(x)) %>%
                        `rownames<-`(rownames(x))
                })
        }))

    c(
        list(xs[[1]]),
        list(tail(xs, -1), reassignments) %>%
            transpose %>%
            reduce(
                function(acc, x) {
                    c(accReassignments, accAssignments) %<-% acc
                    c(assignment, reassignment) %<-% x
                    rownames(reassignment) <- rownames(reassignment) %>%
                        map(~
                            if (. %in% accReassignments) accReassignments[.]
                            else .
                        )
                    reassignmentMap <- apply(
                        reassignment,
                        2,
                        function(x) rownames(reassignment)[which.max(x)]
                    )
                    reassignmentMap[!apply(reassignment, 2, max)] <-
                        setdiff(colnames(reassignment), rownames(reassignment))
                    list(
                        reassignmentMap,
                        c(
                            accAssignments,
                            list(vapply(
                                assignment,
                                function(x) as.numeric(reassignmentMap[x]),
                                numeric(1)
                            ))
                        )
                    )
                },
                .init = list(character(), list())
            ) %>%
            `[[`(2)
    ) %>%
        setNames(names(xs))
}


#' Array pie plot
#'
#' @param scores (n, K) scoring matrix
#' @param coordinates `data.frame` with `rownames` matching those in `scores`
#' and columns `x` and `y` specifying the plotting position of each observation
#' @param image a `grob` to use as background to the plots
#' @param spotScale pie chart size
#' @param spotOpacity pie chart opacity
#' @return `ggplot` object of the pie plot
#' @keywords internal
.arrayPlot <- function(
    scores,
    coordinates,
    image = NULL,
    spotScale = 1,
    spotOpacity = 1
) {
    spots <- intersect(rownames(scores), rownames(coordinates))

    r <- spotScale * min(dist(coordinates[spots, ])) / 2

    c(ymin, ymax) %<-% range(coordinates$y)
    c(xmin, xmax) %<-% range(coordinates$x)
    c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~. - 3 * r) }
    c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~. + 3 * r) }

    if (!is.null(image)) {
        ymin <- max(ymin, 0)
        ymax <- min(ymax, nrow(image$raster))
        xmin <- max(xmin, 0)
        xmax <- min(xmax, ncol(image$raster))

        image$raster <- image$raster[ymin:ymax, xmin:xmax]
        annotation <- ggplot2::annotation_custom(image, -Inf, Inf, -Inf, Inf)
    } else {
        annotation <- NULL
    }

    coordinates$y <- ymax - coordinates$y + ymin

    ggplot() +
        annotation +
        scatterpie::geom_scatterpie(
            mapping = aes(x = x, y = y, r = r),
            data = cbind(scores[spots, ], coordinates[spots, ]),
            cols = colnames(scores),
            alpha = spotOpacity
        ) +
        coord_fixed() +
        scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
        scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
        guides(fill = guide_legend(title = "Cluster")) +
        theme_bw() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()
        )
}


#' Cluster tree
#'
#' @param (n, R) assignment matrix, where R is the number of resolutions
#' @param transitionProportions how to compute the transition proportions.
#' Possible values are:
#' - `"From"`: based on the total number of assignments in the lower-resolution
#' cluster
#' - `"To"`: based on the total number of assignments in the higher-resolution
#' cluster
#' @param transitionLabels show edge labels
#' @param transitionThreshold hide edges with transition proportions below this
#' threshold
#' @return `ggplot` object of the cluster tree
#' @keywords internal
.clusterTree <- function(
    assignments,
    transitionProportions = "To",
    transitionLabels = FALSE,
    transitionThreshold = 0.0
) {
    if (!transitionProportions %in% c("From", "To")) {
        stop(sprintf(
            "Invalid value `transitionProportions`: %s",
            str(transitionProportions)
        ))
    }

    data <- assignments %>% as.data.frame %>%
        rownames_to_column("spot") %>%
        gather(resolution, cluster, -spot, convert = TRUE) %>%
        mutate(node = sprintf("R%dC%d", resolution, cluster))

    graph <- igraph::graph_from_data_frame(
        d = data %>%
            mutate(toResolution = resolution + 1) %>%
            (function(x) inner_join(
                x,
                x %>% select(everything(), toCluster = cluster, toNode = node),
                by = c("spot", "toResolution" = "resolution")
            )) %>%
            group_by(node, toNode, cluster, toCluster) %>%
            summarize(transCount = n()) %>%
            (function(x) {
                groupingVar <-
                    if (transitionProportions == "To") quo(toNode)
                    else quo(node)
                group_by(x, !!groupingVar)
            }) %>%
            mutate(transProp = transCount / sum(transCount)) %>%
            select(node, toNode, everything()),
        vertices = data %>%
            group_by(node, resolution, cluster) %>%
            summarize(size = n())
    )

    vertices <- cbind(
        igraph::layout.reingold.tilford(graph) %>% `colnames<-`(c("x", "y")),
        igraph::get.vertex.attribute(graph) %>%
            as.data.frame(stringsAsFactors = FALSE)
    )

    edges <- c(
        igraph::get.edgelist(graph) %>% array_branch(2) %>%
            `names<-`(c("from", "to")),
        igraph::get.edge.attribute(graph)
    ) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        filter(transProp > transitionThreshold) %>%
        inner_join(
            vertices %>% select(name, x, y),
            by = c("from" = "name")
        ) %>%
        inner_join(
            vertices %>% select(name, xend = x, yend = y),
            by = c("to" = "name")
        )

    ggplot() +
        geom_segment(
            aes(
                x, y,
                xend = xend, yend = yend,
                alpha = transProp
            ),
            col = "black",
            data = edges
        ) +
        geom_point_interactive(
            aes(
                x, y,
                data_id = as.factor(resolution),
                color = as.factor(cluster),
                size = size,
                tooltip = sprintf("%s: %d", name, size)
            ),
            data = vertices %>% filter(name != "R1C1")
        ) +
        {
            if (isTRUE(transitionLabels))
                ggrepel::geom_label_repel(
                    aes(
                        x = (x + xend) / 2,
                        y = (y + yend) / 2,
                        color = as.factor(
                            if (transitionProportions == "To") toCluster
                            else cluster
                        ),
                        label = round(transProp, 2)
                    ),
                    data = edges,
                    show.legend = FALSE
                )
            else NULL
        } +
        labs(alpha = "Proportion", color = "Cluster") +
        scale_size(guide = "none", range = c(2, 7)) +
        theme_bw() +
        theme(
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()
        )
}


#' SpatialCPie server
#'
#' @param distances list of spot-cluster distance matrices, one for each
#' resolution
#' @param colors vector of colors for each cluster label
#' @param image background image for the array plots, passed to
#' `grid::rasterGrob()`
#' @param coordinates `data.frame` with `rownames` matching the `names` in
#' `scores` and columns `x` and `y` specifying the plotting position of each
#' observation
#' @return server function, to be passed to `shiny::shinyApp()`
#' @keywords internal
.makeServer <- function(
    distances,
    colors,
    image,
    coordinates
) {
    assignments <- distances %>%
        map(~apply(., 1, function(x) as.numeric(colnames(.))[which.min(x)])) %>%
        invoke(cbind, .)

    function(input, output, session) {
        ###
        ## INPUTS
        edgeProportions <- reactive({ input$edgeProportions })
        edgeThreshold   <- reactive({ input$edgeThreshold   }) %>% debounce(500)
        edgeLabels      <- reactive({ input$edgeLabels      })
        showImage       <- reactive({ input$showImage       })
        scoreMultiplier <- reactive({ input$simC            }) %>% debounce(500)
        spotOpacity     <- reactive({ input$spotOpacity     }) %>% debounce(500)
        spotSize        <- reactive({ input$spotSize        }) %>% debounce(500)

        setSelection <- function(value) {
            assign("selection", value, pos = parent.frame())
        }
        reactiveSelection <- reactive({ setSelection(input$tree_selected) })
        setSelection(tail(names(distances), -1))

        ###
        ## CLUSTER TREE
        treePlot <- reactive({
            .clusterTree(
                assignments,
                transitionProportions = edgeProportions(),
                transitionLabels = edgeLabels() == "Show",
                transitionThreshold = edgeThreshold()
            ) +
                scale_color_manual(values = colors)
        })

        output$tree <- ggiraph::renderggiraph({
            plot <- ggiraph::ggiraph(
                code = print(treePlot()),
                hover_css = paste(
                    "stroke:#888;",
                    "stroke-width:0.2em;",
                    "stroke-opacity:0.5;"
                ),
                selected_css = paste(
                    "stroke:#000;",
                    "stroke-width:0.2em;",
                    "stroke-opacity:0.5;"
                )
            )

            ## Copy selection from the previous tree
            ## Note: This could also have been done by sending the
            ## `"tree_set"` message on the `"onFlushed"` event. However,
            ## that would create a race condition between the message and
            ## the browser's loading of the ggiraph dependency file; if the
            ## latter is loaded last, the selection would be overwritten.
            ## Thus, we instead modify the dependency file directly to
            ## include the plot's initial selection.
            if (length(selection) > 0) {
                dependency <- paste(
                    plot$dependencies[[1]]$src$file,
                    plot$dependencies[[1]]$script,
                    sep = "/"
                )
                src <- read_file(dependency)

                ## Refresh selection annotations when plot is initialized
                src <- sub(
                    "(function init_prop_[^{]*\\{[^}]*)",
                    sprintf(
                        "\\1refresh_selected('%s', '%s', '%s');",
                        plot$x$sel_array_name,
                        plot$x$selected_class,
                        plot$x$uid
                    ),
                    src
                )

                ## Set initial selection
                selection_array <- sprintf(
                    "[%s]",
                    paste0("'", selection, "'") %>%
                        invoke(paste, sep = ",", .)
                )
                src <- paste0(
                    src,
                    sprintf("%s = %s;", plot$x$sel_array_name, selection_array),
                    sprintf("Shiny.setInputValue('tree_selected', %s);",
                        selection_array)
                )

                write_file(src, dependency)
            }

            ## Remove UI element for lasso selection etc.
            plot$x$ui_html <- ""

            plot
        })

        ###
        ## ARRAY PLOT
        for (d in tail(names(distances), -1)) {
            ## We evaluate the below block in a new frame (with anonymous
            ## function call) in order to protect the value of `d`, which
            ## will have changed when the reactive expressions are
            ## evaluated.
            (function() {
                d_ <- d
                infoName <- sprintf("array_info_%s", d_)
                plotName <- sprintf("array_%s", d_)
                assign(envir = parent.frame(), infoName, reactive(
                    .likeness(distances[[d_]], scoreMultiplier())
                ))
                assign(envir = parent.frame(), plotName, reactive(
                    .arrayPlot(
                        scores = eval(call(infoName)),
                        coordinates = coordinates,
                        image =
                            if (!is.null(image) && !is.null(coordinates) &&
                                    showImage() == "Show")
                                grid::rasterGrob(
                                    image,
                                    width = unit(1, "npc"),
                                    height = unit(1, "npc"),
                                    interpolate = TRUE
                                )
                            else NULL,
                        spotScale = spotSize() / 5,
                        spotOpacity = spotOpacity() / 100
                    ) +
                        scale_fill_manual(values = colors) +
                        ggtitle(sprintf("Resolution %s", d_))
                ))
                output[[paste0("plot", d_)]] <- shiny::renderPlot(
                    {
                        message(sprintf("Loading resolution \"%s\"...", d_))
                        eval(call(plotName))
                    },
                    width = 500, height = 500
                )
            })()
        }

        output$array <- shiny::renderUI({
            sort(as.numeric(input$tree_selected)) %>%
                map(~paste0("plot", .)) %>%
                keep(~. %in% names(shiny::outputOptions(output))) %>%
                map(~shiny::div(
                    style = "display: inline-block;",
                    shiny::div(
                        style = paste(
                            "position: relative;",
                            "width: 500px;",
                            "height: 500px;"
                        ),
                        list(
                            shiny::plotOutput(.),
                            shiny::div(
                                style = paste(
                                    "z-index: -1;",
                                    "position: absolute;",
                                    "top: 50%; left: 50%;"
                                ),
                                shiny::div(
                                    style = paste(
                                        "background: #eee;",
                                        "padding: 1em;",
                                        "position: relative;",
                                        "left: -50%;"
                                    ),
                                    "Loading..."
                                )
                            )
                        )
                    )
                )) %>%
                invoke(shiny::div, style = "text-align:center", .)
        })

        ###
        ## EXPORT
        outputs <- reactive({
            list(
                clusters = assignments[, selection],
                treePlot = treePlot(),
                piePlots = lapply(
                    setNames(nm = selection),
                    function(x) eval(call(sprintf("array_%s", x)))
                ),
                piePlotsInfo = lapply(
                    setNames(nm = selection),
                    function(x) eval(call(sprintf("array_info_%s", x)))
                )
            )
        })
        observeEvent(input$done, shiny::stopApp(returnValue = outputs()))
    }
}


#' SpatialCPie UI
#'
#' @param imageButton show image radio buttons
#' @return SpatialCPie UI, to be passed to `shiny::shinyApp()`
#' @keywords internal
.makeUI <- function(
    imageButton = FALSE
) {
    miniUI::miniPage(
        shiny::tags$head(shiny::tags$style(shiny::HTML(
            ".recalculating { position: relative; z-index: -2 }"))),
        miniUI::gadgetTitleBar("SpatialCPie"),
        miniUI::miniContentPanel(
            shiny::fillPage(
                shiny::sidebarLayout(
                    shiny::sidebarPanel(width = 3,
                        shiny::radioButtons(
                            "edgeLabels",
                            "Edge labels:", c("Show", "Hide")
                        ),
                        shiny::radioButtons(
                            "edgeProportions",
                            "Edge proportions:", c("To", "From")
                        ),
                        shiny::numericInput(
                            "edgeThreshold",
                            "Min proportion:",
                            max = 1, min = 0, value = 0.05, step = 0.01
                        )
                    ),
                    shiny::div(
                        style = "text-align: center",
                        shiny::mainPanel(ggiraph::ggiraphOutput("tree"))
                    )
                ),
                shiny::hr(),
                shiny::sidebarLayout(
                    shiny::sidebarPanel(width = 3,
                        if (isTRUE(imageButton))
                            shiny::radioButtons(
                                "showImage",
                                "HE image:", c("Show", "Hide")
                            )
                        else NULL,
                        shiny::numericInput(
                            "simC",
                            "Score multiplier:",
                            max = 10, min = 0.1, value = 1, step = 0.2
                        ),
                        shiny::numericInput(
                            "spotOpacity",
                            "Opacity:",
                            max = 100, min = 1, value = 100, step = 10
                        ),
                        shiny::numericInput(
                            "spotSize",
                            "Size:",
                            max = 10, min = 1, value = 5, step = 1
                        )
                    ),
                    shiny::mainPanel(shiny::uiOutput("array"))
                )
            )
        )
    )
}


#' SpatialCPie App
#'
#' @usage See `runCPie()`
#' @return SpatialCPie `shiny::shinyApp()` object
.makeApp <- function(
    counts,
    assignments,
    image = NULL,
    spotCoordinates = NULL
) {
    ## Compute coordinates and intersect spots across resolutions
    spots <- assignments %>% map(names) %>% reduce(intersect)

    if (!is.null(spotCoordinates) != 0) {
        spots <- intersect(spots, rownames(spotCoordinates))
        coordinates <- spotCoordinates
    } else {
        c(xcoord, ycoord) %<-% {
            strsplit(spots, "x") %>% transpose %>% map(as.numeric) }
        coordinates <- as.data.frame(cbind(x = xcoord, y = ycoord))
        rownames(coordinates) <- spots
    }

    assignments <- assignments %>% map(~.[spots])
    counts <- counts[, spots]

    ## Relabel the data, making sure that all labels in resolution r are in
    ## [1..r]
    assignments <- lapply(
        assignments,
        function(assignment) {
            oldLabels <- unique(assignment)
            newLabels <- setNames(seq_along(oldLabels), nm = oldLabels)
            setNames(
                newLabels[as.character(assignment)],
                nm = names(assignment)
            )
        }
    )

    ## Sort resolutions
    names(assignments) <- assignments %>% map_int(max)
    assignments <- assignments[
        as.character(sort(as.numeric(names(assignments))))]

    ## Add first resolution (corresponding to the root of the tree)
    if (names(assignments)[1] != "1") {
        assignments <- c(
            list("1" = setNames(rep(1, length(spots)), nm = spots)),
            assignments
        )
    }

    ## Relabel the data to maximize overlap between labels in consecutive
    ## resolutions
    assignments <- .maximizeOverlap(assignments)

    ## Compute sample-centroid distances
    distances <- assignments %>%
        map(function(assignment) {
            labels <- unique(assignment)
            centers <- labels %>%
                map(~rowMeans(
                    counts[, which(assignment == .), drop = FALSE]
                )) %>%
                invoke(cbind, .)
            colnames(centers) <- labels
            .pairwiseDistance(t(centers), t(counts))
        })

    ## Compute colors so that dissimilar clusters are far away in color space
    colors <- assignments %>% unname %>% invoke(c, .) %>%
        (function(x) data.frame(
            spot = names(x),
            cluster = x,
            stringsAsFactors = FALSE
        )) %>%
        inner_join(
            counts %>%
                as.data.frame(stringsAsFactors = FALSE) %>%
                tibble::rownames_to_column("gene") %>%
                gather(spot, counts, -gene),
            by = "spot"
        ) %>%
        group_by(cluster, gene) %>%
        summarize(mean = mean(counts)) %>%
        spread(cluster, mean) %>%
        select(-gene) %>%
        as.matrix %>% t %>%
        stats::prcomp(rank = 2, center = TRUE) %>% `$`("x") %>%
        (function(x) cbind(
            50,
            (2 * (x - min(x)) / (max(x) - min(x)) - 1) * 100
        )) %>%
        (colorspace::LAB) %>%
        (colorspace::hex)(fixup = TRUE)

    shiny::shinyApp(
        ui = .makeUI(!is.null(image)),
        server = .makeServer(
            distances = distances,
            colors = colors,
            image = image,
            coordinates = coordinates
        )
    )
}


#' Run SpatialCPie
#'
#' Runs the SpatialCPie gadget.
#' @param counts gene count matrix.
#' @param assignments list of cluster assignments for each resolution.
#' @param image image to be used as background to the plot.
#' @param spotCoordinates `data.frame` with pixel coordinates. The rows should
#' correspond to the columns (spatial areas) in the count file.
#' @param view `shiny::viewer` object
#' @return a list with the following items:
#' - `"clusters"`: Cluster assignments (may differ from `assignments`)
#' - `"treePlot"`: The cluster tree ggplot object
#' - `"piePots"`: The pie plot ggplot objects
#' - `"piePlotsInfo"`: Likeness scores between spots and clusters in each
#' resolution
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
#'     ## Perform clustering
#'     assignments <- lapply(
#'         2:5, function(x) kmeans(t(counts), centers = x)$cluster)
#'
#'     ## Run SpatialCPie
#'     runCPie(counts, assignments)
#' }
runCPie <- function(
    counts,
    assignments,
    image = NULL,
    spotCoordinates = NULL,
    view = NULL
) {
    shiny::runGadget(
        app = .makeApp(
            counts,
            assignments,
            image,
            spotCoordinates
        ),
        viewer = view %||% shiny::dialogViewer("SpatialCPie")
    )
}
