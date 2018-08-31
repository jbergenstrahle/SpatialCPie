#' Likeness score
#'
#' @keywords internal
.likeness <- function(d, c = 1.0) {
    score <- exp(-c * d)
    score / rowSums(score)
}

#' Pair-wise distances
#'
#' @keywords internal
#' @import purrr
.pairwiseDistance <- function(
    a, b, d = function(xs, y) sqrt(colSums((t(xs) - y) ^ 2))
) {
    apply(a, 1, partial(d, b))
}

#' Maximize overlap
#'
#' Given a list of of lists with labels, relabel the lists so as to maximize the
#' overlap between labels in consecutive lists.
#' @keywords internal
.maximizeOverlap <- function(xs) {
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
#' @keywords internal
#' @import ggplot2 purrr scatterpie zeallot
#' @importFrom stats dist
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
    c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~. - 2 * r) }
    c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~. + 2 * r) }

    if (!is.null(image)) {
        ymin <- max(ymin, 0)
        ymax <- min(ymax, nrow(image$raster))
        xmin <- max(xmin, 0)
        xmax <- min(xmax, ncol(image$raster))

        image$raster <- image$raster[ymin:ymax, xmin:xmax]
        annotation <- annotation_custom(image, -Inf, Inf, -Inf, Inf)
    } else {
        annotation <- NULL
    }

    coordinates$y <- ymax - coordinates$y + ymin

    ggplot() +
        annotation +
        geom_scatterpie(
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
#' @keywords internal
#' @import ggplot2 purrr ggiraph intergraph
#' @importFrom dplyr filter group_by inner_join mutate n select summarise
#' @importFrom dplyr ungroup
#' @importFrom igraph graph_from_data_frame layout.reingold.tilford
#' @importFrom ggrepel geom_label_repel
#' @importFrom tidyr gather
#' @importFrom tidyselect everything
#' @importFrom utils str tail
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

    transitions <- lapply(seq_len(ncol(assignments) - 1), function(i) {
        fromResolution <- colnames(assignments)[i]
        toResolution <- colnames(assignments)[i + 1]
        cbind(
            fromResolution = as.numeric(fromResolution),
            toResolution = as.numeric(toResolution),
            expand.grid(
                fromCluster = sort(unique(assignments[, fromResolution])),
                toCluster = sort(unique(assignments[, toResolution]))
            ) %>%
                transpose %>%
                map(lift(function(fromCluster, toCluster) {
                    isFrom <- assignments[, fromResolution] == fromCluster
                    isTo <- assignments[, toResolution] == toCluster
                    count <- sum(isFrom & isTo)
                    prop <- if (transitionProportions == "From")
                                count / sum(isFrom)
                            else count / sum(isTo)
                    c(
                        fromCluster = fromCluster,
                        toCluster = toCluster,
                        transCount = count,
                        transProp = prop
                    )
                })) %>%
                invoke(rbind, .)
        )
    }) %>%
        invoke(rbind, .) %>%
        as.data.frame %>%
        mutate(fromCluster = factor(fromCluster)) %>%
        mutate(toCluster = factor(toCluster))
    browser()

    graph <- graph_from_data_frame(
        transitions %>%
            filter(transCount > 0) %>%
            mutate(fromNode = paste0("R", fromResolution, "C", fromCluster)) %>%
            mutate(toNode = paste0("R", toResolution, "C", toCluster)) %>%
            select(fromNode, toNode, everything()),
        vertices = assignments %>%
            as.data.frame %>%
            gather(key = resolution, value = cluster) %>%
            group_by(resolution, cluster) %>%
            summarize(size = n()) %>%
            ungroup() %>%
            mutate(resolution = as.numeric(resolution)) %>%
            mutate(cluster = as.numeric(cluster)) %>%
            mutate(node = paste0("R", resolution, "C", cluster)) %>%
            select(node, everything())
    )

    vertices <- cbind(
        layout.reingold.tilford(graph) %>% `colnames<-`(c("x", "y")),
        igraph::get.vertex.attribute(graph) %>%
            as.data.frame(stringsAsFactors = FALSE)
    ) %>%
        mutate(label = sprintf(
            "(resolution %d; cluster %d)",
            resolution,
            cluster
        ))

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
                data_id = resolution,
                color = as.factor(cluster),
                size = size,
                tooltip = label
            ),
            data = vertices %>% filter(name != "R1C1")
        ) +
        {
            if (isTRUE(transitionLabels))
                geom_label_repel(
                    aes(
                        x = (x + xend) / 2,
                        y = (y + yend) / 2,
                        color = if (transitionProportions == "To") toCluster
                                else fromCluster,
                        label = round(transProp, 2)
                    ),
                    data = edges,
                    show.legend = FALSE
                )
            else NULL
        } +
        labs(alpha = "Proportion", color = "Cluster") +
        scale_size(guide = "none", range = c(2, 7))
}

#' SpatialCPie server
#'
#' @keywords internal
.makeServer <- function(
    counts,
    assignments,
    distances,
    colors,
    img,
    coords
) {
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

        ###
        ## CLUSTER TREE
        treePlot <- reactive({
            .clusterTree(
                do.call(cbind, assignments),
                transitionProportions = edgeProportions(),
                transitionLabels = edgeLabels() == "Show",
                transitionThreshold = edgeThreshold()
            ) +
                scale_color_manual(values = colors)
        })

        output$tree <- renderggiraph({
            plot <- ggiraph(
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
            if (length(input$tree_selected) > 0) {
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
                    paste0("'", input$tree_selected, "'") %>%
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
                        coordinates = coords,
                        image =
                            if (!is.null(img) && !is.null(coords) &&
                                    showImage() == "Show")
                                rasterGrob(
                                    img,
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
                output[[paste0("plot", d_)]] <- renderPlot(
                    {
                        message(sprintf("Loading resolution \"%s\"...", d_))
                        eval(call(plotName))
                    },
                    width = 400, height = 400
                )
            })()
        }

        output$array <- renderUI({
            sort(as.numeric(input$tree_selected)) %>%
                map(~paste0("plot", .)) %>%
                keep(~. %in% names(outputOptions(output))) %>%
                map(~div(
                    style = "display: inline-block;",
                    div(
                        style = paste(
                            "position: relative;",
                            "width: 400px;",
                            "height: 400px;"
                        ),
                        list(
                            plotOutput(.),
                            div(
                                style = paste(
                                    "z-index: -1;",
                                    "position: absolute;",
                                    "top: 50%; left: 50%;"
                                ),
                                div(
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
                invoke(div, style = "text-align:center", .)
        })

        ###
        ## EXPORT
        observeEvent(input$done, {
            stopApp(returnValue = list(
                clusters = assignments[input$tree_selected],
                tree = treePlot(),
                piePlots = lapply(
                    setNames(nm = input$tree_selected),
                    function(x) eval(call(sprintf("array_%s", x)))
                ),
                piePlotsInfo = lapply(
                    setNames(nm = input$tree_selected),
                    function(x) eval(call(sprintf("array_info_%s", x)))
                )
            ))
        })

        ###
        ## STARTUP
        session$onFlushed(function() {
            session$sendCustomMessage( "tree_set", tail(names(distances), -1))
        })
    }
}

#' SpatialCPie UI
#'
#' @keywords internal
.ui <- miniPage(
    tags$head(tags$style(HTML(
        ".recalculating { position: relative; z-index: -2 }"))),
    gadgetTitleBar("SpatialCPie"),
    miniContentPanel(
        fillPage(
            sidebarLayout(
                sidebarPanel(width = 3,
                    radioButtons(
                        "edgeLabels", "Edge labels:", c("Show", "Hide")
                    ),
                    radioButtons(
                        "edgeProportions", "Edge proportions:", c("To", "From")
                    ),
                    numericInput(
                        "edgeThreshold", "Min proportion:",
                        max = 1, min = 0, value = 0.05, step = 0.01
                    )
                ),
                div(
                    style = "text-align: center",
                    mainPanel(ggiraphOutput("tree"))
                )
            ),
            hr(),
            sidebarLayout(
                sidebarPanel(width = 3,
                    if (!is.null(img))
                        radioButtons(
                            "showImage", "HE image:", c("Show", "Hide")
                        )
                    else NULL,
                    numericInput(
                        "simC", "Score multiplier:",
                        max = 10, min = 0.1,
                        value = 1, step = 0.2
                    ),
                    numericInput(
                        "spotOpacity", "Opacity:",
                        max = 100, min = 1,
                        value = 100, step = 10
                    ),
                    numericInput(
                        "spotSize", "Size:",
                        max = 10, min = 1,
                        value = 5, step = 1
                    )
                ),
                mainPanel(uiOutput("array"))
            )
        )
    )
)

#' Run SpatialCPie
#'
#' Runs the SpatialCPie gadget
#' @param counts gene count matrix.
#' @param assignments list of cluster assignments for each resolution.
#' @param img image to be used as background to the plot (optional).\cr
#' Note: For the tissue image to be correctly aligned to the spatial areas, the
#' pixel.coords argument also needs to be provided.
#' @param pixel.coords `data.frame` with pixel coordinates (optional).
#' The rows should correspond to the columns (spatial areas) in the count file.
#' @param view Shiny gadgets viewer options.
#' Available options: "dialog" (default), "browser", "pane".
#' @return FIXME
#' @keywords arrayplot arraypieplot clustertree
#' @export
#' @import shiny miniUI ggplot2 grid zeallot grid purrr readr
#' @importFrom colorspace LAB hex
#' @importFrom dplyr summarize
#' @importFrom stats prcomp setNames
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather spread
#' @importFrom utils head tail
#' @examples
#' options(device.ask.default = FALSE)
#'
#' ## Set up coordinate system
#' coords <- as.matrix(expand.grid(1:10, 1:10))
#'
#' ## Generate data set with three distinct genes
#' profiles <- diag(rep(1, 3)) + runif(9)
#' centers <- cbind(c(5, 2), c(2, 8), c(8, 2))
#' mixes <- apply(coords, 1, function(x) {
#'     x <- exp(-colSums((centers-x)^2)/50)
#'     x / sum(x)
#' })
#' means <- 100 * profiles %*% mixes
#' counts <- matrix(rpois(prod(dim(means)), means), nrow=nrow(profiles))
#' colnames(counts) <- apply(
#'     coords, 1, function(x) do.call(paste, c(as.list(x), list(sep="x"))))
#' rownames(counts) <- paste("gene", 1:nrow(counts))
#'
#' ## Perform clustering
#' assignments <- lapply(
#'     2:5, function(x) kmeans(t(counts), centers=x)$cluster)
#'
#' ## Run SpatialCPie
#' runCPie(counts, assignments)
runCPie <- function(
    counts,
    assignments,
    img = NULL,
    view = NULL,
    pixel.coords = NULL
) {
    ## Compute coordinates and intersect spots across resolutions
    spots <- assignments %>% map(names) %>% reduce(intersect)

    if (!is.null(pixel.coords) != 0) {
        spots <- intersect(spots, rownames(pixel.coords))
        coords <- pixel.coords
    } else {
        c(xcoord, ycoord) %<-% {
            strsplit(spots, "x") %>% transpose %>% map(as.numeric) }
        coords <- as.data.frame(cbind(x = xcoord, y = ycoord))
        rownames(coords) <- spots
    }

    assignments <- assignments %>% map(~.[spots])

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
                rownames_to_column("gene") %>%
                gather(spot, counts, -gene),
            by = "spot"
        ) %>%
        group_by(cluster, gene) %>%
        summarize(mean = mean(counts)) %>%
        spread(cluster, mean) %>%
        select(-gene) %>%
        as.matrix %>% t %>%
        prcomp(rank = 2, center = TRUE) %>% `$`("x") %>%
        (function(x) cbind(
            50,
            (2 * (x - min(x)) / (max(x) - min(x)) - 1) * 100
        )) %>%
        LAB %>% hex(fixup = TRUE)

    runGadget(
        viewer = view %||% dialogViewer(),
        server = .makeServer(
            assignments = assignments,
            distances = distances,
            colors = colors,
            img = img,
            coords = coords
        ),
        app = .ui
    )
}
