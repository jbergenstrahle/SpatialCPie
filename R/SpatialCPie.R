#' Likeness score
#'
#' Computes the normalized likeness score from a given distance matrix.
#' @param d distance matrix of dimension (n, k), where n is the number of data
#' points and k is the number of clusters.
#' @param c difference weighting constant. A higher value implies that higher
#' distances will be punished more. In the limit c -> ∞, the closest cluster
#' will be assigned a score of 1 and the other clusters a score of 0.
#' @return (n, k) matrix of likeness scores.
#' @keywords internal
likeness <- function(d, c = 1.0) {
  score <- exp(-c * d)
  score / rowSums(score)
}

#' Pair-wise distances
#'
#' Computes all pair-wise distances between two point clouds.
#' @param a FIXME
#' @param b FIXME
#' @param d FIXME
#' @return FIXME
#' @keywords internal
#' @import purrr
pairwiseDistance <- function(
  a, b, d = function(xs, y) sqrt(colSums((t(xs)-y)^2))
) {
  apply(a, 1, partial(d, b))
}

#' Array pie plot
#'
#' Creates the pie plot.
#' @param scores matrix of dimension (n, k), where n is the number of data
#' points and k is the number of categories that are displayed in the pie
#' charts.
#' @param coordinates data.frame of dimension (n, 2) encoding the plotting
#' positions of the n data points. `coordinates` should contain the columns 'x'
#' and 'y' and its `rownames` should correspond to those in `scores`.
#' @param image a grob to display as a background image to the plot.
#' @param spotScale size of the pie charts.
#' @param spotOpacity opacity of the pie charts.
#' @return ggplot object of the array pie plot.
#' @keywords internal
#' @import ggplot2 purrr scatterpie zeallot
#' @importFrom stats dist
arrayPlot <- function(scores, coordinates, image = NULL, spotScale = 1,
                      spotOpacity = 1) {
  spots <- intersect(rownames(scores), rownames(coordinates))

  r <- spotScale * min(dist(coordinates[spots, ])) / 2

  c(ymin, ymax) %<-% range(coordinates$y)
  c(xmin, xmax) %<-% range(coordinates$x)
  c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~.-2*r) }
  c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~.+2*r) }

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
#' Creates the clustering tree.
#' @param graph graph of the clustering tree, as received from `CPieClust`.
#' @param transitionProportions how to compute and display edge labels. Valid values
#' are:
#' * `"To"`: Proportions are computed based on the number of features in the
#' higher resolution cluster.
#' * `"From"`: Proportions are computed based on the number of features in the
#' lower resolution cluster.
#' @param transitionLabels Show transition proportion labels.
#' @return ggplot object of the clustering tree.
#' @keywords internal
#' @import ggplot2 purrr ggiraph ggnetwork intergraph
#' @importFrom dplyr filter group_by mutate n select summarise ungroup
#' @importFrom igraph graph_from_data_frame layout.reingold.tilford
#' @importFrom ggrepel geom_label_repel
#' @importFrom tidyr gather
#' @importFrom tidyselect everything
#' @importFrom utils tail
clusterTree <- function(
  assignments,
  transitionProportions = "To",
  transitionLabels = FALSE,
  transitionThreshold = 0.0
) {
  if (!transitionProportions %in% c("From", "To")) {
    stop(sprintf("Invalid value `transitionProportions`: %s",
                 str(transitionProportions)))
  }

  ## Compute edges
  # Loop over resolutions
  transitions <- lapply(seq_len(ncol(assignments) - 1), function(i) {
    # Extract two neighbouring assignments
    from.res <- colnames(assignments)[i]
    to.res <- colnames(assignments)[i + 1]

    # Get the cluster names
    from.clusters <- sort(unique(assignments[, from.res]))
    to.clusters <- sort(unique(assignments[, to.res]))

    # Get all possible combinations
    trans.df <- expand.grid(FromClust = from.clusters,
                            ToClust = to.clusters)

    # Loop over the possible transitions
    trans <- apply(trans.df, 1, function(x) {
      from.clust <- x[1]
      to.clust <- x[2]

      # Find the ST-spots from those clusters
      is.from <- assignments[, from.res] == from.clust
      is.to <- assignments[, to.res] == to.clust

      # Count them up
      trans.count <- sum(is.from & is.to)

      # Get the sizes of the two clusters
      from.size <- sum(is.from)
      to.size <- sum(is.to)

      # Get the proportions of ST-spots moving along this edge
      trans.prop <-
        if (transitionProportions == "From") trans.count / from.size
        else trans.count / to.size

      return(c(trans.count, trans.prop))
    })

    # Tidy up the results
    trans.df$FromRes <- as.numeric(gsub("res.", "", from.res))
    trans.df$ToRes <- as.numeric(gsub("res.", "", to.res))
    trans.df$TransCount <- trans[1, ]
    trans.df$TransProp <- trans[2, ]

    return(trans.df)
  })

  # Bind the results from the different resolutions together
  transitions <- do.call("rbind", transitions)

  # Tidy everything up
  levs <- levels(as.factor(sort(transitions$ToClust)))
  transitions <- transitions %>%
    mutate(FromClust = factor(FromClust, levels = levs))  %>%
    mutate(ToClust = factor(ToClust, levels = levs))

  ## Compute nodes
  nodes <- assignments %>%
    as.data.frame %>%
    gather(key = Res, value = Cluster) %>%
    group_by(Res, Cluster) %>%
    summarise(Size = n()) %>%
    ungroup() %>%
    mutate(Res = stringr::str_replace(Res, "res.", "")) %>%
    mutate(Res = as.numeric(Res), Cluster = as.numeric(Cluster)) %>%
    mutate(Node = paste0("R", Res, "C", Cluster)) %>%
    select(Node, everything())

  ## Construct graph
  graph <- transitions %>%
    filter(TransCount > 0) %>%
    # Rename the nodes
    mutate(FromNode = paste0("R", FromRes, "C", FromClust)) %>%
    mutate(ToNode = paste0("R", ToRes, "C", ToClust)) %>%
    # Reorder columns
    select(FromNode, ToNode, everything()) %>%
    # Build a graph using igraph
    graph_from_data_frame(vertices = nodes)

  vertices <- cbind(
    layout.reingold.tilford(graph) %>% `colnames<-`(c("x", "y")),
    igraph::get.vertex.attribute(graph) %>% as.data.frame
  )
  rownames(vertices) <- vertices$name
  vertices <- vertices[, colnames(vertices) != "name"]
  vertices$label <- sprintf(
    "(resolution %d; cluster %d)",
    vertices$Res,
    vertices$Cluster
  )

  edges <- c(
    igraph::get.edgelist(graph) %>% array_branch(2) %>%
      `names<-`(c("from", "to")),
    igraph::get.edge.attribute(graph)
  ) %>%
    transpose %>%
    discard(lift(
      function(TransProp, ...) TransProp < transitionThreshold
    )) %>%
    map(lift(function(from, to, ...) {
      cbind(vertices[from, c("x", "y")], vertices[to, c("x", "y")] %>%
        `colnames<-`(c("xend", "yend"))) %>%
        cbind(...)
    })) %>%
    invoke(rbind, .)

  # Drop root vertex; root should not be selectable
  vertices <- vertices[rownames(vertices) != "R1C1", ]

  ## Construct plot
  ggplot() +
    geom_segment(
      aes(
        x, y,
        xend = xend, yend = yend,
        alpha = TransProp
      ),
      col = "black",
      data = edges
    ) +
    geom_point_interactive(
      aes(
        x, y,
        data_id = Res,
        color = as.factor(Cluster),
        size = Size,
        tooltip = label
      ),
      data = vertices
    ) +
    {
      if (isTRUE(transitionLabels))
        geom_label_repel(
          aes(
            x = (x + xend) / 2,
            y = (y + yend) / 2,
            color = if (transitionProportions == "To") ToClust else FromClust,
            label = round(TransProp, 2)
          ),
          data = edges,
          show.legend = FALSE
        )
      else NULL
    } +
    labs(alpha = "Proportion", color = "Cluster") +
    scale_size(guide = "none", range = c(2, 7)) +
    theme_blank()
}

#' SpatialCPie gadget
#'
#' Runs the SpatialCPie gadget
#' @param counts gene count matrix.
#' @param clusterAssignments list of cluster assignments for each resolution.
#' @param img image to be used as background to the plot (optional).\cr
#' Note: For the tissue image to be correctly aligned to the spatial areas, the
#' pixel.coords argument also needs to be provided.
#' @param pixel.coords `data.frame` with pixel coordinates (optional).
#' The rows should correspond to the columns (spatial areas) in the count file.
#' @param view Shiny gadgets viewer options.
#' Available options: "dialog" (default), "browser", "pane".
#' @param transProp.threshold threshold value for which transitions that are shown in the clustering tree.
#' @return FIXME
#' @keywords arrayplot arraypieplot clustertree
#' @export
#' @import shiny miniUI ggplot2 grid zeallot grid purrr readr
#' @importFrom stats setNames
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
#'   x <- exp(-colSums((centers-x)^2)/50)
#'   x / sum(x)
#' })
#' means <- 100 * profiles %*% mixes
#' counts <- matrix(rpois(prod(dim(means)), means), nrow=nrow(profiles))
#' colnames(counts) <- apply(
#'   coords, 1, function(x) do.call(paste, c(as.list(x), list(sep="x"))))
#' rownames(counts) <- paste("gene", 1:nrow(counts))
#'
#' ## Perform clustering
#' clusterAssignments <- lapply(
#'   2:5, function(x) kmeans(t(counts), centers=x)$cluster)
#'
#' ## Run SpatialCPie
#' runCPie(counts, clusterAssignments)
runCPie <- function(
  counts,
  clusterAssignments,
  img = NULL,
  view = "dialog",
  pixel.coords = NULL
) {
  # Compute coordinates and intersect spots across resolutions
  spots <- clusterAssignments %>%
    map(names) %>%
    reduce(intersect)

  if (!is.null(pixel.coords) != 0) {
    spots <- intersect(spots, rownames(pixel.coords))
    coords <- pixel.coords
  } else {
    c(xcoord, ycoord) %<-% {
       strsplit(spots, "x") %>%
         transpose %>%
         map(as.numeric)
    }
    coords <- as.data.frame(cbind(x = xcoord, y = ycoord))
    rownames(coords) <- spots
  }

  clusterAssignments <- clusterAssignments %>% map(~.[spots])

  # Relabel the data, making sure that all labels in resolution r are in 1:r
  clusterAssignments <- lapply(
    clusterAssignments,
    function(assignments) {
      oldLabels <- unique(assignments)
      newLabels <- setNames(1:length(oldLabels), nm = oldLabels)
      setNames(newLabels[as.character(assignments)], nm = names(assignments))
    }
  )

  # Sort resolutions
  names(clusterAssignments) <- clusterAssignments %>% map_int(max)
  clusterAssignments <- clusterAssignments[
    as.character(sort(as.numeric(names(clusterAssignments))))]

  # Add first resolution (corresponding to the root of the tree)
  if (names(clusterAssignments)[1] != "1") {
    clusterAssignments <- c(
      list("1" = setNames(rep(1, length(spots)), nm = spots)),
      clusterAssignments
    )
  }

  # Relabel the data so as to maximize the number of spots that keep the same
  # label between resolutions
  reassignments <- list(
    unname(head(clusterAssignments, -1)),
    unname(tail(clusterAssignments, -1))
  ) %>%
    transpose %>%
    map(lift(function(prev, cur) {
      setNames(nm = sort(unique(prev))) %>%
        map(
          function(x)
            setNames(nm = sort(unique(cur))) %>%
              map_dbl(function(y) sum(`*`(prev == x, cur == y)))
        ) %>%
        invoke(rbind, .) %>%
        (function(x) {
           nfrom <- nrow(x)
           nto <- ncol(x)
           n <- max(nfrom, nto)
           x_ <- x
           x_ <- do.call(rbind, c(list(x_), rep(list(rep(0, n)), n - nfrom)))
           x_ <- do.call(cbind, c(list(x_), rep(list(rep(0, n)), n - nto)))
           lpSolve::lp.assign(-x_)$solution[1:nfrom, 1:nto, drop = FALSE] %>%
             `colnames<-`(colnames(x)) %>%
             `rownames<-`(rownames(x))
        })
    }))

  clusterAssignments <- c(
    list(clusterAssignments[[1]]),
    list(
      tail(clusterAssignments, -1),
      reassignments
    ) %>%
      transpose %>%
      reduce(
        function(acc, x) {
          c(accReassignments, accAssignments) %<-% acc
          c(assignment, reassignment) %<-% x
          rownames(reassignment) <- rownames(reassignment) %>%
            map(~if (. %in% accReassignments) accReassignments[.] else .)
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
    setNames(names(clusterAssignments))

  # Compute sample-centroid distances
  distances <- clusterAssignments %>%
    map(function(assignments) {
      labels <- unique(assignments)
      centers <- labels %>%
        map(~rowMeans(counts[, which(assignments == .), drop = FALSE])) %>%
        invoke(cbind, .)
      colnames(centers) <- labels
      pairwiseDistance(t(centers), t(counts))
    })

  # Compute colors so that dissimilar clusters are far away in color space
  maxRes <- tail(clusterAssignments, 1)[[1]]
  pca <- sort(unique(maxRes)) %>%
    map(~rowMeans(counts[, maxRes == ., drop = FALSE])) %>%
    invoke(rbind, .) %>%
    prcomp(rank = 2, center = TRUE) %>%
    `$`("x")
  colors <- colorspace::LAB(cbind(
    rep(50, nrow(pca)),
    round((2 * (pca - min(pca)) / (max(pca) - min(pca)) - 1) * 100)
  )) %>%
    colorspace::hex(fixup = TRUE)

  server <- function(input, output, session) {
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
      clusterTree(
        do.call(cbind, clusterAssignments),
        transitionProportions = edgeProportions(),
        transitionLabels = edgeLabels() == "Show",
        transitionThreshold = edgeThreshold()
      ) +
        scale_color_manual(values = colors)
    })

    output$tree <- renderggiraph({
      plot <- ggiraph(
        code = print(treePlot()),
        hover_css = "stroke:#888;stroke-width:0.2em;stroke-opacity:0.5;",
        selected_css = "stroke:#000;stroke-width:0.2em;stroke-opacity:0.5;"
      )

      # Copy selection from the previous tree
      # Note: This could also have been done by sending the `"tree_set"` message
      # on the `"onFlushed"` event. However, that would create a race condition
      # between the message and the browser's loading of the ggiraph dependency
      # file; if the latter is loaded last, the selection would be overwritten.
      # Thus, we instead modify the dependency file directly to include the
      # plot's initial selection.
      if (length(input$tree_selected) > 0) {
        dependency <- paste(
          plot$dependencies[[1]]$src$file,
          plot$dependencies[[1]]$script,
          sep = "/"
        )
        src <- read_file(dependency)

        # Refresh selection annotations when plot is initialized
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

        # Set initial selection
        selection_array <- sprintf("[%s]", lift(paste)(
          paste0("'", input$tree_selected, "'"),
          sep = ", "
        ))
        src <- paste0(
          src,
          sprintf("%s = %s;", plot$x$sel_array_name, selection_array)
        )
        src <- paste0(
          src,
          sprintf("Shiny.setInputValue('tree_selected', %s);", selection_array)
        )

        write_file(src, dependency)
      }

      # Remove UI element for lasso selection etc.
      plot$x$ui_html <- ""

      plot
    })

    ###
    ## ARRAY PLOT
    for (d in tail(names(distances), -1)) {
      # We evaluate the below block in a new frame (with anonymous function
      # call) in order to protect the value of `d`, which will have changed when
      # the reactive expressions are evaluated
      (function() {
        d_ <- d
        infoName <- sprintf("array_info_%s", d_)
        plotName <- sprintf("array_%s", d_)
        assign(envir = parent.frame(), infoName, reactive(
          likeness(distances[[d_]], scoreMultiplier())
        ))
        assign(envir = parent.frame(), plotName, reactive(
          arrayPlot(
            scores = eval(call(infoName)),
            coordinates = coords,
            image =
              if (!is.null(img) &&
                  !is.null(pixel.coords) &&
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
            cat(sprintf("Loading resolution \"%s\"...\n", d_))
            eval(call(plotName))
          },
          width = 400, height = 400
        )
      })()
    }

    output$array <- renderUI({
      lift(div)(
        style = "text-align:center",
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
          ))
      )
    })

    ###
    ## EXPORT
    observeEvent(input$done, {
      stopApp(returnValue = list(
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
      session$sendCustomMessage("tree_set", names(distances))
    })
  }

  if (view == "pane") {
    viewer <- paneViewer()
  } else if (view == "browser") {
    viewer <- browserViewer()
  } else if (view == "dialog") {
    viewer <- dialogViewer("")
  } else {
    viewer <- dialogViewer("")
  }

  runGadget(
    server = server,
    viewer = viewer,
    app = miniPage(
      tags$head(tags$style(HTML(".recalculating { position: relative; z-index: -2 }"))),
      gadgetTitleBar("SpatialCPie"),
      miniContentPanel(
        fillPage(
          sidebarLayout(
            sidebarPanel(width = 3,
              radioButtons("edgeLabels", "Edge labels:", c("Show", "Hide")),
              radioButtons("edgeProportions", "Edge proportions:", c("To", "From")),
              numericInput("edgeThreshold", "Min proportion:", max = 1, min = 0, value = 0.05, step = 0.01)
            ),
            div(style = "text-align: center", mainPanel(ggiraphOutput("tree")))
          ),
          hr(),
          sidebarLayout(
            sidebarPanel(width = 3,
              if (!is.null(img))
                radioButtons("showImage", "HE image:", c("Show", "Hide"))
              else NULL,
              numericInput("simC", "Score multiplier:", max = 10, min = 0.1, value = 1, step = 0.2),
              numericInput("spotOpacity", "Opacity:", max = 100, min = 1, value = 100, step = 10),
              numericInput("spotSize", "Size:", max = 10, min = 1, value = 5, step = 1)
            ),
            mainPanel(uiOutput("array"))
          )
        )
      )
    )
  )
}
