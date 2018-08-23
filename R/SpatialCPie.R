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
#' @param transitionLabels how to compute and display edge labels. Valid values
#' are:
#' * `NULL`: Don't show any edge labels.
#' * `"To"`: Proportions are computed based on the number of features in the
#' higher resolution cluster.
#' * `"From"`: Proportions are computed based on the number of features in the
#' lower resolution cluster.
#' @return ggplot object of the clustering tree.
#' @keywords internal
#' @import ggplot2 purrr ggiraph ggnetwork intergraph
#' @importFrom igraph graph_from_data_frame layout.reingold.tilford
clusterTree <- function(
  assignments,
  transitionLabels = NULL,
  transitionThreshold = 0.0
) {
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
      trans.prop.from <- trans.count / from.size
      trans.prop.to <- trans.count / to.size

      return(c(trans.count, trans.prop.from, trans.prop.to))
    })

    # Tidy up the results
    trans.df$FromRes <- as.numeric(gsub("res.", "", from.res))
    trans.df$ToRes <- as.numeric(gsub("res.", "", to.res))
    trans.df$TransCount <- trans[1, ]
    trans.df$TransPropFrom <- trans[2, ]
    trans.df$TransPropTo <- trans[3, ]

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
    dplyr::select(Node, everything())

  ## Construct graph
  graph <- transitions %>%
    # Remove edges without any cell...
    filter(TransCount > 0) %>%
    # ...or making up only a small proportion of the new cluster
    filter(TransPropTo > transitionThreshold) %>%
    # Rename the nodes
    mutate(FromNode = paste0("R", FromRes, "C", FromClust)) %>%
    mutate(ToNode = paste0("R", ToRes, "C", ToClust)) %>%
    # Reorder columns
    dplyr::select(FromNode, ToNode, everything()) %>%
    # Build a graph using igraph
    graph_from_data_frame(vertices = nodes)

  ## Construct plot
  data <- ggnetwork(
    graph,
    layout = layout.reingold.tilford(graph),
    cell.jitter = 0.75
  )
  data$label <- sprintf("Resolution %s", data$Res)
  isPoint <- data$x == data$xend & data$y == data$yend
  points <- tail(data[isPoint, ], n = -1)
  edges <- data[!isPoint, ]
  ggplot() +
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = as.factor(Cluster),
      data_id = Res,
      tooltip = label
    ) +
    geom_edges(
      col = "black",
      mapping = aes(alpha = TransPropTo),
      data = edges
    ) +
    geom_point_interactive(
      size = 5,
      data = points,
    ) +
    {
      if (!is.null(transitionLabels) && transitionLabels %in% c("To", "From"))
        geom_edgetext_repel(
          mapping = aes(label = round(
            if (transitionLabels == "To") TransPropTo
            else TransPropFrom,
            2
          )),
          data = edges,
          show.legend = FALSE
      )
      else NULL
    } +
      theme_blank() +
      labs(alpha = "Proportion", color = "Cluster")
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
#' @examples
#' options(device.ask.default = FALSE)
#'
#' ## set up coordinate system
#' coords <- as.matrix(expand.grid(1:10, 1:10))
#'
#' ## generate data set with three distinct genes
#' profiles <- diag(rep(1, 3)) + runif(9)
#' centers <- cbind(c(5, 2), c(2, 8), c(8, 2))
#' mixes <- array_branch(coords, 1) %>%
#'   map(~exp(-colSums((centers-.)^2)/50)) %>%
#'   map(~./sum(.)) %>%
#'   invoke(cbind, .)
#' means <- 100 * profiles %*% mixes
#' counts <- matrix(rpois(prod(dim(means)), means), nrow=nrow(profiles))
#' colnames(counts) <- array_branch(coords, 1) %>% map(~lift(paste)(., sep="x"))
#' rownames(counts) <- paste("gene", 1:nrow(counts))
#'
#' ## perform clustering
#' clusterAssignments <- 2:5 %>% map(~kmeans(t(counts), centers=.)$cluster)
#'
#' ## run SpatialCPie
#' runCPie(counts, clusterAssignments)
runCPie <- function(counts, clusterAssignments, img = NULL, view = "dialog",
                pixel.coords = NULL, transProp.threshold=0.02
) {
  ui <- miniPage(
    tags$head(tags$style(HTML(".recalculating { position: relative; z-index: -2 }"))),
    gadgetTitleBar("SpatialCPie"),
    miniTabstripPanel(
      miniTabPanel("Visualize", icon = icon("area-chart"),
                   miniContentPanel(
                     fillPage(
                       uiOutput("array"),
                       ggiraphOutput("tree")
                     )
                   ),
                   miniButtonBlock(
                     radioButtons("showTrans", "Edge labels:", c("None", "From", "To")),
                     if (!is.null(img)) radioButtons("imgButton", "HE image:", c("Show", "Hide"))
                     else NULL,
                     numericInput("simC", "C", max=10, min=0.1, value=1, step=0.2),
                     numericInput("OpacityButton", "Opacity", max=100, min=1, value=100, step=10),
                     numericInput("SpotSizeButton", "Size", max=10, min=1, value=5, step=1)
                   ))
    )
  )

  server <- function(input, output, session) {
    treeDf <- lift(cbind)(clusterAssignments)
    colnames(treeDf) <- apply(treeDf, 2, max)
    treeDf <- treeDf[, sort(colnames(treeDf))]

    # make sure that there's a root node
    if (colnames(treeDf)[1] > 1) {
      treeDf <- cbind(rep(1, nrow(treeDf)), treeDf)
      colnames(treeDf) <- c(1, colnames(treeDf)[2:ncol(treeDf)])
    }

    # set coordinates of pies
    if (!is.null(pixel.coords) && length(img) != 0) {
      grobHE <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"),
                           interpolate = TRUE)

      coord_df <- pixel.coords
      colnames(coord_df) <- c("x", "y")
    } else {
      grobHE <- NULL

      spots <- clusterAssignments %>% map(names) %>% reduce(intersect)
      c(xcoord, ycoord) %<-% {
         strsplit(spots, "x") %>% transpose %>% map(as.numeric) }
      coord_df <- as.data.frame(cbind(x = xcoord, y = ycoord))
      rownames(coord_df) <- spots
    }

    # compute sample-centroid distances
    distances <- apply(treeDf[, 2:ncol(treeDf)], 2, function(assignments) {
      labels <- unique(assignments)
      centers <- labels %>%
        map(~rowMeans(counts[, which(assignments == .), drop = FALSE])) %>%
        invoke(cbind, .)
      colnames(centers) <- labels
      pairwiseDistance(t(centers), t(counts))
    })

    arrayImgInput <- reactive({
      !is.null(img) &&
        !is.null(pixel.coords) &&
        input$imgButton == "Show"
    })

    spotOpacity <- reactive({
      input$OpacityButton
    })

    # -- set up cluster tree output
    tree_plot <- reactive({
      clusterTree(treeDf, input$showTrans, transProp.threshold)
    })

    init_tree <- TRUE
    output$tree <- renderggiraph({
      plot <- ggiraph(
        code = print(tree_plot()),
        hover_css = "stroke:#888;stroke-width:0.2em;stroke-opacity:0.5;",
        selected_css = "stroke:#000;stroke-width:0.2em;stroke-opacity:0.5;"
      )

      # Get previous selection
      if (isTRUE(init_tree)) {
        selection <- names(distances)
        init_tree <<- FALSE
      } else {
        selection <- input$tree_selected
      }

      # Set selection
      # Note: This could also have been done by sending the `'tree_set'` message
      # on the `'onFlushed'` event. However, that would create a race condition
      # between the message and the browser's loading of the ggiraph dependency
      # file; if the latter is loaded last, the selection would be overwritten.
      # Thus, we instead modify the dependency file directly to include the
      # plot's initial selection.
      if (length(selection) > 0) {
        dep_filename <- paste(
          plot$dependencies[[1]]$src$file,
          plot$dependencies[[1]]$script,
          sep="/"
        )
        dep_src <- read_file(dep_filename)

        # Refresh selection annotations when plot is initialized
        dep_src <- sub(
          "(function init_prop_[^{]*\\{[^}]*)",
          sprintf(
            "\\1refresh_selected('%s', '%s', '%s');",
            plot$x$sel_array_name,
            plot$x$selected_class,
            plot$x$uid
          ),
          dep_src
        )

        # Set initial selection
        selection_array <- sprintf("[%s]", lift(paste)(
          paste0("'", selection, "'"),
          sep = ", "
        ))
        dep_src <- paste0(
          dep_src,
          sprintf("%s = %s;", plot$x$sel_array_name, selection_array)
        )
        dep_src <- paste0(
          dep_src,
          sprintf("Shiny.setInputValue('tree_selected', %s);", selection_array)
        )

        write_file(dep_src, dep_filename)
      }

      # Remove UI element for lasso selection etc.
      plot$x$ui_html <- ""

      plot
    })

    # -- set up array plot output
    for (d in names(distances)) {
      # We evaluate the below block in a new frame (by calling an anonymous
      # function) in order to protect the value of `d`, which will have changed
      # when the reactive expressions are evaluated
      (function() {
        d_ <- d
        info_name <- sprintf("array_info_%s", d_)
        plot_name <- sprintf("array_%s", d_)
        assign(envir = parent.frame(), info_name, reactive(
          likeness(distances[[d_]], input$simC)
        ))
        assign(envir = parent.frame(), plot_name, reactive(
          arrayPlot(
            scores = eval(call(info_name)),
            coordinates = coord_df,
            image = if (arrayImgInput()) grobHE else NULL,
            spotScale = input$SpotSizeButton / 5,
            spotOpacity = spotOpacity() / 100
          ) +
            ggtitle(sprintf("Resolution %s", d_))
        ))
        output[[paste0("plot", d_)]] <- renderPlot(
          {
            cat(sprintf("Loading resolution \"%s\"...\n", d_))
            eval(call(plot_name))
          },
          width = 600, height = 600
        )
      })()
    }

    output$array <- renderUI({
      lift(div)(
        style = "text-align:center",
        sort(input$tree_selected) %>%
          map(~paste0("plot", .)) %>%
          keep(~. %in% names(outputOptions(output))) %>%
          map(~div(
            style = "display: inline-block;",
            div(
              style = paste(
                "position: relative;",
                "width: 600px;",
                "height: 600px;"
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

    observeEvent(input$done, {
      stopApp(returnValue = list(
        tree = tree_plot(),
        piePlots = lapply(
          setNames(input$tree_selected, input$tree_selected),
          function(x) eval(call(sprintf("array_%s", x)))
        ),
        piePlotsInfo = lapply(
          setNames(input$tree_selected, input$tree_selected),
          function(x) eval(call(sprintf("array_info_%s", x)))
        )
      ))
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

  runGadget(ui, server, viewer = viewer)
}
