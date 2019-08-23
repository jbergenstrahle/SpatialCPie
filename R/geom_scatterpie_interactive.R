#' @importFrom ggforce StatPie
#' @importFrom ggplot2 layer
#' @importFrom utils modifyList
geom_scatterpie_interactive <- function(
    mapping = NULL,
    data = NULL,
    ...,
    stat = "pie",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
    layer(
        data = data,
        mapping = modifyList(mapping, list(r0 = 0)),
        stat = stat,
        geom = GeomInteractiveScatterpie,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}


#' @importFrom ggiraph GeomInteractivePolygon
#' @importFrom ggplot2 ggproto
GeomInteractiveScatterpie <- ggproto(
    "GeomInteractiveScatterpie",
    GeomInteractivePolygon,
    draw_panel = function(data, panel_scales, coord) {
        g <- GeomInteractivePolygon$draw_panel(data, panel_scales, coord)
        g$id <- match(g$id, unique(g$id))
        g
    }
)
